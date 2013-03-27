/*
Revision of last commit: $Rev$
Author of last commit: $Author$
Date of last commit: $Date$
*/

#include "../MOEAErrorValueCalculator.h"
#include "../MOEAFitterInterface.h"
#include "../DataTrace.h"
#include <algorithm>
#include <math.h>

inline static double dabs(double c)
{
	if (c<0)return -c;
	return c;
}


/* Should normally never be needed */
inline static double d_sqrt(double diff)
{
	if (diff >0)
		return sqrt(diff);
	else
		return 1.0;
}



static double getDistance(const DataTrace& a, const DataTrace& b)
{
	if (b.getValidLength() <  a.getValidLength())
		throw std::runtime_error("Traces differ in number of points: " + str(b.getLength()) + "<"+str(a.getLength()));

	double sum_i =0;
	for (int i =0 ;i < a.getValidLength(); i++)
	{
		double diff = (b.get(i) - a.get(i))/2.0;
		sum_i += 1.0 / (1.0+diff*diff);
	}
	sum_i /= (a.getLength()/a.getSamplingFrequency());
	return  sum_i;
}





/* Wrapper for the Phase plane error function */
struct matrix_error_function : public MOEAErrorValueCalculator:error_function
{

	int totalExpTraces =0;
	vector< VdVdtMatrix* > expVdVdtMatrices;
	VdVdtMatrix* 	modelVdVdtMatrix;


	void readExperimentalTrace(DataTrace& trace)
	{
		if (fixedParams["VdVdtMatrixType"] == "Distance")
			expVdVdtMatrices[totalExpTraces] = new DistanceVdVdtMatrix(expData[nTrace], FixedParameters(fixedParams["VdVdtMatrixParameters"],fixedParams.getGlobals()));
		
		else
			expVdVdtMatrices[totalExpTraces] = new DirectVdVdtMatrix(expData[nTrace], FixedParameters(fixedParams["VdVdtMatrixParameters"],fixedParams.getGlobals()));

		cout << "Exp matrix " << totalExpTraces<<": \n" << expVdVdtMatrices[nTrace]->toString() << endl;
		totalExpTraces++;
	}


	void readModelTrace(DataTrace& tr)
	{
		modelVdVdtMatrix->readFrom(tr);
	}

	virtual vector<float> getErrorValues();

	void initialize(ModelResults& expData, FixedParameters& params)
	{

		this->expVdVdtMatrices = vector< VdVdtMatrix *> (expData.getLength(), (VdVdtMatrix*)NULL);

		if (params["VdVdtMatrixType"] == "Distance")
			modelVdVdtMatrix = new DistanceVdVdtMatrix(FixedParameters(fixedParams["VdVdtMatrixParameters"],fixedParams.getGlobals()));
		else
			modelVdVdtMatrix = new DirectVdVdtMatrix(FixedParameters(fixedParams["VdVdtMatrixParameters"],fixedParams.getGlobals()));
	}


	virtual ~matrix_error_function() {
		delete modelVdVdtMatrix;
		for (int i = 0; i < (int)expVdVdtMatrices.size(); i++) {
			delete expVdVdtMatrices[i];
		}
	}

};


struct distance_error_function : public MOEAErrorValueCalculator::error_function 
{
	void initialize(ModelResults& expData, FixedParameters& params)
	{
		
	}

	void readModelTrace();


	virtual ~distance_error_function() {}
};



/* Collectively, returns spike-feature-related errors */
struct spike_feature_error_functions: public MOEAErrorValueCalculator::error_function
{
	double spikeThreshold;	/* Voltage threshold for detecting spikes */
	double spikeTimeWindow;	/* Time to scan for overshoot + hyperpolarization depth values after spike onset */
	int compareAverages;	/* How to calculate objective values: if 1 => calculate the difference in averages between experimental and model traces - if 0 => calculate the differences for each spike */

	vector< feature_values> expSpikeFeatureValues;

	void scanForSpikes(DataTrace& dtr, vector<spike_struct>& spikes);

	void initialize(ModelResults& expData, FixedParameters& params)
	{
		this->spikeThreshold  = toDouble(params["spikeThreshold"]);
		this->spikeTimeWindow = toDouble(params["spikeTimeWindow"]);
		this->compareAverages = toInt(params["CompareAverages"]);
		this->expSpikes.push_back(vector<spike_struct>()); 

		vector<spike_struct>& spikes = this->expSpikes.back();

		for (int nTrace =0; nTrace < expData.getLength(); nTrace++)
		{
			DataTrace& dtr = expData[nTrace];
			this->expSpikes.push_back(vector<spike_struct>()); 
			vector<spike_struct>& spikes = this->expSpikes.back();

			this->expSpikeFeatureValues.push_back(feature_values()); 
			feature_values& fv = this->expSpikeFeatureValues.back();

			this->scanForSpikes(dtr, spikes);
			fv.readSpikes(spikes);
			fv.finalize(dtr.getStopTime() - dtr.getStartTime());
		}
	}


	vector<double> getErrorValues(ModelResults&  modelResult)
	{
		
		for (unsigned nTrace =0; nTrace < modeResult.getLength(); nTrace++)
		{
			DataTrace& modelTrace = results[i][nTrace];
			vector<spike_struct > spikes;

			this->scanForSpikes(modelTrace, spikes);
			cerr << " model trace #"<<i<< " has "<< spikes.size() << " spikes." << endl;

			if (this->compareAverages)
			{
				/* Calculate spike objectives by averaging (recommended method) */
				feature_values mResults;
				mResults.readSpikes(spikes);
				mResults.finalize(modelTrace.getStopTime() - modelTrace.getStartTime());

				feature_values& expResult = this->expSpikeFeatureValues.at(nTrace);

				/* In the FINS paper, the suggestion is to return errors divided by the experimental stdev. However since we only
				 * (usually) have 1 (?) experimental trace, we will use absolute values */
				allObjectives[MEASURE_WIDTH] += results[i][nTrace].getWeight() * (dabs(mResults.width - expResult.width)); // / expResult.widthSd);
				allObjectives[MEASURE_APO] += results[i][nTrace].getWeight() * (dabs(mResults.apo - expResult.apo)); // / expResult.apoSd);
				allObjectives[MEASURE_AHD] += results[i][nTrace].getWeight() * (dabs(mResults.ahd - expResult.ahd)); // / expResult.ahdSd);
				allObjectives[MEASURE_SRATE] += results[i][nTrace].getWeight() * (dabs(mResults.sRate - expResult.sRate)); // / expResult.sRateSd);
				allObjectives[MEASURE_AINDEX] += results[i][nTrace].getWeight() * (dabs(mResults.aIndex - expResult.aIndex)); // / expResult.aIndexSd);
			}
			else
			{
				/* Calculate spike differences comparing spikes one-by-one 
				 * This will lead to unreasonable errors when the number of spikes is not the same
				 */
				vector<spike_struct>& expSpikesVector = this->expSpikes.at(nTrace);
				spike_struct dummy;
				dummy.width = 100.0;

				for (unsigned k=0; k < expSpikesVector.size(); k++)
				{
					spike_struct& mySpike = (k >=spikes.size()) ? dummy :   spikes.at(k);
					spike_struct& expSpike = expSpikesVector.at(k);

					allObjectives[MEASURE_WIDTH] += results[i][nTrace].getWeight() * (dabs(mySpike.width - expSpike.width));
					allObjectives[MEASURE_APO] += results[i][nTrace].getWeight() * (dabs(mySpike.overshoot - expSpike.overshoot));
					allObjectives[MEASURE_AHD] += results[i][nTrace].getWeight() * (dabs(mySpike.ahd - expSpike.ahd));
					allObjectives[MEASURE_SRATE] += results[i][nTrace].getWeight() * (dabs(mySpike.isi - expSpike.isi));
					allObjectives[MEASURE_AINDEX] += results[i][nTrace].getWeight() * (dabs(mySpike.startTime - expSpike.startTime));
				}
			}
		}
	}

	virtual ~spike_feature_error_functions();
};




/* Accumulate spike information */
void MOEAErrorValueCalculator::feature_values::readSpikes(vector<spike_struct>& spikes) 
{
	bool latencySet=false;

	for (unsigned i=0; i < spikes.size(); i++)
	{
		spike_struct& s = spikes.at(i);
		width_x += s.width;
		width_x2 += s.width*s.width;

		apo_x += s.overshoot;
		apo_x2 += s.overshoot * s.overshoot;
		ahd_x += s.ahd;
		ahd_x2 += s.ahd*s.ahd;

		if (i>0)
		{
			double s_isi1 = s.isi - spikes.at(i-1).isi;
			double s_isi2 = s.isi + spikes.at(i-1).isi;

			if (s_isi2 >0 && s.isi >0.0)
			{
				double fr = s_isi1/s_isi2;
				aIndex_x += fr;
				aIndex_x2 += fr*fr;
				this->totalISI++;

				s_isi1 = 1.0/s.isi;
				sRate_x += s_isi1;
				sRate_x2 += s_isi1*s_isi1;
			}
		}

		if (!latencySet)
		{
		    latencySet = true;
		    latency = s.startTime;
		}

		this->totalSpikes++;
	}
}

/* Once finished reading spike information, call finalize() to 
 * calculate standard deviations */
void MOEAErrorValueCalculator::feature_values::finalize(double totalTime)
{
	if (this->totalSpikes < 1)
		return;
	double totalV = this->totalSpikes;

	if (totalTime>0.0)
		sRate = totalV / totalTime;

	apo = apo_x / totalV;
	apoSd = d_sqrt((apo_x2/totalV) - apo*apo);

	ahd = ahd_x / totalV;
	ahdSd = d_sqrt((ahd_x2/totalV) - ahd*ahd);

	width = width_x / totalV;
	widthSd = d_sqrt((width_x2/totalV) - width*width);

	if (totalISI>0)
	{
		sRate = sRate_x / (double)totalISI;
		sRateSd = d_sqrt((sRate_x2/(double)totalISI) - sRate*sRate);

		aIndex = this->aIndex_x / ((double)this->totalISI);
		aIndexSd = d_sqrt((aIndex_x2/((double)this->totalISI)) - aIndex*aIndex);
	}
	else
	{
		sRateSd = aIndexSd =1.0;/* Just to avoid division/0 */
	}
}





void MOEAErrorValueCalculator::feature_values::printOn(ostream& out)
{
	out << " Width = " << width << " [" << widthSd << "] " << endl;
	out << " APO = " << apo << " [" << apoSd << "] " << endl;
	out << " AHD = " << ahd << " [" << ahdSd << "] " << endl;
	out << " Rate = " << sRate << " [" << sRateSd << "] " << endl;
	out << " Acc.Index = " << aIndex << " [" << aIndexSd << "] " << endl;
	out << " Latency = " << latency << endl;
}





MOEAErrorValueCalculator::MOEAErrorValueCalculator(ModelInterface * interface, ExperimentInterface * experiment, FixedParameters params) 
	: ErrorValueCalculator(interface), FixedParamObject(params)
{
	if (toInt(fixedParams["enableFileExport"]) > 0) {
		this->enableFileExport(fixedParams["exportFile"]);
	}

	cout << "Reading experimental data " << endl;
	ModelResults expData = experiment->getData();	
	if (expData.size() ==0)
		throw std::runtime_error("Ntrace==0! No experimental data found/read, aborting");
	cout << "Read."<< endl;

	istringstream paramstring(params["Objectives"]);
	string calcType;
	double calcWeight;
	cerr <<  "Specified error values: " << params["Objectives"] << endl;
	int ind =0;

	while ((paramstring >> calcType))
	{
		if (!(calcType.length()>1 && ((paramstring>>calcWeight))))
			break;

		error_function* err = 0;


		if (calcType == "Matrix") 
		{
			err = new matrix_error_function();
		}

		else if (calcType == "SpikeFeatures")
		{
			err = new spike_feature_error_functions();
		}
		else if ( calcType == "Distance")
		{
			err = new distance_error_function();
		}
		else 
		{ 
			throw new std::runtime_error("No matching Calculator type: " + calcType);
		}


		err->initialize(expData, fixedParams);
		this->errorFunctions->push_back(err);
	}

}



MOEAErrorValueCalculator::~MOEAErrorValueCalculator()
{
	exportFileStream.close();
	delete modelVdVdtMatrix;
	for (int i = 0; i < (int)expVdVdtMatrices.size(); i++) {
		delete expVdVdtMatrices[i];
	}
}


void MOEAErrorValueCalculator::calculateErrorValue(ModelTuningParameters & params) {

	vector< ModelTuningParameters > paramList(1);
	paramList[0] = params;
	
	calculateParallelErrorValue(paramList);

	/// This is necessary because otherwise the value is not transferred 
	/// since no reference is passed to calculateParallelErrorValue
	params.setErrorValue(paramList[0].getErrorValue());
	params.setErrorValuesVector(paramList[0].getErrorValuesVector());

}


void MOEAErrorValueCalculator::calculateParallelErrorValue(vector<ModelTuningParameters> & paramList)
{
	vector< ModelResults > results = model->runParallelModel(paramList);

	vector< double > errorValues(paramList.size());

	for (unsigned int i = 0; i < paramList.size(); i++) 
	{
		
		double allObjectives[MEASURE_MAX] = {0};


		if (!results.at(i).getLength())
		{
			/* FIXME: this will happen if NEURON encounters division by zero errors. We have adapted the code
			 * so that neurofitter won't crash when it runs unattended, but we must also disregard those values.
			 */
			cerr << "!!It seems like we got bad data in model output at run #" << i << " setting to highest error values" << endl;
			paramList[i].setErrorValue(999999999.0);

			vector<double>& ov = paramList.at(i).getErrorValuesVector();
			for (int k=0; k < MEASURE_MAX; k++)
				ov.push_back(999999999.9);
			continue;
		}



		for (unsigned ob=0; ob < this->objFunctions.size(); i++)
		{
			this->objFunctions[i]->readTrace(modelTrace);
		}


		/* Evaluate all traces */
		for (int nTrace = 0; nTrace < results[i].getLength(); nTrace++) 
		{
			DataTrace& modelTrace = results[i][nTrace];

			for (unsigned ob=0; ob < this->objFunctions.size(); i++)
			{
				this->objFunctions[i]->readTrace(modelTrace);
			}



			vector<spike_struct > spikes;

			this->scanForSpikes(modelTrace, spikes);
			cerr << " model trace #"<<i<< " has "<< spikes.size() << " spikes." << endl;

			// Calculate phase plane error
			if (this->objIndex[MEASURE_MATRIX] != -1)
			{
				modelVdVdtMatrix->readFrom(results[i][nTrace]);
				allObjectives[MEASURE_MATRIX] += results[i][nTrace].getWeight() * expVdVdtMatrices.at(nTrace)->compare(*modelVdVdtMatrix);
				//showMessage(modelVdVdtMatrix->toString() + "\n",5,fixedParams);        	
			}

			if (this->objIndex[MEASURE_EUCLIDEAN] != -1)
			{
				allObjectives[MEASURE_EUCLIDEAN] += results[i][nTrace].getWeight() * getDistance(this->expTraces.at(nTrace), modelTrace);
			}

			if (this->compareAverages)
			{
				/* Calculate spike objectives by averaging (recommended method) */
				feature_values mResults;
				mResults.readSpikes(spikes);
				mResults.finalize(modelTrace.getStopTime() - modelTrace.getStartTime());

				feature_values& expResult = this->expSpikeFeatureValues.at(nTrace);

				/* In the FINS paper, the suggestion is to return errors divided by the experimental stdev. However since we only
				 * (usually) have 1 (?) experimental trace, we will use absolute values */
				allObjectives[MEASURE_WIDTH] += results[i][nTrace].getWeight() * (dabs(mResults.width - expResult.width)); // / expResult.widthSd);
				allObjectives[MEASURE_APO] += results[i][nTrace].getWeight() * (dabs(mResults.apo - expResult.apo)); // / expResult.apoSd);
				allObjectives[MEASURE_AHD] += results[i][nTrace].getWeight() * (dabs(mResults.ahd - expResult.ahd)); // / expResult.ahdSd);
				allObjectives[MEASURE_SRATE] += results[i][nTrace].getWeight() * (dabs(mResults.sRate - expResult.sRate)); // / expResult.sRateSd);
				allObjectives[MEASURE_AINDEX] += results[i][nTrace].getWeight() * (dabs(mResults.aIndex - expResult.aIndex)); // / expResult.aIndexSd);
			}
			else
			{
				/* Calculate spike differences comparing spikes one-by-one 
				 * This will lead to unreasonable errors when the number of spikes is not the same
				 */
				vector<spike_struct>& expSpikesVector = this->expSpikes.at(nTrace);
				spike_struct dummy;
				dummy.width = 100.0;

				for (unsigned k=0; k < expSpikesVector.size(); k++)
				{
					spike_struct& mySpike = (k >=spikes.size()) ? dummy :   spikes.at(k);
					spike_struct& expSpike = expSpikesVector.at(k);

					allObjectives[MEASURE_WIDTH] += results[i][nTrace].getWeight() * (dabs(mySpike.width - expSpike.width));
					allObjectives[MEASURE_APO] += results[i][nTrace].getWeight() * (dabs(mySpike.overshoot - expSpike.overshoot));
					allObjectives[MEASURE_AHD] += results[i][nTrace].getWeight() * (dabs(mySpike.ahd - expSpike.ahd));
					allObjectives[MEASURE_SRATE] += results[i][nTrace].getWeight() * (dabs(mySpike.isi - expSpike.isi));
					allObjectives[MEASURE_AINDEX] += results[i][nTrace].getWeight() * (dabs(mySpike.startTime - expSpike.startTime));

				}
			}

			double diff =0; 

			for (int k=0; k < MEASURE_MAX;k++)
			{
				if (this->objIndex[k]>=0)
					diff += allWeights[k] * allObjectives[k];
			}

			errorValues[i] += diff;

		}

		numberOfEvaluations++;

		errorValueHistory.push_back(pair< ModelTuningParameters, double >(paramList[i], errorValues[i]));

		//showMessage("Eval: " + str(numberOfEvaluations) + " Generation: " + str(numberOfGenerations) + " Calculated ErrorValue of: " + paramList[i].toString() + ": " + str(errorValues[i]) + "\n",3,fixedParams);

		if (exportFileStream.is_open()) 
		{
			exportFileStream << numberOfGenerations << " "<< numberOfEvaluations ;
			cerr << numberOfGenerations << ":"<< numberOfEvaluations << ":{";
			for (int j = 0; j < paramList[i].getLength(); j++) 
			{
				exportFileStream <<" " << (paramList[i][j]);
				cerr << " " << (paramList[i][j]);
			}

			cerr << "} => {";
			for (int j =0; j < MEASURE_MAX; j++)
			{
				if (this->objIndex[j]>=0)
				{
					exportFileStream << " " << allObjectives[this->objIndex[j]];
					cerr << " " << allObjectives[this->objIndex[j]];
				}
			}

			exportFileStream << " " << errorValues[i];
			exportFileStream << endl;
		}

		cerr << "} sum=" << errorValues[i] <<endl;

		/* Return the objective values back in the correct order (as specified in the config file) */
		vector<double>& objectivesValues = paramList.at(i).getErrorValuesVector();
		objectivesValues.clear();

		for (int j =0; j < MEASURE_MAX; j++)
		{
			if (this->objIndex[j]>=0)
			{
				objectivesValues.push_back( allObjectives[this->objIndex[j]]);
			}
		}

		/* Sets the weighted error for non-moea algorithms */
		paramList[i].setErrorValue(errorValues[i]);
	}
	numberOfGenerations++;
}


vector< pair< ModelTuningParameters, double > > MOEAErrorValueCalculator::getErrorValueHistory() {
	return errorValueHistory;
}


void MOEAErrorValueCalculator::enableFileExport(const string fileName) {
	exportFileStream.open(fileName.c_str(), ios::out);
	
	showMessage("MOEAErrorValueCalculator: Enabled export to file: " + fileName + "\n",3,fixedParams);        	
}

void MOEAErrorValueCalculator::disableFileExport() {
	exportFileStream.close();
}
   
string MOEAErrorValueCalculator::getExportFile() const {
	return exportFile;
}




void MOEAErrorValueCalculator::scanForSpikes(DataTrace& dtr, vector<spike_struct>& spikes)
{
	spike_struct spike;

	if (dtr.getLength() < 3)
		return;
	
	double timePassed = dtr.getStartTime();
	//double timePassed =0.0;
	double timeStep = 1.0 / dtr.getSamplingFrequency();
	bool inSpike = false;
	double lastV = dtr.get(0);
	//cerr <<" Trace start time = "<< timePassed<< " length = "<< dtr.getLength()<< endl;
	int spikeStartIndex=0;
	for (int i=1; i < dtr.getValidLength(); i++)
	{
		double thisV = dtr.get(i);
		if (thisV > this->spikeThreshold && lastV <= this->spikeThreshold)
		{
			/* Start of a new spike, store the previous one (if any) */
			//cerr<<"T="<<timePassed <<"V = " << thisV <<"lastV="<<lastV<< endl;
			if (inSpike)
			{

				spike.isi = timePassed - (spike.startTime + spike.width);

				//cerr <<"Spike  #"<< str((int)spikes.size()) <<" t=" <<str(spike.startTime) <<" ov="<< str(spike.overshoot) << " ahd=" << str(spike.ahd) << " w=" << str(spike.width) << " isi="<< spike.isi <<  endl;

				spikes.push_back(spike);
				spike.reset();
			}
		 	inSpike = true;
			spike.startTime = timePassed;
			spikeStartIndex =i;
		}


		/* constrain ahd detection in a time range after onset */
		if (inSpike && timePassed < spike.startTime + this->spikeTimeWindow)
		{
			if (spike.ahd > thisV)
				spike.ahd = thisV;
			if (spike.overshoot < thisV)
				spike.overshoot = thisV;
		}


		if (inSpike && thisV < this->spikeThreshold && lastV >= this->spikeThreshold)
		{
			/* End of spike. Determine spike width at half maximum */ 
			double halfMax = this->spikeThreshold + (spike.overshoot - this->spikeThreshold)/2.0;
			bool inHM=false;
			int stStartIndex=0;
			for (int k = spikeStartIndex+1; k < i ; k++)
			{
				if (dtr.get(k-1) < halfMax && dtr.get(k) >=halfMax)
				{
					inHM = true;
					stStartIndex=k;
				}
				else if (dtr.get(k-1) >= halfMax && dtr.get(k) <halfMax)
				{
					if (inHM)
						spike.width = (k-stStartIndex)*timeStep;
					inHM = false;
				}
			}
		}
		lastV = thisV;
		timePassed += timeStep;
	}


	/* Append the last spike */
	if (inSpike)
	{
		spike.isi = timePassed - (spike.startTime + spike.width);
		spikes.push_back(spike);
	}
}



