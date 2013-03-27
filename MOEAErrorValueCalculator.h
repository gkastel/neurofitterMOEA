/*
Revision of last commit: $Rev$
Author of last commit: $Author$
Date of last commit: $Date: 2006-12-11 18:44:13 +0900 (Mon, 11 Dec 2006) $
*/


#ifndef NEUROFITTER_MOEAErrorValueCALCULATOR_H
#define  NEUROFITTER_MOEAErrorValueCALCULATOR_H 1

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#include "ModelInterface.h"
#include "ExperimentInterface.h"
#include "ModelTuningParameters.h"
#include "ErrorValueCalculator.h"
#include "MOEAFitterInterface.h"
#include "FixedParamObject.h"
#include "../DirectVdVdtMatrix.h"
#include "../DistanceVdVdtMatrix.h"



/** 
 * MOEA error calculator - calculate multiple objective (error) values.
 * Based on the reports on on http://frontiersin.org/neuroscience/abstract/10.3389/neuro.01/1.1.001.2007 
 * This error value calculator is used to combine many spike features at once. 
 *
 * When the fitter is MOEA Fitter, the error values will be stored in an array, and used to approach
 * the best pareto front. For all the other Fitters, the value returned will be a weighted sum of the objective 
 * functions (errors) calculated.
 *
	<ErrorValueCalculatorType>MOEA</ErrorValueCalculatorType>
	<ErrorValueCalculatorParameters>
	<Objectives>

<!--    
	Error (objective) functions to use. To omit an error value, delete the line. The weights are only used
	in the case a non-MOOP fitter is selected, and are used to sum up the particular error value. For MOEA, 
	the weights are irrelevant.

-->
<!--    Error Value Function            Weight --> 

	Matrix				1.0 <!-- 'Matrix' is the Phase-plane matrix error value -->

	APWidth 			1.0 <!-- 'APWidth' is the action potential width error, measured from the threshold voltage value (see below) -->

	APOvershoot  			1.0  <!-- 'APOvershoot' is the amplitude-of-spike error -->

	AfterHyperpolarizationDepth 	0.6  <!-- 'AfterHyperpolarizationDepth' is the AHD, measured as the minimum of voltage between two spikes -->

	SpikeRate   			1.0 <!-- 'SpikeRate' is the error in the rate of spiking, measured from the experimental mean -->

	AccomodationIndex  		1.0  <!-- 'AccomodationIndex' is the error in model accomodation index -->

	Distance  			1.0  <!-- 'Distance' is a measure of absolute difference between traces   -->

	</Objectives>

	<CompareAverages>1</compareAverages><!-- if enabled, will calculate differences between (average values of objectives in experimental trace) - (average values in model trace), if disabled the sum of differences in objectives in each spike will be calculated -->

	<enableFileExport>1</enableFileExport><!-- if enabled, file exporting will write to the file 
						the format is:  GenerationNo IndividualNo ModelParam1 ModelParam2 ... Objective1Error Objective2Error ... -->
	<exportFile>neurofitter.log</exportFile>

	<spikeThreshold>-10</spikeThreshold>      <!-- Threshold (in whatever units you use in your data files) to start measuring spikes from.  -->

	<spikeTimeWindow>0.5</spikeTimeWindow>    <!-- Max duration of a spike (in whatever units your data files have) -->

	<!-- If the  matrix measure is used, these are the parameters of the phase-plane analyzer -->
	<VdVdtMatrixType>Direct</VdVdtMatrixType> 
	<VdVdtMatrixParameters>
		<vLength>30</vLength>
		<dVdtLength>30</dVdtLength>
		<minimalV>-100</minimalV>
		<maximalV>100</maximalV>
		<comparePrecision>1e-15</comparePrecision>
		<numericOutputFormat>0</numericOutputFormat>
		<SumOfSquareRoots>1</SumOfSquareRoots>
	</VdVdtMatrixParameters>
	</ErrorValueCalculatorParameters>

 *
 */
class MOEAErrorValueCalculator : public ErrorValueCalculator, public FixedParamObject {

public:

	/* A single spike features
	 */
	struct spike_struct
	{
		float startTime, width, overshoot, ahd, isi;
		void reset()
		{
			startTime = width = overshoot = ahd = isi = 0;
		}

		spike_struct() {reset(); }

	};

	/* One for each parameter set */
	struct feature_values {

		int totalSpikes, totalISI; /* Total number of spikes encountered, total number if interspike intervals */
		double aIndex, aIndexSd, aIndex_x, aIndex_x2; /* Accomodation index - mean, variance, sum_i(x), sum_i(x^2) */
		double sRate, sRateSd, sRate_x, sRate_x2;     /* Spike rate */
		double apo, apoSd, apo_x, apo_x2; 	      /* Action potential overshoot (amplitude) */
		double ahd, ahdSd, ahd_x, ahd_x2;             /* After hyperpolarization depth */
		double width, widthSd, width_x, width_x2;     /* Spike width, counted from onset threshold */
		double latency, latencySd;		      /* First spike latency - not used */
		bool latencySet;

		feature_values()
		{
			reset();
		}

		void reset()
		{
			memset(this, 0, sizeof(*this));
		}

		/* Accumulate data from a spikes vector */
		void readSpikes(vector<spike_struct>& spikes);
		/* Calculate SDs and rates */
		void finalize(double totalTime);

		void printOn(ostream& stream);
	};



	MOEAErrorValueCalculator(ModelInterface * interface, ExperimentInterface * experiment, FixedParameters); 
	~MOEAErrorValueCalculator();

	virtual void calculateErrorValue(ModelTuningParameters & param);

	virtual void calculateParallelErrorValue(vector< ModelTuningParameters > & params);

	virtual void enableFileExport(const string fileName);
	virtual void disableFileExport();

	virtual string getExportFile() const;
	virtual vector< pair< ModelTuningParameters, double > > getErrorValueHistory();


	struct error_function 
	{
		virtual void readTrace(DataTrace& tr) =0;
		/* Finalize the errors and print them out */
		virtual vector<float>& getErrorValues();
			
		void initialize(FixedParameters& params) =0;
	};


private:

	/* Stores the experimental data traces for later comparison */
	vector<DataTrace> expTraces;

	vector<error_function*> errorFunctions;


	vector<double> errorWeights;  /* Weight of each error function */
	vector<int> errorValueOrder;  /* Used to order the error values in the log files in the same order they were specified in the configuration file */

	ofstream exportFileStream;
};

#endif
