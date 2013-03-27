/*
Revision of last commit: $Rev: 244 $
Author of last commit: $Author: werner $
Date of last commit: $Date: 2007-10-22 17:44:25 +0900 (Mon, 22 Oct 2007) $
*/

#include "../ModelTuningParameters.h"

ModelTuningParameters::ModelTuningParameters() : 
	tuningParameters(), bounds(), errorValueIsValid(false) {}

ModelTuningParameters::ModelTuningParameters(const int newTParamsLength) 
	: tuningParameters(newTParamsLength), bounds(), errorValueIsValid(false) {
}

ModelTuningParameters::ModelTuningParameters(const vector< double > newTParams, const int newTParamsLength, const vector < double > newBounds) : 
	tuningParameters(newTParamsLength), bounds(), errorValueIsValid(false) {

	ModelTuningParameters::setTuningParameters(newTParams);
	ModelTuningParameters::setBounds(newBounds);
}

ModelTuningParameters::ModelTuningParameters(const string paramString, const int newTParamsLength, const string newBounds) : 
	tuningParameters(newTParamsLength), bounds(), errorValueIsValid(false) {

	ModelTuningParameters::setTuningParameters(paramString, newTParamsLength);
	ModelTuningParameters::setBounds(newBounds, 2*newTParamsLength); 
}

ModelTuningParameters::ModelTuningParameters(const int newTParamsLength, const string newBounds) : 
	tuningParameters(newTParamsLength), bounds(), errorValueIsValid(false) {

	ModelTuningParameters::setBounds(newBounds, 2*newTParamsLength); 
}

double ModelTuningParameters::getLowerBound(const int subscript) const {
	if (subscript < 0 || subscript >= (int)tuningParameters.size()) {
		crash("ModelTuningParameters","Lower bound subscript out of range: "+subscript);
	}
	//Check sanity of bound
	if (bounds[2*subscript] - bounds[2*subscript+1] > __DBL_EPSILON__) crash("ModelTuningParameters","Lower bound " + str(bounds[2*subscript]) + " is larger than upper bound " + str(bounds[2*subscript+1]));
	return bounds[2*subscript];
}

double ModelTuningParameters::getUpperBound(const int subscript) const {
	if (subscript < 0 || subscript >= (int)tuningParameters.size()) {
		crash("ModelTuningParameters","Upper bound subscript out of range: "+subscript);
	}
	//Check sanity of bound
	if (bounds[2*subscript] - bounds[2*subscript+1] > __DBL_EPSILON__) crash("ModelTuningParameters","Lower bound " + str(bounds[2*subscript]) + " is larger than upper bound " + str(bounds[2*subscript+1]));
	return bounds[2*subscript+1];	
}

int ModelTuningParameters::getLength() const {
	return tuningParameters.size();	
}


void ModelTuningParameters::setTuningParameters(const vector< double > newTParams) {

	tuningParameters = newTParams;	

}

void ModelTuningParameters::setTuningParameters(const string paramString, const int newTParamsLength) {
	vector< double > newTParams(newTParamsLength);
	
	// "Parse" the string
	if (paramString != "") {
		istringstream stream(paramString);	
		for (int i = 0; i < newTParamsLength; i++) {
			if (!stream.good()) crash("ModelTuningParameters","Error while converting string into parameters: "+paramString);
			stream >> newTParams[i];
		}
	}

	ModelTuningParameters::setTuningParameters(newTParams);

}

void ModelTuningParameters::setBounds(const vector< double > newBounds) {

	bounds = newBounds;

}

void ModelTuningParameters::setBounds(const string boundString, const int newBoundsLength) {
	if (boundString != "") {
		vector< double > newBounds(newBoundsLength);
		
		// "Parse" the string
    	istringstream stream (boundString);
    	for (int i = 0; i < newBoundsLength; i++) {
			if (!stream.good()) crash("ModelTuningParameters","Error while converting string into bounds: "+boundString);
        	stream >> newBounds[i];
    	}   
    
    	ModelTuningParameters::setBounds(newBounds);
	}
	else {
		ModelTuningParameters::setBounds(vector< double >(0));	
	}
}

double &ModelTuningParameters::operator[]( int subscript ) {
	if (subscript < 0 || subscript >= (int)tuningParameters.size()) {crash("ModelTuningParameters","Subscript out of range: "+subscript);}
	return tuningParameters[subscript];
}

const double &ModelTuningParameters::operator[]( int subscript ) const {
	if (subscript < 0 || subscript >= (int)tuningParameters.size()) {crash("ModelTuningParameters","Subscript out of range: "+subscript);}
	return tuningParameters[subscript];
}

void ModelTuningParameters::setErrorValue(const double newValue) {
	errorValue = newValue;
	errorValueIsValid = true;
}

void ModelTuningParameters::resetErrorValue() {
	errorValueIsValid = false;
}

double ModelTuningParameters::getErrorValue() const {

	if (!errorValueIsValid) {
		crash("ModelTuningParameters","Getting error value which is uninitialized");
	}
	return errorValue;

}

string ModelTuningParameters::toString() const {
	ostringstream o;
	o << "{ ";
	for (int i = 0; i < getLength(); i++) {
   		o << (*this)[i] << " ";		
	}
	o << "}";
	return o.str();
}


void ModelTuningParameters::printOn(OutputChannel & output) const {

	int length = tuningParameters.size();
	output << length; 
	for (int i = 0; i < length; i++) {	
		output << tuningParameters[i];
	}
	int boundsLength = bounds.size();
	output << boundsLength; 
	for (int i = 0; i < boundsLength; i++) {	
		output << bounds[i];
	}
	output << (int)errorValueIsValid;
	output << errorValue;

#ifdef WITH_MOEA
	output << (int) errorValuesVector.size();
	for (unsigned i =0; i < errorValuesVector.size(); i++)
		output << errorValuesVector.at(i);
#endif

}

void ModelTuningParameters::readFrom(InputChannel & input) {

	int length;
	input >> length;
	tuningParameters = vector< double >(length);
	for (int i = 0; i < length; i++) {
		input >> tuningParameters[i];
	}
	
	int boundsLength;
	input >> boundsLength;
	bounds = vector< double >(boundsLength);
	for (int i = 0; i < boundsLength; i++) {
		input >> bounds[i];
	}
	
	int fValid;
	input >> fValid; errorValueIsValid = (bool)fValid;
	input >> errorValue;

#ifdef WITH_MOEA
	input >> fValid; // Size
	if (fValid>0)
	{	
		double dummy;
		this->errorValuesVector.resize(fValid);
		for (int i =0; i < fValid; i++)
		{
			input >> dummy;
			this->errorValuesVector.push_back(dummy);
		}
	}
#endif

}




#ifdef WITH_XMLRPC
void ModelTuningParameters::readFromXmlrpc(xmlrpc_c::value_array& array)
{
	int p = 0;
	int total = paramList.getInt(p++); // Total run parameters

	vector<double > modelParameters(total);
	for (int i =0; i < total; i++)
	{
		modelParameters.push_back(paramsList.getDouble(p++));
	}
	this->setTuningParameters(modelParameters);

	/* Bounds */
	total = paramsList.getInt(p++); // Total bounds
	vector<double> bounds(total);
	for (int i=0; i < total; i++)
	{
		bounds.push_back(paramsList.getDouble(p++));
	}
	this->setBounds(bounds);

	total = paramsList.getInt(p++); // fValid;
	this->setErrorvalue(paramsList.getDouble(p++)); // Error value

#ifdef WITH_MOEA
	total = paramsList.getInt(p++); // Read moea objectives

	modelParameters.clear();
	modelParameters.resize(total);
	for (int i=0; i < total; i++)
		modelParameters.push_back(paramList.getDouble(p++));
	this->setErrorValuesVector(modelParameters);
#endif 

}


void ModelTuningParameters::writeToXmlrpc(xmlrpc_c::value_array& array)
{
	vector<xmlrpc_c::value> arr;
	int length = tuningParameters.size();
	arr.push_back(xmlrpc_c::value_int(length));
	
	for (int i = 0; i < length; i++) {	
		arr.push_back(xmlrpc_c::value_double(tuningParameters[i]));
	}

	int boundsLength = bounds.size();
	arr.push_back(xmlrpc_c::value_int(boundsLength));
	for (int i = 0; i < boundsLength; i++) {	
		arr.push_back(xmlrpc_c::value_double(bounds[i]));
	}
	arr.push_back(xmlrpc_c::value_int((int)errorValueIsValid));

	arr.push_back(xmlrpc_c::value_double(errorValue));

#ifdef WITH_MOEA
	arr.push_back(xmlrpc_c::value_int((int)errorValuesVector.size()));
	for (unsigned i =0; i < errorValuesVector.size(); i++)
		arr.push_back(xmlrpc_c::value_double((double)errorValuesVector.at(i)));
#endif

	
}
#endif


