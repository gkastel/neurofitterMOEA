/*
Revision of last commit: $Rev: 222 $
Author of last commit: $Author: werner $
Date of last commit: $Date: 2007-10-04 14:03:41 +0900 (Thu, 04 Oct 2007) $
*/

#include "../XMLRPCErrorValueCalculator.h"
#ifdef WITH_MOEA
#include "../MOEAErrorValueCalculator.h"
#endif



XMLRPCErrorValueCalculator::XMLRPCErrorValueCalculator(ModelInterface* model, ExperimentInterface* experiment, FixedParameters params) 
	: ErrorValueCalculator(NULL), FixedParamObject(params) 
{
	FixedParameters fitFixedParams(fixedParams["ErrorValueCalculatorParameters"],fixedParams.getGlobals());
	FixedParameters paramHosts(fixedParams["XMLRPCServers"],fixedParams.getGlobals());

	if (fixedParams["ErrorValueCalculatorType"] == "Matrix") {
		localErrorValue = new MatrixErrorValueCalculator(model,experiment,fitFixedParams);
	}
#ifdef WITH_MOEA
	else if (fixedParams["ErrorValueCalculatorType"] == "MOEA")
		localErrorValue = new MOEAErrorValueCalculator(model,experiment,fitFixedParams);
#endif
	else throw std::runtime_error("No matching error value calculator type");
	
	if (rank != 0) startSlave();
	
	if (rank == 0) {                                
		if (toInt(fixedParams["enableFileExport"]) > 0) {
			this->enableFileExport(fixedParams["exportFile"]);
		}
	}
}

XMLRPCErrorValueCalculator::~XMLRPCErrorValueCalculator() {
	exportFileStream.close();
	if (rank == 0) {
		for (int i = 1; i < ntasks; ++i) {
			showMessage("Sending kill signal to slave " + str(i) + "\n",4,fixedParams);
			//mpiChannel.setMessageId(dietag);
			//mpiChannel.setMessageRank(i);
			//mpiChannel << 0;
		}
	}
	if (localErrorValue != NULL) delete localErrorValue;	                                                                        
}



void MPIErrorValueCalculator::calculateErrorValue(ModelTuningParameters & params) {

	vector< ModelTuningParameters > paramList(1);
	paramList[0] = params;
	
	calculateParallelErrorValue(paramList);

	/// This is necessary because otherwise the error value is not transferred
	/// since no reference is passed to calculateParallelErrorValue
	params.setErrorValue(paramList[0].getErrorValue());
}


void MPIErrorValueCalculator::calculateParallelErrorValue(vector< ModelTuningParameters > & paramList) {

	showMessage("Running a total of "+ str((int)paramList.size()) + " jobs on " + str(ntasks-1) + " parallel processors\n",3,fixedParams);

	int nSubmitted = 0;
	int nReceived = 0;
	int taskRank = 1;

	//Run models on all available slaves
	while (nSubmitted < (int)ntasks-1 && nSubmitted < (int)paramList.size()) {
		runErrorValueOnSlave(taskRank++, nSubmitted, paramList[nSubmitted]);
		nSubmitted++;
	}

	//There are more jobs than slaves
	while (nSubmitted < (int)paramList.size()) {
		receiveErrorValueFromSlave(taskRank, paramList);
		nReceived++;
		runErrorValueOnSlave(taskRank, nSubmitted ,paramList[nSubmitted]);
		nSubmitted++;
	}

	//Receive the remainder of the results
	while (nReceived < nSubmitted) {
		receiveErrorValueFromSlave(taskRank, paramList);
		nReceived++;
	}

	numberOfGenerations++;

}


void MPIErrorValueCalculator::runErrorValueOnSlave(int slaveNumber, int resultNumber, const ModelTuningParameters params) 
{
    
    showMessage("Sending parameters to slave: " + str(slaveNumber) + "... ",4,fixedParams);

    
    mpiChannel.setMessageRank(slaveNumber);
    mpiChannel << resultNumber;
    params.printOn(mpiChannel);
    showMessage(" Parameters sent \n",4,fixedParams);

}

void MPIErrorValueCalculator::receiveErrorValueFromSlave(int & taskRank, vector< ModelTuningParameters > & paramList) {

    int resultNumber;

    showMessage("Waiting for error value from slave ... \n",4,fixedParams);
    mpiChannel.setMessageRank(MPI_ANY_SOURCE);
    mpiChannel >> resultNumber;
    taskRank = mpiChannel.getMessageRank();

    showMessage("Receiving result " + str(resultNumber) + " from slave " + str(taskRank) + "... ",4,fixedParams);

    paramList[resultNumber].readFrom(mpiChannel);

	numberOfEvaluations++;

	if (exportFileStream.is_open()) {
			exportFileStream << numberOfGenerations << " " << numberOfEvaluations << " " << paramList[resultNumber].getErrorValue() << " ";
			for (int j = 0; j < paramList[resultNumber].getLength(); j++) {
				exportFileStream << (paramList[resultNumber][j]) << " ";
			}
			exportFileStream << endl;
   	}
                                                                        
    showMessage(" Error value received\n",4,fixedParams);

}









class runModelMethod : public xmlrpc_c::method {
public:
    runModelMethod() {
        // signature and help strings are documentation -- the client
        // can query this information with a system.methodSignature and
        // system.methodHelp RPC.
        this->_signature = "A:A"; 
        this->_help = "Runs the model given the parameters";
    }

	void execute(xmlrpc_c::paramList const& paramList, xmlrpc_c::value *   const  retvalP) 
	{
		ModelTuningParameters parameters;
		readModelParametersFromArray(paramList, parameters);

		cerr << "Slave " <<  (rank)  <<  " running model with parameters: "  <<  parameters.toString() <<  endl;
	    
		localErrorValue->calculateErrorValue(parameters);

		cerr << "Slave " << (rank) << " sending error value back to master ... " << endl;

		writeModelParametersIntoArray(parameters, resultsArray);
			
		mpiChannel.setMessageRank(0);
		mpiChannel << resultNumber;
		parameters.printOn(mpiChannel);
		
		showMessage("Slave " + str(rank) + " has sent error value back\n",4,fixedParams);            
		*
		
        int const addend(paramList.getInt(0));
        int const adder(paramList.getInt(1));
        
        paramList.verifyEnd(2);
        
        *retvalP = xmlrpc_c::value_int(addend + adder);
    }
};


static void
setupSignalHandlers(void) {

/* In UNIX, when you try to write to a socket that has been closed
 *  from the other end, your write fails, but you also get a SIGPIPE
 *  signal.  That signal will kill you before you even have a chance
 *  to see the write fail unless you catch, block, or ignore it.
 *  If a client should connect to us and then disconnect before we've
 *  sent our response, we see this socket-closed behavior.  We
 *  obviously don't want to die just because a client didn't complete
 *  an RPC, so we ignore SIGPIPE.
 **/
    struct sigaction mysigaction;
    sigemptyset(&mysigaction.sa_mask);
    mysigaction.sa_flags = 0;
    mysigaction.sa_handler = SIG_IGN;
    sigaction(SIGPIPE, &mysigaction, NULL);
}



void MPIErrorValueCalculator::startSlave() 
{
    xmlrpc_c::registry myRegistry;
    xmlrpc_c::methodPtr const runModelMethodP(new runModelMethod);
    myRegistry.addMethod("neurofitter.runModel", runModelMethod);

    xmlrpc_c::serverAbyss neurofitterSlave(
        myRegistry,
        1980,              // TCP port on which to listen
        "/tmp/xmlrpc_log"  // Log file
        );

	cerr << "Starting neurofitter slave " << rank <<  endl;

    while (true) {

	cout << "Slave "  << (this->rank)  << " waiting for next rpc call " << endl;

        neurofitterSlave.runOnce();
	/* This waits for the next connection, accepts it, reads the
	 * HTTP POST request, executes the indicated RPC, and closes
	 * the connection.
	 */

    }
	showMessage("Slave " + str(this->rank) + " exiting ended\n",4,fixedParams);
    return 0;

	/* 
        ModelTuningParameters parameters;
        parameters.readFrom(mpiChannel);

        showMessage("Slave " + str(rank) + " running model with parameters: " + parameters.toString() + "\n",3,fixedParams);
    
        localErrorValue->calculateErrorValue(parameters);

        showMessage("Slave " + str(rank) + " sending error value back to master ... ",4,fixedParams);
                
        mpiChannel.setMessageRank(0);
        mpiChannel << resultNumber;
        parameters.printOn(mpiChannel);
        
        showMessage("Slave " + str(rank) + " has sent error value back\n",4,fixedParams);            
	*/
    
}




vector< pair< ModelTuningParameters, double > > MPIErrorValueCalculator::getErrorValueHistory() {
	return errorValueHistory;
}


void MPIErrorValueCalculator::enableFileExport(const string fileName) {
	exportFileStream.open(fileName.c_str(), ios::out);
	
	showMessage("MPIErrorValueCalculator: Enabled export to file: " + fileName + "\n",3,fixedParams);            
}

void MPIErrorValueCalculator::disableFileExport() {
	exportFileStream.close();
}
   
string MPIErrorValueCalculator::getExportFile() const {
	return exportFile;
}
