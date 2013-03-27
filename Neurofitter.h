/*
Revision of last commit: $Rev: 241 $     
Author of last commit: $Author: werner $  
Date of last commit: $Date: 2007-10-22 16:17:53 +0900 (Mon, 22 Oct 2007) $    
*/

#ifndef NEUROFITTER_H
#define NEUROFITTER_H

#define WITH_MOEA 1

#ifdef WITH_MPI
	#include <mpi.h>
	#include "MPIModelInterface.h"
	#include "MPIErrorValueCalculator.h"
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

#include "DataTrace.h"
#include "XMLString.h"
#include "FixedParameters.h"

#include "TracesReader.h"
#include "NormalTracesReader.h"

#include "ModelTuningParameters.h"
#include "ModelResults.h"
#include "ExperimentInterface.h"
#include "ModelInterface.h"
#include "ErrorValueCalculator.h"
#include "FitterInterface.h"

#include "GenesisModelInterface.h"
#include "NeuronModelInterface.h"
#include "ExecutableModelInterface.h"
#include "FakeExperimentInterface.h"
#include "FileExperimentInterface.h"
#include "FileListExperimentInterface.h"
#include "VdVdtMatrix.h"
#include "MapVdVdtMatrix.h"
#include "DirectVdVdtMatrix.h"
#include "DistanceVdVdtMatrix.h"
#include "MatrixErrorValueCalculator.h"
#include "EasyFitterInterface.h"
#include "MeshFitterInterface.h"
#include "RandomFitterInterface.h"
#include "FileFitterInterface.h"
#include "SwarmFitterInterface.h"

#ifdef WITH_EO
	#include "EOErrorValueCalculator.h"
	#include "EOFitterInterface.h"
#endif

#ifdef WITH_NOMAD
	#include "NOMADFitterInterface.h"
#endif

#ifdef WITH_EO
	#ifdef WITH_NOMAD
		#include "EONOMADFitterInterface.h"
	#endif
#endif

#ifdef WITH_MOEA
#include "MOEAErrorValueCalculator.h"
#include "MOEAFitterInterface.h"
#endif





#endif
