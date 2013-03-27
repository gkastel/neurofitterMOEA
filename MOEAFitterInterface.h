/*
 *
Revision of last commit: $Rev$
Author of last commit: $Author$
Date of last commit: $Date$
*/

#ifndef NEUROFITTER_MOEAFITTERINTERFACE_H
#define NEUROFITTER_MOEAFITTERINTERFACE_H

#include "FitterResults.h"
#include "FitterInterface.h"
#include "ErrorValueCalculator.h"
#include "FixedParamObject.h"
#include "MOEAErrorValueCalculator.h"
#include "MersenneTwister.h"

using namespace std;

/* 
 * Multiobjective optimization evolutionary algorithm
 * This optimization combines a number of error measures, calculated with the assorted 
 * MOEAErrorValueCalculator  and tries to approach  the pareto-optimal solutions of the
 * problem space. There is no single objective function here. The results of running the
 * algorithm will be a set of points describing the pareto-optimal sets of parameters. 
 * The pareto-optimal set that is closest to (0,0) point is given back as an optimal 
 * solution. This solution is the best compromise between objectives. However, all the 
 * other solutions in the final generation of the evolutionary algorithm are also optimal 
 * for some objectives, but not for all.  The solutions generated in the final generation 
 * of the algorithm are written in a file. You can examine the solutions by plotting the 
 * error values in pairs.
 *
 * The method of Pareto-optimization is Nondominated-Sorting II described by Deb et.al., 1995. 
 * The evolution algorithm is eoSGA, the simple genetic algorithm
 * Settings (neurofitter.xml): 
 
  <FitterType>MOEA</FitterType>
  <FitterParameters>

        <PopulationSize>100</PopulationSize>              <!-- Size of points in the pareto front to optimize -->

        <NumberOfOffspring>90</NumberOfOffspring> 	<!-- Offspring generated in each generation --->

        <MutationStep>0.01</MutationStep>		<!-- Multiplication factor for mutating parameters in each generation -->

        <MutationProbability>70</MutationProbability>   <!-- Probability of mutation for parameters  Warning: value is in percentage (0-100.0%) (i.e. a model parameter will be mutated with a probability=MutationProbability% and new_value=old_value*MutationStep*random(0,1.0) -->

        <MaxGenerations>2000</MaxGenerations>		<!-- How many generations to run -->

        <OutputFile>3</OutputFile>		<!-- Write the results of optimization  in an output file  -->

	<!--  The format of the output file is: 
	Param-1  Param-2 Param-3 .... ErrorValue1 ErrorValue2 ErrorValue3 ...
	Param-1  Param-2 Param-3 .... ErrorValue1 ErrorValue2 ErrorValue3 ...
	with one line for each individual in the final population (each individual represents a pareto-optimal 
	solution 
	-->
  </FitterParameters>



 *
 * MOEA fitter *must* be used only with MOEA error value calculator.
 */
class MOEAFitterInterface : public FitterInterface, FixedParamObject
{

	public:
	MOEAFitterInterface(ErrorValueCalculator * fit, FixedParameters params):FitterInterface(fit), FixedParamObject(params) {};
	//Inherited from FitterInterface
	virtual FitterResults runFitter(ModelTuningParameters * resolution);

	private:
	/* Settings for this fitter */
	int numberOfOffspring, populationSize, maxGenerations;
	double mutationStep, mutationProbability;


};

#endif



