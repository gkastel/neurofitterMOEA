/*
 * Multiobjective evolutionary optimization implementation
 * This is based on t-eoPareto.cpp from the tests directory of EO 1.0
 */
 

#include <eo>

//#include <utils/eoMOFitnessStat.h>
#include <eoNDSorting.h>
#include <eoParetoFitness.h>
#include "../MOEAFitterInterface.h"


using namespace std;

// Look: overloading the maximization without overhead (thing can be inlined)
class MinimizingFitnessTraits : public eoParetoFitnessTraits
{
  public :
  static bool maximizing(int) { return false; }
};

typedef eoParetoFitness<MinimizingFitnessTraits> fitness_type;


struct individual : public EO<fitness_type>
{
  vector<double> paramValues;
};


/* Mutate an individual */	
struct Mutator : public eoMonOp<individual>
{
	
	ModelTuningParameters* modelParameters;
	double mutationFactor;
	double mutationProbability;
	Mutator(ModelTuningParameters* parameters, double mf, double mp)
	: modelParameters(parameters), mutationFactor(mf), mutationProbability(mp)
	{ }

	
	bool operator() (individual& p)
	{
	    double minval, maxval;

	    for (unsigned i = 0; i < p.paramValues.size(); ++i)
	    {
		minval = modelParameters->getLowerBound(i);
		maxval = modelParameters->getUpperBound(i);

		if (rng.flip(mutationProbability))
		{
			p.paramValues[i] += rng.normal() * mutationFactor * p.paramValues[i];
		}

		if (p.paramValues[i] < minval)
		{
		        p.paramValues[i] = minval;
			cerr << "*Warning: Parameter " << i <<" exceeded lower bound in mutation, consider a smaller mutation step" << endl;
		}
		else if (p.paramValues[i] > maxval)
		{
		        p.paramValues[i] = maxval;
			cerr << "*Warning: Parameter " << i <<" exceeded upper bound in mutation, consider a smaller mutation step" << endl;
		}
	    }
	    return true;
	}
};


struct BatchEvaluator : public eoPopEvalFunc<individual>
{
	ErrorValueCalculator* calculator;

	BatchEvaluator(ErrorValueCalculator* calc) : calculator(calc) 
	{
	}

	void operator() (eoPop<individual>& parents, eoPop<individual>& pop)
	{
		assert(calculator != NULL);

		vector< vector<double> > errorValues;

		vector< ModelTuningParameters > runQueue;
		vector<int> runQueueIds;
		for (unsigned j =0; j < pop.size(); j++)
		{
			if (pop.at(j).invalid())
			{
				ModelTuningParameters mtp;
				mtp.setTuningParameters(pop.at(j).paramValues);
				runQueue.push_back(mtp);
				runQueueIds.push_back(j);
			}
		}


		cerr <<"Evaluating " << runQueue.size() << " individuals in parallel" << endl;
		calculator->calculateParallelErrorValue(runQueue);
		
		for (unsigned j =0; j < runQueue.size(); j++)
		{
			int index = runQueueIds.at(j);
			fitness_type errors;
			vector<double>& objectiveValues = runQueue.at(j).getErrorValuesVector();

			if (!objectiveValues.size())
			{
				throw std::runtime_error("Error values vector is empty! MOEAFitterInterface Only works with MOEAErrorValueCalculator (or any multi-objective calculator, for that matter)!");
			}

			errors.resize(objectiveValues.size());
			for (unsigned k =0; k  < objectiveValues.size(); k++)
				errors.at(k) = objectiveValues.at(k); 

			/*
			cerr <<index <<" objectives size= " << objectiveValues.size() << " errors size=" << errors.size();
			for (unsigned k =0; k  < errors.size(); k++)
				cerr << " " << errors.at(k);
			cerr  << endl;
			*/
			pop.at(index).fitness(errors);
		}
	}
};


struct Initializer : public eoInit<individual>
{
	ModelTuningParameters * modelParameters;	

	Initializer(ModelTuningParameters* params) : modelParameters(params)
	{}


	void operator() (individual& p)
	{
		p.paramValues.resize(modelParameters->getLength());
		for (int k =0; k < modelParameters->getLength(); k++)
		{
			float v = modelParameters->getLowerBound(k) + (rng.uniform(1.0))*(modelParameters->getUpperBound(k) - modelParameters->getLowerBound(k));
			p.paramValues.at(k) = v;
		}
		p.invalidate();
	}
};


eoPerf2Worth<individual, double>& make_perf2worth(eoState& state)
{
	return state.storeFunctor(new eoNDSorting_II<individual>());
}


/**
 * Trying out an elitist non-dominated sorted replacement scheme.
 */
class eoNDPlusReplacement : public eoReplacement<individual>
{
public:

    // using eoNDPlusReplacement< EOT, WorthT >::first;

    eoNDPlusReplacement(eoPerf2Worth<individual, double>& _perf2worth)
        : perf2worth(_perf2worth)
        {}

    struct WorthPair : public pair<double, const individual*>
    {
        bool operator<(const WorthPair& other) const
            { return other.first < this->first; }
    };


  void operator()(eoPop<individual>& _parents, eoPop<individual>& _offspring)
  {
  	cerr << " Replacing " << _parents.size() <<" with " << _offspring.size();

    unsigned sz = _parents.size();
    _parents.reserve(_parents.size() + _offspring.size());
    std::copy(_offspring.begin(), _offspring.end(), back_inserter(_parents));

    // calculate worths
    perf2worth(_parents);
    perf2worth.sort_pop(_parents);
    cerr <<" Sorted population " << endl;
    for (unsigned k=0; k < _parents.size(); k++)
    {
    	fitness_type ft = _parents.at(k).fitness();
	cerr << k <<" { ";
	for (unsigned l =0; l < ft.size(); l++)
	{
		cerr << ft.at(l) << " " ;
	}
	cerr << "} " << endl;
    }

	cerr << "Cutoff at " << sz <<endl;
    perf2worth.resize(_parents, sz);

    _offspring.clear();
  }

private :
    eoPerf2Worth<individual, double>& perf2worth;
};



FitterResults MOEAFitterInterface::runFitter(ModelTuningParameters * startingPoint) 
{
	this->numberOfOffspring =  toInt(fixedParams["NumberOfOffspring"]);
	this->populationSize =  toInt(fixedParams["PopulationSize"]);
	this->mutationStep =  toDouble(fixedParams["MutationStep"]);
	this->mutationProbability =  toDouble(fixedParams["MutationProbability"]);
	this->maxGenerations =  toInt(fixedParams["MaxGenerations"]);

	try 
	{
		//MOEAErrorValueCalculator* calculator = dynamic_cast<MOEAErrorValueCalculator*>(this->errorValue);
		//if (!calculator)
		//	throw std::runtime_error(__FILE__": MOEAFitterInterface can only be used with the MOEAErrorValueCalculator, sorry");

		eoState state;
		Initializer init(startingPoint);
		BatchEvaluator batchEval(this->errorValue);
		Mutator mutate(startingPoint, this->mutationStep, this->mutationProbability);

		eoPop<individual> pop(this->populationSize, init);

		eoPerf2Worth<individual,double>& perf2worth = make_perf2worth(state);
		
		eoStochTournamentWorthSelect<individual,double> select(perf2worth, (double)0.95);

		// One general operator
		eoProportionalOp<individual> opsel;
		opsel.add(mutate, 1.0);

		// the breeder
		eoGeneralBreeder<individual> breeder(select, opsel);

		// replacement
		eoNDPlusReplacement replace(perf2worth);

		unsigned long generation = 0;
		eoGenContinue<individual> gen(this->maxGenerations, generation);
		eoCheckPoint<individual> cp(gen);

		//eoMOFitnessStat<individual> fitness0(0, "FirstObjective");
		//eoMOFitnessStat<individual> fitness1(1, "SecondObjective");

		//cp.add(fitness0);
		//cp.add(fitness1);

		//eoGnuplot1DSnapshot snapshot("pareto");
		//snapshot.with(eoGnuplot::Points(3));
		//cp.add(snapshot);
		//snapshot.add(fitness0);
		//snapshot.add(fitness1);

		eoEasyEA<individual> ea(cp, batchEval, breeder, replace);
		
		// turns out, i like short names
		eoPop<individual> dummies;
		batchEval(dummies, pop);
		ea(pop);


		// Final population reached 
		string& s = fixedParams["OutputFile"];
		if (s.length())
		{
			ofstream out(s.c_str(), ios::out);
			for (unsigned i=0; i < pop.size(); i++)
			{
				for (unsigned k=0;k < pop.at(i).paramValues.size(); k++)
					out << pop.at(i).paramValues.at(k) << " " ;

				fitness_type f = pop.at(i).fitness();
				for (unsigned k=0;k < f.size(); k++)
					out << f.at(k) << " " ;
				out << endl;

				for (unsigned k=0;k < f.size(); k++)
					out << f.at(k) << " " ;


			}
			out.close();
		}

		//We will  be returning the individual with all error values closest to zero
		individual bestFit;
		double minDistance = 999999999.0;
		int minIndex =-1;

		for (unsigned i=0; i < pop.size(); i++)
		{
			individual& p = pop.at(i);
			double d =0;
			for (unsigned k=0; k<p.fitness().size(); k++)
				d += p.fitness().at(k)*p.fitness().at(k);
			if (d < minDistance)
			{
				minDistance = d;
				minIndex = i;
			}
		}

		cerr<< "Warning: the best fit returned from MOEAFitter will not be the smallest weighted sum, but the smallest distance from axes' Zero"<<endl;
		if (minIndex != -1)
		{
			ModelTuningParameters tmp;
			tmp.setTuningParameters(pop.at(minIndex).paramValues);
			FitterResults fr;
			fr.setBestFit(tmp, minDistance);
			return fr;
		}
		else
			return FitterResults();
	}
	catch (std::exception& the_ball)
	{
	    std::cerr << "Exceptions happen. " << the_ball.what() << std::endl;
	    throw the_ball;
	}
}



