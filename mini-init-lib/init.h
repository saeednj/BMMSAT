#ifndef _INITIALIZER_H_
#define _INITIALIZER_H_

#include <vector>
#include <map>
#include "core/SolverTypes.h"
#include "core/Solver.h"

using namespace std;
using namespace Minisat;

struct BetaDist {
    double a, b;
};

#define BayesianWeight(v) (max(parameters[v].a, parameters[v].b) / (parameters[v].a + parameters[v].b))


/** Polarity and Activity Initializer for MiniSat based SAT solvers
 * Current methods:
 * 1. Default (All-zero / All-False)
 * 2. Bayesian Moment Matching
 * 3. Jeroslow-Wang
 * 4. Random
 */

class SearchInitializer {
    public:

        enum InitMethod {
            DEFAULT,
            BMM,
            JW,
            RANDOM,
        };

        SearchInitializer(Solver* s, int polarity, int activity, int initEpoch, int updateEpoch)
        {
            srand(time(NULL));

            solver = s;

            polarity_init_method = (InitMethod)polarity;
            activity_init_method = (InitMethod)activity;

            init_epochs = initEpoch;
            update_epochs = updateEpoch;

            if ( polarity_init_method == BMM || activity_init_method == BMM ) {
                init_bayesian();
                bayesian();
            }

            if ( polarity_init_method == JW || activity_init_method == JW )
                jeroslow_wang();

            if ( polarity_init_method == RANDOM )
                for(int v=0; v<solver->nVars(); v++)
                    solver->setPolarity(v, rand() % 2 ? false : true);

            if ( activity_init_method == RANDOM )
                for(int v=0; v<solver->nVars(); v++)
                    solver->setActivity(v, rand() / RAND_MAX * 0.00001);
        }

        virtual ~SearchInitializer() {}

        // TODO update interface

    private:
        void bayesian();
        template<typename T>
        void bayesian_update(T& c);
        void init_bayesian();

        void jeroslow_wang();
        
        Solver* solver;

        InitMethod polarity_init_method, activity_init_method;

        // BMM
        vector<BetaDist> parameters;
        int init_epochs;
        int update_epochs;
};
#endif
