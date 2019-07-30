#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <sstream>
using namespace std;

typedef vector<int> Clause;

void parse(vector<Clause>& cs, int &nVars)
{
    int nc;
    char line[2005];
    while( fgets(line, 2000, stdin) != NULL )
    {
        if ( line[0] == 'c' ) continue;

        if ( line[0] == 'p' )
            sscanf(line, "p cnf %d %d", &nVars, &nc);
        else
        {
            istringstream iss(line);
            int x;
            Clause c;
            while( iss >> x )
            {
                if ( x == 0 ) break;
                c.push_back(x);
            }
            cs.push_back(c);
        }
    }
    assert(nc == cs.size());
}

int main()
{
    int nVars;

    int epochs = 100;

    srand(time(NULL));

    vector<Clause> cs;
    parse(cs, nVars);

    vector<double> dummy(2, 0.0);
    vector<vector<double>> parameters(nVars, dummy);
    vector<vector<double>> updatedParams(nVars, dummy);

    for( auto &x : parameters )
    {
        x[0] = (double)rand() / RAND_MAX + 1;
        x[1] = (double)rand() / RAND_MAX + 1;
    }
//    parameters[0][0] = 1.05;
//    parameters[1][0] = 10;

//    printf("%lf\n", lgamma(1000));

    vector<int> assign(nVars, 0);
    int best = cs.size() + 1;
    vector<vector<double>> bestSolution;
    int bestK;

//    for( double x=1; x<101; x*=10)
//        printf("gamma %.2lf: %.6lf\n", x, tgamma(x));
    for( int k=0; k<epochs; k++ )
    {
        for( Clause &c : cs )
        {
            double coeff_product = 1.0;
            for( int lit : c )
            {
                int v = abs(lit) - 1;
                int sign = (lit > 0);

                updatedParams[v][0] = parameters[v][0] + (!sign);
                updatedParams[v][1] = parameters[v][1] + (sign);

                double coeff = (sign ? parameters[v][1] : parameters[v][0]) / (parameters[v][0] + parameters[v][1]);
                coeff_product *= coeff;
            }
            double normalization_constant = 1 - coeff_product;

            for( int lit : c )
            {
                int v = abs(lit) - 1;
                int sign = (lit > 0);

                double sumP = parameters[v][0] + parameters[v][1];
                double sumUP = updatedParams[v][0] + updatedParams[v][1];

                for( int i=0; i<=1; i++ )
                {
                    double moment1 = parameters[v][i] / sumP - coeff_product * updatedParams[v][i] / sumUP;
                    double moment2 = parameters[v][i] * (parameters[v][i] + 1) / ((sumP) * (sumP+1)) - coeff_product * updatedParams[v][i] * (updatedParams[v][i] + 1) / ((sumUP) * (sumUP+1));

                    moment1 /= normalization_constant;
                    moment2 /= normalization_constant;

                    parameters[v][i] = moment1 * (moment1 - moment2) / (moment2 - moment1*moment1);
                }
            }
        }

        for( int v=0; v<nVars; v++ )
        {
            assign[v] = (parameters[v][0] > parameters[v][1]) ? 1 : -1;
        }
        
        int unsat = 0;
        for( Clause& c : cs )
        {
            bool sat = false;
            for( int lit : c )
                if ( assign[abs(lit)-1] * (lit) > 0 ) { sat = true; break; }
            if ( !sat ) unsat++;
        }

        printf("epoch %d : unsat clauses = %d\n", k, unsat);

        if ( unsat < best )
        {
            best = unsat;
            bestSolution = parameters;
            bestK = k + 1;
        }
        if ( unsat == 0 ) break;
        
    }

/*    for( int t=0; t<(1<<nVars); t++ )
    {
        for( double x=0.05; x<0.99; x+=0.05 )
        { 
            double total_log_prob = 0.0;
            double total_prob = 1.0;
            for( int v=0; v<nVars; v++ )
            {
                double val = ( (t & (1 << v)) > 0 ) ? x : 1-x;
                double a = parameters[v][0];
                double b = parameters[v][1];
                double logp = lgamma(a + b) - lgamma(a) - lgamma(b) + (a-1)*log(val) + (b-1)*log(1-val);
                double p = (tgamma(a + b) / (tgamma(a) * tgamma(b))) * pow(val, a-1) * pow(1-val, b-1);
                total_log_prob += logp;
                total_prob *= p;
                //            printf("(DBG: a=%.6lf, b=%.6lf, p=%.6lf) ", a, b, p);
                printf("%3d ", ( (t & (1 << v)) > 0 ) ? v+1 : -(v+1));
            }
            //printf("%.6lf\n", exp(total_log_prob));
            printf("%.2lf %.6lf\n", x, total_log_prob);
        }
    }*/

    if ( best == 0 )
    {
        printf("SAT! converged after %d steps:\n", bestK);
        for( int v=0; v<nVars; v++ )
            printf("%d ", (v+1) * assign[v]);
        printf("\n");
    }
    else
    {
        printf("UNDET! converged after %d steps:\n", bestK);
        printf("at least %d clauses remained unsat (out of %d)\n", best, cs.size());
    }

//    for( int v=0; v<nVars; v++ )
//        printf("%d : %.6lf\n", v+1, bestSolution[v][0] / (bestSolution[v][0] + bestSolution[v][1]));
//    printf("Parameters:\n");
//    for( int v=0; v<nVars; v++ )
//        printf("%d : %.6lf\n", v+1, parameters[v][0] / (parameters[v][0] + parameters[v][1]));

    return 0;
}

