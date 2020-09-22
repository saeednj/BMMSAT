#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

#define COUNT_UNSAT true
//#define MULTI_SOLUTION true
#define MAIN_BAYESIAN true

typedef vector<int> Clause;
struct BetaDist { double a, b; };

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

class Bayesian {
    public:
        Bayesian(int n, vector<Clause>& cs)
            :nvars(n), clauses(cs)
        {
            parameters.resize(nvars);
            updatedParams.resize(nvars);
        }
        ~Bayesian(){}

        void init(int type = 0);
        void update(Clause& c);
        int main(int epochs);
        void getAssignment(vector<int>& assign);

        vector<Clause>& clauses;
        int nvars;
        vector<BetaDist> parameters;
        vector<BetaDist> updatedParams;
};

void Bayesian::init(int type)
{
    if ( type == 0 )
    {
        for( int i=0; i<nvars; i++ )
        {
            parameters[i].a = (double)rand() / RAND_MAX + 1;
            parameters[i].b = (double)rand() / RAND_MAX + 1;
        }
    }
    else if ( type == 1 )
    {
        vector<int> cnt(nvars, 0);

        for( Clause& c : clauses )
        {
            for( int j=0; j<c.size(); j++ )
            {
                int v = abs(c[j]) - 1;
                if ( c[j] > 0 )
                    cnt[v]++;
                else
                    cnt[v]--;
            }
        }

        for( int i=0; i<nvars; i++ )
        {
            if ( cnt[i] > 0 )
            {
                parameters[i].a = 10 - (double)rand() / RAND_MAX;
                parameters[i].b =  0 + (double)rand() / RAND_MAX;
            }
            else
            {
                parameters[i].a =  0 + (double)rand() / RAND_MAX;
                parameters[i].b = 10 - (double)rand() / RAND_MAX;
            }
        }
    }
    else
    {
        fprintf(stderr, "Wrong init mode\n");
        abort();
    }
}

void Bayesian::update(Clause &c)
{
    double coeff_product = 1.0;
    for( int lit : c )
    {
        int v = abs(lit) - 1;
        int sgn = (lit > 0);

        updatedParams[v].a = parameters[v].a + (!sgn);
        updatedParams[v].b = parameters[v].b + (sgn);

        double coeff = (sgn ? parameters[v].b : parameters[v].a) / (parameters[v].a + parameters[v].b);
        coeff_product *= coeff;
    }
    double normalization_constant = 1 - coeff_product;

    for( int lit : c )
    {
        int v = abs(lit) - 1;

        double sumP = parameters[v].a + parameters[v].b;
        double sumUP = updatedParams[v].a + updatedParams[v].b;

        double *p[2], *up[2];
        p[0] = &parameters[v].a;
        p[1] = &parameters[v].b;
        up[0] = &updatedParams[v].a;
        up[1] = &updatedParams[v].b;
        for( int k=0; k<=1; k++ )
        {
            double moment1 = *p[k] / sumP - coeff_product * *up[k] / sumUP;
            double moment2 = *p[k] * (*p[k] + 1) / ((sumP) * (sumP+1)) - coeff_product * *up[k] * (*up[k] + 1) / ((sumUP) * (sumUP+1));

            moment1 /= normalization_constant;
            moment2 /= normalization_constant;

            *p[k] = moment1 * (moment1 - moment2) / (moment2 - moment1*moment1);
        }
    }
}



int Bayesian::main(int epochs)
{
#ifdef COUNT_UNSAT
    vector<int> assign(nvars, 0);
    int best = clauses.size() + 1;
    vector<BetaDist> bestSolution;
    int bestK;
#endif

    for( int k=0; k<epochs; k++ )
    {
        for( Clause& c : clauses )
            update(c);

#ifdef COUNT_UNSAT
        getAssignment(assign);
        int unsat = 0;
        for( Clause& c : clauses )
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
#endif
    }

#ifdef COUNT_UNSAT
    if ( best == 0 )
    {
        printf("SAT! converged after %d steps:\n", bestK);
        for( int v=0; v<nvars; v++ )
            printf("%d ", (v+1) * assign[v]);
        printf("\n");
        return 0;
    }
    else
    {
        printf("UNDET! converged after %d steps:\n", bestK);
        printf("at least %d clauses remained unsat (out of %lu)\n", best, clauses.size());
        return 1;
    }
#endif
}

void Bayesian::getAssignment(vector<int>& assign)
{
    for( int v=0; v<nvars; v++ )
    {
        assign[v] = (parameters[v].a > parameters[v].b) ? 1 : -1;
    }
}

int main(int argc, char **argv)
{
    int epochs;
    if ( argc > 1 ) epochs = atoi(argv[1]);
    else epochs = 100;

    srand(time(NULL));

    int nVars;
    vector<Clause> cs;
    parse(cs, nVars);

#ifdef MULTI_SOLUTION
    int solcnt = 10;
    vector<vector<int>> sol(solcnt, vector<int>(nVars, 0));
    int cnt = 0, res;
    do {
//        printf("Solution %d\n", cnt+1);
        Bayesian b(nVars, cs);
        b.init();
        res = b.main(epochs);
        b.getAssignment(sol[cnt]);
        cnt++;
        random_shuffle(cs.begin(), cs.end());
    } while( res == 1 && cnt < solcnt );

    if ( res == 0 )
    {
        return 0; // solution found
    }

    for( int j=0; j<sol[0].size(); j++ )
        for( int i=1; i<solcnt; i++ )
        {
            sol[0][j] += sol[i][j];
        }

//    printf("====\n");
    for( int i=0; i<sol[0].size(); i++ )
        printf("%d %d\n", i, sol[0][i]);


    return 1;
#endif

#ifdef MAIN_BAYESIAN
    Bayesian b(nVars, cs);
    b.init();
    int res = b.main(epochs);
    if ( res == 0 )
    {
        vector<int> s(nVars);
        b.getAssignment(s);
        for( int i=0; i<nVars; i++ )
            printf("%d ", (i+1)*s[i]);
        printf("\n");
    }
    else
    {
        printf("Not satisfied\n");
    }
    return 0;
#endif
}

