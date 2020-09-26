#include "init.h"

template <typename T>
void SearchInitializer::bayesian_update(T& c)
{
    double coeff_product = 1.0;
    for( int j=0; j<c.size(); j++ )
    {
        Lit l = c[j];
        if ( solver->value(l) == l_True ) return; // this clause is already satisfied
        if ( solver->value(l) == l_False ) continue; // this literal is already falsified
        int v = var(l);
        int sgn = !sign(l);

        double coeff = (sgn ? parameters[v].b : parameters[v].a) / (parameters[v].a + parameters[v].b);
        coeff_product *= coeff;
    }
    double normalization_constant = 1 - coeff_product;

    for( int j=0; j<c.size(); j++ )
    {
        Lit l = c[j];
        if ( solver->value(l) == l_False ) continue; // this literal is already falsified
        int v = var(l);
        int sgn = !sign(l);

        double *p[2], up[2];
        p[0] = &parameters[v].a;
        p[1] = &parameters[v].b;

        // updated parameters
        up[0] = parameters[v].a + (!sgn);
        up[1] = parameters[v].b + (sgn);

        double sumP = *p[0] + *p[1];
        double sumUP = up[0] + up[1];

        for( int k=0; k<=1; k++ )
        {
            double moment1 = *p[k] / sumP - coeff_product * up[k] / sumUP;
            double moment2 = *p[k] * (*p[k] + 1) / ((sumP) * (sumP+1)) - coeff_product * up[k] * (up[k] + 1) / ((sumUP) * (sumUP+1));

            moment1 /= normalization_constant;
            moment2 /= normalization_constant;

            *p[k] = moment1 * (moment1 - moment2) / (moment2 - moment1*moment1);
        }
    }
}

void SearchInitializer::bayesian()
{
    int n = solver->nVars();

    for( int k=0; k<init_epochs; k++ )
    {
        for( int i=0; i<solver->nClauses(); i++ )
        {
            Clause& c = solver->getClause(i);
            bayesian_update(c);
        }

/*        for( int i=0; i<learnts_core.size(); i++ )
        {
            Clause& c = ca[learnts_core[i]];
            bayesian_update(c);
        }

        for( int i=0; i<learnts_tier2.size(); i++ )
        {
            Clause& c = ca[learnts_tier2[i]];
            bayesian_update(c);
        }

        for( int i=0; i<learnts_local.size(); i++ )
        {
            Clause& c = ca[learnts_local[i]];
            bayesian_update(c);
        }*/
    }

    if ( polarity_init_method == BMM )
        for( int v=0; v<n; v++ )
            solver->setPolarity(v, (parameters[v].a > parameters[v].b) ? false : true);

    if ( activity_init_method == BMM )
        for( int v=0; v<n; v++ )
            solver->setActivity(v, BayesianWeight(v));
}

void SearchInitializer::init_bayesian()
{
    int n = solver->nVars();
    parameters.resize(n);
    for( int i=0; i<n; i++ )
    {
        parameters[i].a = (double)rand() / RAND_MAX + 10;
        parameters[i].b = (double)rand() / RAND_MAX + 10;
    }
}

void SearchInitializer::jeroslow_wang()
{
    int n = solver->nVars();
    vector<double> cnt(2 * n, 0.0);

    for( int i=0; i<solver->nClauses(); i++ )
    {
        Clause& c = solver->getClause(i);
        if ( c.size() >= 50 ) continue;
        double sc = pow(2, -c.size());
        for( int j=0; j<c.size(); j++ )
        {
            Lit q = c[j];
            if ( sign(q) )
                cnt[var(q) + n] += sc;
            else
                cnt[var(q)] += sc;
        }
    }

    if ( polarity_init_method == JW )
    {
        for( int v=0; v<n; v++ )
            solver->setPolarity(v, (cnt[v] > cnt[v+n]) ? false : true);
    }

    if ( activity_init_method == JW )
    {
        int m = solver->nClauses();
        if ( m == 0 ) m = 1;
        for( int v=0; v<n; v++ )
            solver->setActivity(v, (cnt[v] + cnt[v+n]) / m);
    }
}

