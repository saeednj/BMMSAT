#ifndef Minisat_WDimacs_h
#define Minisat_WDimacs_h

#include <stdio.h>

#include "utils/ParseUtils.h"
#include "core/SolverTypes.h"

namespace Minisat {

//=================================================================================================
// WDIMACS Parser:

template<class B, class Solver>
static void readWClause(B& in, Solver& S, vec<Lit>& lits, uint64_t& weight) {
    int     parsed_lit, var;
    weight = parseInt64(in);
    assert(weight >= 1);
    lits.clear();
    for (;;){
        parsed_lit = parseInt(in);
        if (parsed_lit == 0) break;
        var = abs(parsed_lit)-1;
        while (var >= S.nVars()) S.newVar();
        lits.push( (parsed_lit > 0) ? mkLit(var) : ~mkLit(var) );
    }
}

template<class B, class Solver>
static void parse_WDIMACS_main(B& in, Solver& S) {
    WeightedClause wc;
    //vec<Lit> lits;
    int vars    = 0;
    int clauses = 0;
    int cnt     = 0;
    uint64_t top     = 0;
    //uint64_t weight;
    for (;;){
        skipWhitespace(in);
        if (*in == EOF) break;
        else if (*in == 'p'){
            if (eagerMatch(in, "p wcnf")){
                vars    = parseInt(in);
                clauses = parseInt(in);
                top     = parseInt64(in);
                // SATRACE'06 hack
                // if (clauses > 4000000)
                //     S.eliminate(true);
            }else{
                printf("PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
            }
        } else if (*in == 'c' || *in == 'p')
            skipLine(in);
        else{
            cnt++;
            readWClause(in, S, wc.lits, wc.weight);
            if ( wc.weight == top ) // Hard clause?
                S.addClause_(wc.lits);
            else
                S.addSoftClause(wc);
        }
    }
    if (vars != S.nVars())
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of variables.\n");
    if (cnt  != clauses)
        fprintf(stderr, "WARNING! DIMACS header mismatch: wrong number of clauses.\n");
}

// Inserts problem into solver.
//
template<class Solver>
static void parse_WDIMACS(gzFile input_stream, Solver& S) {
    StreamBuffer in(input_stream);
    parse_WDIMACS_main(in, S); }

//=================================================================================================
}

#endif
