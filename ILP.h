//
// Created by Igor on 11/11/2026.
//

#ifndef PLI_ILP_H
#define PLI_ILP_H

#include "TAD_AGraph.h"
#include "TAD_UELGraph.h"

/* Resolve a formulação PLI da imagem usando CPLEX */
void solve_ilp_cplex(const AGraph *ag, const UELGraph *g1, const UELGraph *g2);

#endif //PLI_ILP_H
