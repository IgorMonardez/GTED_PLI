//
// Created by Igor on 11/11/2026.
//

#include "ILP.h"
#include <ilcplex/cplex.h>
#include <stdio.h>
#include <stdlib.h>

void solve_ilp_cplex(const AGraph *ag, const UELGraph *g1, const UELGraph *g2) {
    if (!ag || !g1 || !g2) return;

    int status;
    CPXENVptr env = CPXopenCPLEX(&status);
    if (env == NULL) {
        fprintf(stderr, "Erro ao abrir ambiente CPLEX.\n");
        return;
    }

    CPXLPptr lp = CPXcreateprob(env, &status, "GTED_LowerBound");
    if (lp == NULL) {
        fprintf(stderr, "Erro ao criar problema CPLEX.\n");
        CPXcloseCPLEX(&env);
        return;
    }

    int ncols = ag->n_edges;
    int nrows = ag->n_vertices + g1->n_edges + g2->n_edges;

    /* --- 1. Definir Variáveis (x_e) --- */
    double *obj = (double *)malloc(ncols * sizeof(double));
    double *lb  = (double *)malloc(ncols * sizeof(double));
    double *ub  = (double *)malloc(ncols * sizeof(double));
    char   *ctype = (char *)malloc(ncols * sizeof(char));

    for (int i = 0; i < ncols; i++) {
        obj[i] = ag->edges[i].cost; // Custo delta(e)
        lb[i]  = 0.0;
        ub[i]  = 1.0;
        ctype[i] = CPX_BINARY;      // Variável binária (0 ou 1)
    }

    status = CPXnewcols(env, lp, ncols, obj, lb, ub, ctype, NULL);
    if (status) fprintf(stderr, "Erro ao adicionar colunas.\n");

    /* --- 2. Construir Matriz de Restrições --- */
    int *matbeg = (int *)malloc(nrows * sizeof(int));
    int *matind = (int *)malloc(4 * ncols * sizeof(int)); // Max ~4 non-zeros por coluna
    double *matval = (double *)malloc(4 * ncols * sizeof(double));
    double *rhs = (double *)malloc(nrows * sizeof(double));
    char *sense = (char *)malloc(nrows * sizeof(char));

    int *row_counts = (int *)calloc(nrows, sizeof(int));

    // Primeira passagem: contar não-zeros por linha
    for (int i = 0; i < ncols; i++) {
        AGEdge *e = &ag->edges[i];

        // Restrição de Fluxo (Ax = 0)
        row_counts[e->from_v]++;
        row_counts[e->to_v]++;

        // Restrições de Projeção (Cobertura)
        int p1 = agraph_project_e1(e);
        int p2 = agraph_project_e2(e);

        if (p1 != -1) row_counts[ag->n_vertices + p1]++;
        if (p2 != -1) row_counts[ag->n_vertices + g1->n_edges + p2]++;
    }

    // Preencher matbeg e rhs/sense
    matbeg[0] = 0;
    for (int i = 1; i < nrows; i++) {
        matbeg[i] = matbeg[i - 1] + row_counts[i - 1];
    }
    int nnz = matbeg[nrows - 1] + row_counts[nrows - 1];

    for (int i = 0; i < nrows; i++) {
        sense[i] = CPX_EQ; // Todas as restrições são igualdades (=)
        if (i < ag->n_vertices) rhs[i] = 0.0; // Fluxo = 0
        else rhs[i] = 1.0;                    // Cobertura = 1
    }

    // Segunda passagem: preencher matind e matval
    int *current_pos = (int *)malloc(nrows * sizeof(int));
    for (int i = 0; i < nrows; i++) current_pos[i] = matbeg[i];

    for (int i = 0; i < ncols; i++) {
        AGEdge *e = &ag->edges[i];

        // Fluxo: sai do vértice (-1), entra no vértice (+1)
        int r_from = e->from_v;
        int r_to = e->to_v;

        matind[current_pos[r_from]] = i;
        matval[current_pos[r_from]] = -1.0;
        current_pos[r_from]++;

        matind[current_pos[r_to]] = i;
        matval[current_pos[r_to]] = 1.0;
        current_pos[r_to]++;

        // Projeção
        int p1 = agraph_project_e1(e);
        int p2 = agraph_project_e2(e);

        if (p1 != -1) {
            int r = ag->n_vertices + p1;
            matind[current_pos[r]] = i;
            matval[current_pos[r]] = 1.0;
            current_pos[r]++;
        }
        if (p2 != -1) {
            int r = ag->n_vertices + g1->n_edges + p2;
            matind[current_pos[r]] = i;
            matval[current_pos[r]] = 1.0;
            current_pos[r]++;
        }
    }

    status = CPXaddrows(env, lp, 0, nrows, nnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if (status) fprintf(stderr, "Erro ao adicionar linhas.\n");

    /* --- 3. Resolver o PLI --- */
    status = CPXmipopt(env, lp);

    if (status == 0) {
        double objval;
        CPXgetobjval(env, lp, &objval);
        printf("  [CPLEX] Valor Otimo do PLI (Lower Bound): %.4f\n", objval);

        // Opcional: Imprimir as arestas selecionadas (x_e = 1)
        /*
        double *x_sol = (double *)malloc(ncols * sizeof(double));
        CPXgetx(env, lp, x_sol, 0, ncols - 1);
        printf("  Arestas selecionadas:\n");
        for (int i = 0; i < ncols; i++) {
            if (x_sol[i] > 0.5) {
                printf("    -> Aresta AG %d (cost %.1f)\n", i, ag->edges[i].cost);
            }
        }
        free(x_sol);
        */
    } else {
        printf("  [CPLEX] Falha ao resolver (Status: %d)\n", status);
    }

    /* --- 4. Limpeza --- */
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    free(obj); free(lb); free(ub); free(ctype);
    free(matbeg); free(matind); free(matval);
    free(rhs); free(sense);
    free(row_counts); free(current_pos);
}