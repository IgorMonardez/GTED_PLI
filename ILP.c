#include "ILP.h"
#include <ilcplex/cplex.h>
#include <stdio.h>
#include <stdlib.h>

// Estrutura temporária para ajudar a ordenar os elementos da matriz
typedef struct {
    int row;
    int col;
    double val;
} Triple;

// Função de comparação para o qsort (ordena por linha, depois por coluna)
static int compare_triples(const void *a, const void *b) {
    Triple *ta = (Triple *)a;
    Triple *tb = (Triple *)b;
    if (ta->row != tb->row) return ta->row - tb->row;
    return ta->col - tb->col;
}

void solve_ilp_cplex(const AGraph *ag, const UELGraph *g1, const UELGraph *g2) {
    if (!ag || !g1 || !g2) return;

    if (ag->n_edges == 0) {
        printf("  [CPLEX] Grafo de alinhamento sem arestas. Custo = 0.0\n");
        return;
    }

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
        obj[i] = ag->edges[i].cost;
        lb[i]  = 0.0;
        ub[i]  = 1.0;
        ctype[i] = CPX_BINARY; // Variável binária (0 ou 1)
    }

    status = CPXnewcols(env, lp, ncols, obj, lb, ub, ctype, NULL);
    if (status) fprintf(stderr, "Erro ao adicionar colunas.\n");

    /* --- 2. Construir Matriz de Restrições --- */
    // Alocamos espaço para no máximo 4 não-zeros por coluna (from, to, proj1, proj2)
    Triple *triples = (Triple *)malloc(4 * ncols * sizeof(Triple));
    int t_idx = 0;

    for (int i = 0; i < ncols; i++) {
        AGEdge *e = &ag->edges[i];

        // Restrição de Fluxo (Ax = 0)
        if (e->from_v == e->to_v) {
            // Self-loop: -1 + 1 = 0, então não adicionamos nada (evita problemas)
        } else {
            triples[t_idx].row = e->from_v; triples[t_idx].col = i; triples[t_idx].val = -1.0; t_idx++;
            triples[t_idx].row = e->to_v;   triples[t_idx].col = i; triples[t_idx].val =  1.0; t_idx++;
        }

        // Restrições de Projeção (Cobertura = 1)
        int p1 = agraph_project_e1(e);
        int p2 = agraph_project_e2(e);

        if (p1 != -1) {
            triples[t_idx].row = ag->n_vertices + p1; triples[t_idx].col = i; triples[t_idx].val = 1.0; t_idx++;
        }
        if (p2 != -1) {
            triples[t_idx].row = ag->n_vertices + g1->n_edges + p2; triples[t_idx].col = i; triples[t_idx].val = 1.0; t_idx++;
        }
    }

    int nnz = t_idx;

    // Ordenar os triplos por linha e depois por coluna (Requisito do CPLEX)
    qsort(triples, nnz, sizeof(Triple), compare_triples);

    // Alocar arrays usando CPXINT (garante compatibilidade com a arquitetura do CPLEX)
    CPXINT *matbeg = (CPXINT *)malloc(nrows * sizeof(CPXINT));
    CPXINT *matind = (CPXINT *)malloc(nnz * sizeof(CPXINT));
    double *matval = (double *)malloc(nnz * sizeof(double));
    double *rhs = (double *)malloc(nrows * sizeof(double));
    char *sense = (char *)malloc(nrows * sizeof(char));

    // Preencher matbeg e a matriz (matind e matval)
    int t = 0;
    for (int r = 0; r < nrows; r++) {
        matbeg[r] = t;
        while (t < nnz && triples[t].row == r) {
            matind[t] = triples[t].col;
            matval[t] = triples[t].val;
            t++;
        }
    }

    for (int i = 0; i < nrows; i++) {
        sense[i] = 'E'; // 'E' = Equality (Igualdade)
        if (i < ag->n_vertices) rhs[i] = 0.0; // Fluxo = 0
        else rhs[i] = 1.0;                    // Cobertura = 1
    }

    status = CPXaddrows(env, lp, 0, nrows, nnz, rhs, sense, matbeg, matind, matval, NULL, NULL);
    if (status) fprintf(stderr, "Erro ao adicionar linhas.\n");

    /* --- 3. Resolver o PLI --- */
    // Desabilita os logs do CPLEX para não poluir o terminal
    CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);

    status = CPXmipopt(env, lp);

    if (status == 0) {
        double objval;
        CPXgetobjval(env, lp, &objval);
        printf("  [CPLEX] Valor Otimo do PLI (Lower Bound): %.4f\n", objval);
    } else {
        printf("  [CPLEX] Falha ao resolver (Status: %d)\n", status);
    }

    /* --- 4. Limpeza --- */
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env);

    free(obj); free(lb); free(ub); free(ctype);
    free(matbeg); free(matind); free(matval);
    free(rhs); free(sense);
    free(triples);
}