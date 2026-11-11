#include "TAD_AGraph.h"
#include <stdio.h>
#include <stdlib.h>

/* ---------- indexação ---------- */

int agraph_vertex_index(const AGraph *ag, int u1, int u2)
{
    return u1 * ag->n2 + u2;
}

void agraph_vertex_coords(const AGraph *ag, int idx, int *u1, int *u2)
{
    *u1 = idx / ag->n2;
    *u2 = idx % ag->n2;
}

/* ---------- adição de aresta ---------- */

int agraph_add_edge(AGraph *ag, int from_v, int to_v,
                    AGEdgeKind kind, double cost,
                    int proj_e1, int proj_e2)
{
    if (!ag) return -1;
    if (ag->n_edges == ag->capacity) {
        size_t newcap = ag->capacity ? ag->capacity * 2 : 64;
        AGEdge *tmp = (AGEdge *)realloc(ag->edges, newcap * sizeof(AGEdge));
        if (!tmp) return -1;
        ag->edges    = tmp;
        ag->capacity = newcap;
    }
    AGEdge *e = &ag->edges[ag->n_edges++];
    e->from_v = from_v;
    e->to_v   = to_v;
    e->kind   = kind;
    e->cost   = cost;
    e->proj_e1 = proj_e1;
    e->proj_e2 = proj_e2;
    return 0;
}

/* ---------- construção completa ---------- */

AGraph *agraph_create(const UELGraph *g1, const UELGraph *g2)
{
    if (!g1 || !g2) return NULL;

    AGraph *ag = (AGraph *)calloc(1, sizeof(AGraph));
    if (!ag) return NULL;

    ag->g1         = g1;
    ag->g2         = g2;
    ag->n1         = g1->n_vertices;
    ag->n2         = g2->n_vertices;
    ag->n_vertices = ag->n1 * ag->n2;

    /* (1) arestas verticais:
     *     para cada (u1,v1) em E1 e cada u2 em V2          */
    for (size_t e1 = 0; e1 < g1->n_edges; e1++) {
        int u1 = g1->edges[e1].from;
        int v1 = g1->edges[e1].to;
        for (int u2 = 0; u2 < ag->n2; u2++) {
            int from = agraph_vertex_index(ag, u1, u2);
            int to   = agraph_vertex_index(ag, v1, u2);
            if (agraph_add_edge(ag, from, to, AG_EDGE_VERTICAL, 1.0,
                                (int)e1, -1) != 0) goto oom;
        }
    }

    /* (2) arestas horizontais:
     *     para cada u1 em V1 e cada (u2,v2) em E2          */
    for (int u1 = 0; u1 < ag->n1; u1++) {
        for (size_t e2 = 0; e2 < g2->n_edges; e2++) {
            int u2 = g2->edges[e2].from;
            int v2 = g2->edges[e2].to;
            int from = agraph_vertex_index(ag, u1, u2);
            int to   = agraph_vertex_index(ag, u1, v2);
            if (agraph_add_edge(ag, from, to, AG_EDGE_HORIZONTAL, 1.0,
                                -1, (int)e2) != 0) goto oom;
        }
    }

    /* (3) arestas diagonais:
     *     para cada (u1,v1) em E1 e cada (u2,v2) em E2.
     *     custo = 0 se l1 == l2 (match), 1 se mismatch.     */
    for (size_t e1 = 0; e1 < g1->n_edges; e1++) {
        int  u1 = g1->edges[e1].from;
        int  v1 = g1->edges[e1].to;
        char l1 = g1->edges[e1].label;
        for (size_t e2 = 0; e2 < g2->n_edges; e2++) {
            int  u2 = g2->edges[e2].from;
            int  v2 = g2->edges[e2].to;
            char l2 = g2->edges[e2].label;
            int from = agraph_vertex_index(ag, u1, u2);
            int to   = agraph_vertex_index(ag, v1, v2);
            double cost = (l1 == l2) ? 0.0 : 1.0;
            if (agraph_add_edge(ag, from, to, AG_EDGE_DIAGONAL, cost,
                                (int)e1, (int)e2) != 0) goto oom;
        }
    }
    return ag;

    oom:
    agraph_free(ag);
    return NULL;
}

void agraph_free(AGraph *ag)
{
    if (!ag) return;
    free(ag->edges);
    free(ag);
}

/* ---------- projeção π_i ---------- */

int agraph_project_e1(const AGEdge *e)
{
    return (e->kind == AG_EDGE_HORIZONTAL) ? -1 : e->proj_e1;
}

int agraph_project_e2(const AGEdge *e)
{
    return (e->kind == AG_EDGE_VERTICAL) ? -1 : e->proj_e2;
}

/* ---------- debug ---------- */

void agraph_print(const AGraph *ag)
{
    if (!ag) return;
    printf("AGraph: |V1|=%d, |V2|=%d, |V|=%d, |E|=%zu\n",
           ag->n1, ag->n2, ag->n_vertices, ag->n_edges);
    for (size_t i = 0; i < ag->n_edges; i++) {
        const AGEdge *e = &ag->edges[i];
        int u1, u2, v1, v2;
        agraph_vertex_coords(ag, e->from_v, &u1, &u2);
        agraph_vertex_coords(ag, e->to_v,   &v1, &v2);
        const char *k = (e->kind == AG_EDGE_VERTICAL)   ? "V"
                                                        : (e->kind == AG_EDGE_HORIZONTAL) ? "H" : "D";
        printf("  %3zu: [%s] (%d,%d) -> (%d,%d)  cost=%.1f  π1=%d  π2=%d\n",
               i, k, u1, u2, v1, v2, e->cost, e->proj_e1, e->proj_e2);
    }
}