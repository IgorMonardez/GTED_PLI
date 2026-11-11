//
// Created by Igor on 09/11/2026.
//

#ifndef PLI_TAD_AGRAPH_H
#define PLI_TAD_AGRAPH_H

#include "TAD_UELGraph.h"
#include <stddef.h>

/* Tipos de aresta do grafo de alinhamento (Definição 2). */
typedef enum {
    AG_EDGE_VERTICAL,    /* [(u1,u2),(v1,u2)]  -> projeta em G1          */
    AG_EDGE_HORIZONTAL,  /* [(u1,u2),(u1,v2)]  -> projeta em G2          */
    AG_EDGE_DIAGONAL     /* [(u1,u2),(v1,v2)]  -> projeta em G1 e G2     */
} AGEdgeKind;

/* Aresta de A(G1,G2). Como E é multiconjunto, arestas paralelas
 * são entradas distintas no vetor.                            */
typedef struct {
    int        from_v;   /* índice plano (u1*n2 + u2) */
    int        to_v;
    AGEdgeKind kind;
    double     cost;     /* δ(e)                       */
    int        proj_e1;  /* índice em g1->edges ou -1  */
    int        proj_e2;  /* índice em g2->edges ou -1  */
} AGEdge;

/* Grafo de alinhamento A(G1,G2) = (V, E, δ). */
typedef struct {
    int     n1, n2;      /* |V1|, |V2|                       */
    int     n_vertices;  /* |V| = n1 * n2                    */
    size_t  n_edges;
    size_t  capacity;
    AGEdge *edges;

    const UELGraph *g1;  /* referência NÃO-owning            */
    const UELGraph *g2;
} AGraph;

/* ---- ciclo de vida ---- */
AGraph *agraph_create(const UELGraph *g1, const UELGraph *g2);
void    agraph_free(AGraph *ag);

/* ---- indexação de vértices (V1 x V2 <-> índice plano) ---- */
int     agraph_vertex_index(const AGraph *ag, int u1, int u2);
void    agraph_vertex_coords(const AGraph *ag, int idx, int *u1, int *u2);

/* ---- construção de baixo nível ---- */
int     agraph_add_edge(AGraph *ag, int from_v, int to_v,
                        AGEdgeKind kind, double cost,
                        int proj_e1, int proj_e2);

/* ---- projeção π_i (retorna -1 quando não se aplica) ---- */
int     agraph_project_e1(const AGEdge *e);
int     agraph_project_e2(const AGEdge *e);

/* ---- util ---- */
void    agraph_print(const AGraph *ag);

#endif //PLI_TAD_AGRAPH_H
