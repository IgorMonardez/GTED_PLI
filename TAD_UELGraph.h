
#ifndef PLI_TAD_UELGRAPH_H
#define PLI_TAD_UELGRAPH_H
#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>

typedef struct {
    int from;    /* vértice de origem  */
    int to;      /* vértice de destino */
    char label;  /* rótulo de um único caractere */
} UELEdge;


typedef struct {
    int n_vertices;
    size_t n_edges;
    size_t capacity;
    UELEdge *edges;  /* vetor de tamanho n_edges: multiconjunto de arestas */
    char *Sigma;
    size_t sigma_size;

} UELGraph;

/* Inicialização simples */
UELGraph *uelgraph_init(void);

/* Liberação de memória, assumindo que edges e Sigma são owned pelo grafo */
void uelgraph_free(UELGraph *g);

/* modificadores */
int  uelgraph_add_edge(UELGraph *g, int from, int to, char label);
int  uelgraph_set_alphabet(UELGraph *g, const char *sigma);

/* I/O */
UELGraph *uelgraph_parse_gfa(const char *path);
void      uelgraph_print(const UELGraph *g);

/* validação de invariantes */
bool uelgraph_is_unidirectional(const UELGraph *g);
bool uelgraph_has_eulerian_trail(const UELGraph *g, int *start, int *end);

#endif //PLI_TAD_UELGRAPH_H
