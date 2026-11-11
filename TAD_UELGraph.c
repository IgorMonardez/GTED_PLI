//
// Created by Igor on 08/11/2026.
//

#include "TAD_UELGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINE_MAX 8192

/* Inicialização simples */
UELGraph *uelgraph_init(void) {
    UELGraph *g = (UELGraph *)calloc(1, sizeof(UELGraph));
    if (!g) return NULL;
    g->capacity = 16;
    g->edges = (UELEdge *)malloc(g->capacity * sizeof(UELEdge));
    if (!g->edges) {
        free(g);
        return NULL;
    }
    return g;
}

/* Liberação de memória, assumindo que edges e Sigma são owned pelo grafo */
void uelgraph_free(UELGraph *g) {
    free(g->edges);
    free(g->Sigma);
    g->edges = NULL;
    g->n_vertices = 0;
    g->n_edges = 0;
}

/* ---------- adição de aresta ---------- */
static int grow_edges(UELGraph *g, size_t need) {
    if (need <= g->capacity) return 0;
    size_t cap = g->capacity;
    while (cap < need) cap *= 2;
    UELEdge *tmp = (UELEdge *)realloc(g->edges, cap * sizeof(UELEdge));
    if (!tmp) return -1;
    g->edges = tmp;
    g->capacity = cap;
    return 0;
}

int uelgraph_add_edge(UELGraph *g, int from, int to, char label) {
    if (!g || from < 0 || to < 0) return -1;
    if (grow_edges(g, g->n_edges + 1) != 0) return -1;
    g->edges[g->n_edges].from  = from;
    g->edges[g->n_edges].to    = to;
    g->edges[g->n_edges].label = label;
    g->n_edges++;
    if (from >= g->n_vertices) g->n_vertices = from + 1;
    if (to   >= g->n_vertices) g->n_vertices = to   + 1;
    return 0;
}

int uelgraph_set_alphabet(UELGraph *g, const char *sigma) {
    if (!g) return -1;
    free(g->Sigma);
    g->Sigma = NULL;
    g->sigma_size = 0;
    if (!sigma) return 0;
    size_t n = strlen(sigma);
    g->Sigma = (char *)malloc(n + 1);
    if (!g->Sigma) return -1;
    memcpy(g->Sigma, sigma, n + 1);
    g->sigma_size = n;
    return 0;
}


/* ---------- parser GFA ---------- */
UELGraph *uelgraph_parse_gfa(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) return NULL;

    UELGraph *g = uelgraph_init();
    if (!g) {
        fclose(fp);
        return NULL;
    }

    bool seen[256] = {false};
    int max_vertex = -1;
    char line[LINE_MAX];

    while (fgets(line, sizeof(line), fp)) {
        size_t len = strlen(line);
        while (len && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';
        if (len == 0) continue;

        if (line[0] == 'S' && (line[1] == '\t' || line[1] == '\0')) {
            char *tab = strchr(line, '\t');
            if (!tab) continue;
            int sid = atoi(tab + 1);
            if (sid > max_vertex) max_vertex = sid;
        } else if (line[0] == 'L' && (line[1] == '\t' || line[1] == '\0')) {
            /* L <from> <forient> <to> <toorient> <overlap> <label> */
            char *f[7] = {0};
            int nf = 0;
            char *p = line;
            f[nf++] = p;
            for (; *p && nf < 7; p++) {
                if (*p == '\t') {
                    *p = '\0';
                    f[nf++] = p + 1;
                }
            }
            if (nf < 7) continue;
            int from = atoi(f[1]);
            int to = atoi(f[3]);
            char label = f[6][0];
            if (uelgraph_add_edge(g, from, to, label) != 0) {
                uelgraph_free(g);
                fclose(fp);
                return NULL;
            }
            if (from > max_vertex) max_vertex = from;
            if (to > max_vertex) max_vertex = to;
            seen[(unsigned char) label] = true;
        }
    }
    fclose(fp);

    if (max_vertex + 1 > g->n_vertices) g->n_vertices = max_vertex + 1;

    /* Sigma = rótulos efetivamente usados */
    size_t k = 0;
    for (int c = 0; c < 256; c++) if (seen[c]) k++;
    if (k) {
        char *sigma = (char *) malloc(k + 1);
        if (!sigma) {
            uelgraph_free(g);
            return NULL;
        }
        size_t idx = 0;
        for (int c = 0; c < 256; c++) if (seen[c]) sigma[idx++] = (char) c;
        sigma[idx] = '\0';
        if (uelgraph_set_alphabet(g, sigma) != 0) {
            free(sigma);
            uelgraph_free(g);
            return NULL;
        }
        free(sigma);
    }
    return g;
}

void uelgraph_print(const UELGraph *g) {
    if (!g) return;
    printf("UELGraph: n_vertices=%d, n_edges=%zu\n", g->n_vertices, g->n_edges);
    if (g->Sigma) printf("Sigma (%zu): %s\n", g->sigma_size, g->Sigma);
    for (size_t i = 0; i < g->n_edges; i++)
        printf("  %3zu: %3d -> %3d  [%c]\n",
               i, g->edges[i].from, g->edges[i].to, g->edges[i].label);
}

bool uelgraph_is_unidirectional(const UELGraph *g) {
    if (!g) return false;
    for (size_t i = 0; i < g->n_edges; i++)
        for (size_t j = i + 1; j < g->n_edges; j++)
            if (g->edges[i].from == g->edges[j].to &&
                g->edges[i].to   == g->edges[j].from)
                return false;
    return true;
}

bool uelgraph_has_eulerian_trail(const UELGraph *g, int *start, int *end) {
    if (!g) return false;
    int *indeg  = (int *)calloc(g->n_vertices, sizeof(int));
    int *outdeg = (int *)calloc(g->n_vertices, sizeof(int));
    if (!indeg || !outdeg) { free(indeg); free(outdeg); return false; }
    for (size_t i = 0; i < g->n_edges; i++) {
        outdeg[g->edges[i].from]++;
        indeg [g->edges[i].to  ]++;
    }
    int n_plus = 0, n_minus = 0, s = -1, t = -1;
    for (int v = 0; v < g->n_vertices; v++) {
        int d = outdeg[v] - indeg[v];
        if (d == 1)  { n_plus++;  s = v; }
        if (d == -1) { n_minus++; t = v; }
        if (d != 0 && d != 1 && d != -1) {
            free(indeg); free(outdeg); return false;
        }
    }
    bool ok = (n_plus == 0 && n_minus == 0) ||
              (n_plus == 1 && n_minus == 1);
    if (start) *start = s;
    if (end)   *end   = t;
    free(indeg); free(outdeg);
    return ok;
}
