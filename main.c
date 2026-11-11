#define _POSIX_C_SOURCE 200809L
#include "TAD_UELGraph.h"
#include "TAD_AGraph.h"
#include "ILP.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>

/* ------------------------------------------------------------------ */
/* Array dinâmico de grafos carregados                                */
/* ------------------------------------------------------------------ */

typedef struct {
    char     *path;
    UELGraph *graph;
} LoadedGraph;

typedef struct {
    LoadedGraph *items;
    size_t       count;
    size_t       capacity;
} GraphCollection;

static void gc_init(GraphCollection *gc)
{
    gc->items    = NULL;
    gc->count    = 0;
    gc->capacity = 0;
}

/* gc_push assume ownership de "g" e duplica "path".
 * Em caso de falha de alocação, "g" NÃO é liberado (fica a cargo do caller). */
static int gc_push(GraphCollection *gc, const char *path, UELGraph *g)
{
    if (gc->count == gc->capacity) {
        size_t newcap = gc->capacity ? gc->capacity * 2 : 16;
        LoadedGraph *tmp = (LoadedGraph *)realloc(
                gc->items, newcap * sizeof(LoadedGraph));
        if (!tmp) return -1;
        gc->items    = tmp;
        gc->capacity = newcap;
    }
    char *p = strdup(path);
    if (!p) return -1;
    gc->items[gc->count].path  = p;
    gc->items[gc->count].graph = g;
    gc->count++;
    return 0;
}

static void gc_free(GraphCollection *gc)
{
    for (size_t i = 0; i < gc->count; i++) {
        free(gc->items[i].path);
        uelgraph_free(gc->items[i].graph);
    }
    free(gc->items);
    gc->items    = NULL;
    gc->count    = 0;
    gc->capacity = 0;
}

/* ------------------------------------------------------------------ */
/* Carregamento (arquivo único e diretório)                           */
/* ------------------------------------------------------------------ */

static int load_file(const char *path, GraphCollection *gc)
{
    UELGraph *g = uelgraph_parse_gfa(path);
    if (!g) {
        fprintf(stderr, "[FALHA] %s\n", path);
        return -1;
    }
    if (gc_push(gc, path, g) != 0) {
        uelgraph_free(g);
        return -1;
    }
    return 0;
}

static int cmp_str(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static int load_directory(const char *dir, const char *ext, GraphCollection *gc)
{
    DIR *d = opendir(dir);
    if (!d) { perror("opendir"); return -1; }

    size_t extlen = ext ? strlen(ext) : 0;

    char **names = NULL;
    size_t n = 0, cap = 0;
    struct dirent *ent;

    while ((ent = readdir(d)) != NULL) {
        if (ent->d_name[0] == '.') continue;
        if (ext) {
            size_t l = strlen(ent->d_name);
            if (l < extlen ||
                strcmp(ent->d_name + l - extlen, ext) != 0) continue;
        }
        char full[4096];
        snprintf(full, sizeof(full), "%s/%s", dir, ent->d_name);
        struct stat st;
        if (stat(full, &st) != 0 || !S_ISREG(st.st_mode)) continue;

        if (n == cap) {
            size_t newcap = cap ? cap * 2 : 32;
            char **tmp = (char **)realloc(names, newcap * sizeof(char *));
            if (!tmp) {
                for (size_t k = 0; k < n; k++) free(names[k]);
                free(names);
                closedir(d);
                return -1;
            }
            names = tmp;
            cap   = newcap;
        }
        names[n] = strdup(ent->d_name);
        if (!names[n]) {
            for (size_t k = 0; k < n; k++) free(names[k]);
            free(names);
            closedir(d);
            return -1;
        }
        n++;
    }
    closedir(d);

    if (n == 0) {
        fprintf(stderr, "nenhum arquivo%s encontrado em %s\n",
                ext ? ext : "", dir);
        free(names);
        return 0;
    }
    qsort(names, n, sizeof(char *), cmp_str);

    int errors = 0;
    for (size_t i = 0; i < n; i++) {
        char full[4096];
        snprintf(full, sizeof(full), "%s/%s", dir, names[i]);
        if (load_file(full, gc) != 0) errors++;
        free(names[i]);
    }
    free(names);

    printf("carregados %zu arquivo(s), %d erro(s)\n", n, errors);
    return errors ? -1 : 0;
}

/* ------------------------------------------------------------------ */
/* Uso do array: imprimir todos os grafos carregados                  */
/* ------------------------------------------------------------------ */

static void print_collection(const GraphCollection *gc)
{
    printf("Total de grafos carregados: %zu\n\n", gc->count);
    for (size_t i = 0; i < gc->count; i++) {
        printf("=== [%zu] %s ===\n", i, gc->items[i].path);
        uelgraph_print(gc->items[i].graph);
        printf("\n");
    }
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */


int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "uso: %s <arquivo|pasta> [extensao .txt]\n", argv[0]);
        return 1;
    }

    struct stat st;
    if (stat(argv[1], &st) != 0) { perror("stat"); return 1; }

    GraphCollection gc;
    gc_init(&gc);

    int rc;
    if (S_ISDIR(st.st_mode)) {
        const char *ext = (argc >= 3) ? argv[2] : ".txt";
        rc = load_directory(argv[1], ext, &gc);
    } else {
        rc = load_file(argv[1], &gc);
    }

    if (rc != 0)
        fprintf(stderr, "aviso: houve falha(s) durante o carregamento\n");

    print_collection(&gc);

    for(size_t i = 0; i < gc.count; i++) {
        for(size_t j = 0; j < gc.count; j++) {
            if(i != j) {
                printf("--- AG entre G1: %s e G2: %s ---\n", gc.items[i].path, gc.items[j].path);
                AGraph *ag = agraph_create(gc.items[i].graph, gc.items[j].graph);

                if (ag) {
                    // agraph_print(ag); // Descomente se quiser ver o grafo de alinhamento

                    // Chama o resolvedor do CPLEX
                    solve_ilp_cplex(ag, gc.items[i].graph, gc.items[j].graph);

                    agraph_free(ag); // <--- IMPORTANTE: Evita vazamento de memória
                }
            }
        }
    }

    gc_free(&gc);
    return rc == 0 ? 0 : 2;
}