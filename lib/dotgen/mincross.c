/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: Details at https://graphviz.org
 *************************************************************************/


/* 
 * dot_mincross(g) takes a ranked graphs, and finds an ordering
 * that avoids edge crossings.  clusters are expanded.
 * N.B. the rank structure is global (not allocated per cluster)
 * because mincross may compare nodes in different clusters.
 */

#include <assert.h>
#include <cgraph/agxbuf.h>
#include <cgraph/alloc.h>
#include <cgraph/cgraph.h>
#include <cgraph/exit.h>
#include <dotgen/dot.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

struct adjmatrix_t {
  size_t nrows;
  size_t ncols;
  char *data;
};

/* #define DEBUG */
#define MARK(v)		(ND_mark(v))
#define saveorder(v)	(ND_coord(v)).x
#define flatindex(v)	((size_t)ND_low(v))

	/* forward declarations */
static bool medians(graph_t * g, int r0, int r1);
static int get_edge_panelty_with_type(edge_t *e, int type);
static int get_edge_penalty(edge_t *e);
static void rename_virtual_nodes(Agraph_t *g);
static void print_ranks(graph_t * g);
static char* get_name(void * obj);
static int get_type(void * obj);
static void examine_order(graph_t * g);
static bool is_left_subgraph(node_t *n);
static bool is_right_subgraph(node_t *n);
static bool is_context(node_t *n);
static int nodeposcmpf(node_t ** n0, node_t ** n1);
static int edgeidcmpf(edge_t ** e0, edge_t ** e1);
static void flat_breakcycles(graph_t * g);
static void flat_reorder(graph_t * g);
static void flat_search(graph_t * g, node_t * v);
static void init_mincross(graph_t * g);
static void merge2(graph_t * g);
static void init_mccomp(graph_t *g, size_t c);
static void cleanup2(graph_t * g, int nc);
static int mincross_clust(graph_t * g, int);
static int mincross(graph_t * g, int startpass, int endpass, int);
static void mincross_step(graph_t * g, int pass);
static void mincross_options(graph_t * g);
static void save_best(graph_t * g);
static void restore_best(graph_t * g);
static adjmatrix_t *new_matrix(size_t i, size_t j);
static void free_matrix(adjmatrix_t * p);
static int ordercmpf(int *i0, int *i1);
#ifdef DEBUG
#if DEBUG > 1
static int gd_minrank(Agraph_t *g) {return GD_minrank(g);}
static int gd_maxrank(Agraph_t *g) {return GD_maxrank(g);}
static rank_t *gd_rank(Agraph_t *g, int r) {return &GD_rank(g)[r];}
static int nd_order(Agnode_t *v) { return ND_order(v); }
#endif
void check_rs(graph_t * g, int null_ok);
void check_order(void);
void check_vlists(graph_t * g);
void node_in_root_vlist(node_t * n);
#endif


	/* mincross parameters */
static int MinQuit;
static double Convergence;

static graph_t *Root;
static int GlobalMinRank, GlobalMaxRank;
static edge_t **TE_list;
static int *TI_list;
static bool ReMincross;

#if defined(DEBUG) && DEBUG > 1
static void indent(graph_t* g)
{
  if (g->parent) {
    fprintf (stderr, "  ");
    indent(g->parent);
  }
}

static char* nname(node_t* v)
{
        static char buf[1000];
	if (ND_node_type(v)) {
		if (ND_ranktype(v) == CLUSTER)
			snprintf(buf, sizeof(buf), "v%s_%p", agnameof(ND_clust(v)), v);
		else
			snprintf(buf, sizeof(buf), "v_%p", v);
	} else
		snprintf(buf, sizeof(buf), "%s", agnameof(v));
	return buf;
}
static void dumpg (graph_t* g)
{
    int j, i, r;
    node_t* v;
    edge_t* e;

    fprintf (stderr, "digraph A {\n");
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	fprintf (stderr, "  subgraph {rank=same  ");
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	  v = GD_rank(g)[r].v[i];
          if (i > 0)
 	    fprintf (stderr, " -> %s", nname(v));
          else
 	    fprintf (stderr, "%s", nname(v));
        }
        if (i > 1) fprintf (stderr, " [style=invis]}\n");
        else fprintf (stderr, " }\n");
    }
    for (r = GD_minrank(g); r < GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	  v = GD_rank(g)[r].v[i];
	  for (j = 0; (e = ND_out(v).list[j]); j++) {
             fprintf (stderr, "%s -> ", nname(v));
             fprintf (stderr, "%s\n", nname(aghead(e)));
          }
        }
    }
    fprintf (stderr, "}\n");
}
static void dumpr (graph_t* g, int edges)
{
    int j, i, r;
    node_t* v;
    edge_t* e;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	fprintf (stderr, "[%d] ", r);
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	  v = GD_rank(g)[r].v[i];
 	  fprintf (stderr, "%s(%.02f,%d) ", nname(v), saveorder(v),ND_order(v));
        }
	fprintf (stderr, "\n");
    }
    if (edges == 0) return;
    for (r = GD_minrank(g); r < GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	  v = GD_rank(g)[r].v[i];
	  for (j = 0; (e = ND_out(v).list[j]); j++) {
             fprintf (stderr, "%s -> ", nname(v));
             fprintf (stderr, "%s\n", nname(aghead(e)));
          }
        }
    }
}
#endif

typedef struct {
    Agrec_t h;
    int x, lo, hi;
    Agnode_t* np;
} info_t;

#define ND_x(n) (((info_t*)AGDATA(n))->x)
#define ND_lo(n) (((info_t*)AGDATA(n))->lo)
#define ND_hi(n) (((info_t*)AGDATA(n))->hi)
#define ND_np(n) (((info_t*)AGDATA(n))->np)
#define ND_idx(n) (ND_order(ND_np(n)))

static void
emptyComp (graph_t* sg)
{
    Agnode_t* n;
    Agnode_t* nxt;

    for (n = agfstnode(sg); n; n = nxt) {
	nxt = agnxtnode (sg, n);
	agdelnode(sg,n);
    }
}

#define isBackedge(e) (ND_idx(aghead(e)) > ND_idx(agtail(e)))

static Agnode_t*
findSource (Agraph_t* g, Agraph_t* sg)
{
    Agnode_t* n;

    for (n = agfstnode(sg); n; n = agnxtnode(sg, n))
	if (agdegree(g,n,1,0) == 0) return n;
    return NULL;
}

static int
topsort (Agraph_t* g, Agraph_t* sg, Agnode_t** arr)
{
    Agnode_t* n;
    Agedge_t* e;
    Agedge_t* nxte;
    int cnt = 0;

    while ((n = findSource(g, sg))) {
	arr[cnt++] = ND_np(n); 
	agdelnode(sg, n);
	for (e = agfstout(g, n); e; e = nxte) {
	    nxte = agnxtout(g, e); 
	    agdeledge(g, e);
	}
    }
    return cnt;
}

static int
getComp (graph_t* g, node_t* n, graph_t* comp, int* indices)
{
    int backedge = 0;
    Agedge_t* e;

    ND_x(n) = 1;
    indices[agnnodes(comp)] = ND_idx(n);
    agsubnode(comp, n, 1);
    for (e = agfstout(g,n); e; e = agnxtout(g,e)) {
	if (isBackedge(e)) backedge++;
	if (!ND_x(aghead(e)))
	    backedge += getComp(g, aghead(e), comp, indices);
    }
    for (e = agfstin(g,n); e; e = agnxtin(g,e)) {
	if (isBackedge(e)) backedge++;
	if (!ND_x(agtail(e)))
	    backedge += getComp(g, agtail(e), comp, indices);
    }
    return backedge;
}

/* fixLabelOrder:
 * For each pair of nodes (labels), we add an edge 
 */
static void
fixLabelOrder (graph_t* g, rank_t* rk)
{
    int cnt;
    bool haveBackedge = false;
    Agraph_t* sg;
    Agnode_t* n;
    Agnode_t* nxtp;
    Agnode_t* v;

    for (n = agfstnode(g); n; n = nxtp) {
	v = nxtp = agnxtnode(g, n);
	for (; v; v = agnxtnode(g, v)) {
	    if (ND_hi(v) <= ND_lo(n)) { 
		haveBackedge = true;
		agedge(g, v, n, NULL, 1);
	    }
	    else if (ND_hi(n) <= ND_lo(v)) { 
		agedge(g, n, v, NULL, 1);
	    }
	}
    }
    if (!haveBackedge) return;
    
    sg = agsubg(g, "comp", 1);
    Agnode_t **arr = gv_calloc(agnnodes(g), sizeof(Agnode_t*));
    int *indices = gv_calloc(agnnodes(g), sizeof(int));

    for (n = agfstnode(g); n; n = agnxtnode(g,n)) {
	if (ND_x(n) || agdegree(g,n,1,1) == 0) continue;
	if (getComp(g, n, sg, indices)) {
	    int i, sz = agnnodes(sg);
	    cnt = topsort (g, sg, arr);
	    assert (cnt == sz);
	    qsort(indices, cnt, sizeof(int), (qsort_cmpf)ordercmpf);
	    for (i = 0; i < sz; i++) {
		ND_order(arr[i]) = indices[i];
		rk->v[indices[i]] = arr[i];
	    }
       }
       emptyComp(sg);
    }
    free (arr);
}

/* checkLabelOrder:
 * Check that the ordering of labels for flat edges is consistent. 
 * This is necessary because dot_position will attempt to force the label
 * to be between the edge's vertices. This can lead to an infeasible problem.
 *
 * We check each rank for any flat edge labels (as dummy nodes) and create a
 * graph with a node for each label. If the graph contains more than 1 node, we
 * call fixLabelOrder to see if there really is a problem and, if so, fix it.
 */ 
void
checkLabelOrder (graph_t* g)
{
    int j, r, lo, hi;
    graph_t* lg = NULL;
    agxbuf buf = {0};
    rank_t* rk;
    Agnode_t* u;
    Agnode_t* n;
    Agedge_t* e;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	rk = GD_rank(g)+r;
	for (j = 0; j < rk->n; j++) {
	    u = rk->v[j];
	    if ((e = ND_alg(u))) {
		if (!lg) lg = agopen ("lg", Agstrictdirected, 0);
		agxbprint(&buf, "%d", j);
		n = agnode(lg, agxbuse(&buf), 1);
		agbindrec(n, "info", sizeof(info_t), true);
		lo = ND_order(aghead(ND_out(u).list[0]));
		hi = ND_order(aghead(ND_out(u).list[1]));
		if (lo > hi) {
		    int tmp;
		    tmp = lo;
		    lo = hi;
		    hi = tmp;
		}
		ND_lo(n) = lo;
		ND_hi(n) = hi;
		ND_np(n) = u;
	    }
	}
	if (lg) {
	    if (agnnodes(lg) > 1) fixLabelOrder (lg, rk);
	    agclose(lg);
	    lg = NULL;
	}
    }
    agxbfree(&buf);
}

/* dot_mincross:
 * Minimize edge crossings
 * Note that nodes are not placed into GD_rank(g) until mincross()
 * is called.
 */
void dot_mincross(graph_t * g, int doBalance)
{
    int nc;
    char *s;

    /* check whether malformed input has led to empty cluster that the crossing
     * functions will not anticipate
     */
    {
	size_t i;
	for (i = 1; i <= (size_t)GD_n_cluster(g); ) {
	    if (agfstnode(GD_clust(g)[i]) == NULL) {
	      agwarningf("removing empty cluster\n");
	      memmove(&GD_clust(g)[i], &GD_clust(g)[i + 1],
	        ((size_t)GD_n_cluster(g) - i) * sizeof(GD_clust(g)[0]));
	      --GD_n_cluster(g);
	    } else {
	      ++i;
	    }
	}
    }

    init_mincross(g);

	// if (GD_comp(g).size > 1) {
	// 	// Do not support multiple components, print and return
	// 	fprintf(stderr, "[sy]  - Do not support multiple components\n");
	// 	return;
	// }

    size_t comp;
    for (nc = 0, comp = 0; comp < GD_comp(g).size; comp++) {
	init_mccomp(g, comp);
	nc += mincross(g, 0, 2, doBalance);
    }

    merge2(g);

    /* run mincross on contents of each cluster */
    for (int c = 1; c <= GD_n_cluster(g); c++) {
	nc += mincross_clust(GD_clust(g)[c], doBalance);
#ifdef DEBUG
	check_vlists(GD_clust(g)[c]);
	check_order();
#endif
    }

    if (GD_n_cluster(g) > 0 && (!(s = agget(g, "remincross")) || mapbool(s))) {
	mark_lowclusters(g);
	ReMincross = true;
	nc = mincross(g, 2, 2, doBalance);
#ifdef DEBUG
	for (int c = 1; c <= GD_n_cluster(g); c++)
	    check_vlists(GD_clust(g)[c]);
#endif
    }
    cleanup2(g, nc);
}

static adjmatrix_t *new_matrix(size_t i, size_t j) {
    adjmatrix_t *rv = gv_alloc(sizeof(adjmatrix_t));
    rv->nrows = i;
    rv->ncols = j;
    rv->data = gv_calloc(i * j, sizeof(char));
    return rv;
}

static void free_matrix(adjmatrix_t * p)
{
    if (p) {
	free(p->data);
	free(p);
    }
}

#define ELT(M,i,j)		(M->data[((i)*M->ncols)+(j)])

static void init_mccomp(graph_t *g, size_t c) {
    int r;

    GD_nlist(g) = GD_comp(g).list[c];
    if (c > 0) {
	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	    GD_rank(g)[r].v = GD_rank(g)[r].v + GD_rank(g)[r].n;
	    GD_rank(g)[r].n = 0;
	}
    }
}

static int betweenclust(edge_t * e)
{
    while (ED_to_orig(e))
	e = ED_to_orig(e);
    return (ND_clust(agtail(e)) != ND_clust(aghead(e)));
}

static void do_ordering_node(graph_t *g, node_t *n, bool outflag) {
    int i, ne;
    node_t *u, *v;
    edge_t *e, *f, *fe;
    edge_t **sortlist = TE_list;

    if (ND_clust(n))
	return;
    if (outflag) {
	for (i = ne = 0; (e = ND_out(n).list[i]); i++)
	    if (!betweenclust(e))
		sortlist[ne++] = e;
    } else {
	for (i = ne = 0; (e = ND_in(n).list[i]); i++)
	    if (!betweenclust(e))
		sortlist[ne++] = e;
    }
    if (ne <= 1)
	return;
    /* write null terminator at end of list.
       requires +1 in TE_list alloccation */
    sortlist[ne] = 0;
    qsort(sortlist, ne, sizeof(sortlist[0]), (qsort_cmpf) edgeidcmpf);
    for (ne = 1; (f = sortlist[ne]); ne++) {
	e = sortlist[ne - 1];
	if (outflag) {
	    u = aghead(e);
	    v = aghead(f);
	} else {
	    u = agtail(e);
	    v = agtail(f);
	}
	if (find_flat_edge(u, v))
	    return;
	fe = new_virtual_edge(u, v, NULL);
	ED_edge_type(fe) = FLATORDER;
	flat_edge(g, fe);
    }
}

static void do_ordering(graph_t *g, bool outflag) {
    /* Order all nodes in graph */
    node_t *n;

    for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
	do_ordering_node (g, n, outflag);
    }
}

static void do_ordering_for_nodes(graph_t * g)
{
    /* Order nodes which have the "ordered" attribute */
    node_t *n;
    const char *ordering;

    for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
	if ((ordering = late_string(n, N_ordering, NULL))) {
	    if (streq(ordering, "out"))
		do_ordering_node(g, n, true);
	    else if (streq(ordering, "in"))
		do_ordering_node(g, n, false);
	    else if (ordering[0])
		agerr(AGERR, "ordering '%s' not recognized for node '%s'.\n", ordering, agnameof(n));
	}
    }
}

/* ordered_edges:
 * handle case where graph specifies edge ordering
 * If the graph does not have an ordering attribute, we then
 * check for nodes having the attribute.
 * Note that, in this implementation, the value of G_ordering
 * dominates the value of N_ordering.
 */
static void ordered_edges(graph_t * g)
{
    char *ordering;

    if (!G_ordering && !N_ordering)
	return;
    if ((ordering = late_string(g, G_ordering, NULL))) {
	if (streq(ordering, "out"))
	    do_ordering(g, true);
	else if (streq(ordering, "in"))
	    do_ordering(g, false);
	else if (ordering[0])
	    agerr(AGERR, "ordering '%s' not recognized.\n", ordering);
    }
    else
    {
	graph_t *subg;

	for (subg = agfstsubg(g); subg; subg = agnxtsubg(subg)) {
	    /* clusters are processed by separate calls to ordered_edges */
	    if (!is_cluster(subg))
		ordered_edges(subg);
	}
	if (N_ordering) do_ordering_for_nodes (g);
    }
}

static int mincross_clust(graph_t * g, int doBalance)
{
    int c, nc;

    expand_cluster(g);
    ordered_edges(g);
    flat_breakcycles(g);
    flat_reorder(g);
    nc = mincross(g, 2, 2, doBalance);

    for (c = 1; c <= GD_n_cluster(g); c++)
	nc += mincross_clust(GD_clust(g)[c], doBalance);

    save_vlist(g);
    return nc;
}

static bool left2right(graph_t *g, node_t *v, node_t *w) {
    adjmatrix_t *M;

    /* CLUSTER indicates orig nodes of clusters, and vnodes of skeletons */
    if (!ReMincross) {
	if (ND_clust(v) != ND_clust(w) && ND_clust(v) && ND_clust(w)) {
	    /* the following allows cluster skeletons to be swapped */
	    if (ND_ranktype(v) == CLUSTER && ND_node_type(v) == VIRTUAL)
		return false;
	    if (ND_ranktype(w) == CLUSTER && ND_node_type(w) == VIRTUAL)
		return false;
	    return true;
	}
    } else {
	if (ND_clust(v) != ND_clust(w))
	    return true;
    }
    M = GD_rank(g)[ND_rank(v)].flat;
    if (M == NULL)
	return false;
    if (GD_flip(g)) {
	node_t *t = v;
	v = w;
	w = t;
    }
    return ELT(M, flatindex(v), flatindex(w)) != 0;
}

static int in_cross(node_t * v, node_t * w)
{
    edge_t **e1, **e2;
    int inv, cross = 0, t;

    for (e2 = ND_in(w).list; *e2; e2++) {
	int cnt = get_edge_penalty(*e2);		
		
	inv = ND_order(agtail(*e2));

	for (e1 = ND_in(v).list; *e1; e1++) {
	    t = ND_order(agtail(*e1)) - inv;
	    if (t > 0 || (t == 0 && ED_tail_port(*e1).p.x > ED_tail_port(*e2).p.x))
		cross += get_edge_penalty(*e1) * cnt;
	}
    }
    return cross;
}

static int out_cross(node_t * v, node_t * w)
{
    edge_t **e1, **e2;
    int inv, cross = 0, t;

    for (e2 = ND_out(w).list; *e2; e2++) {
	int cnt = get_edge_penalty(*e2);
	inv = ND_order(aghead(*e2));

	for (e1 = ND_out(v).list; *e1; e1++) {
	    t = ND_order(aghead(*e1)) - inv;
	    if (t > 0 || (t == 0 && (ED_head_port(*e1)).p.x > (ED_head_port(*e2)).p.x))
		cross += get_edge_penalty(*e1) * cnt;
	}
    }
    return cross;

}

static void exchange(node_t * v, node_t * w)
{
    int vi, wi, r;

	fprintf(stderr, "[sy]    - exchange node (%s, %s)\n", agnameof(v), agnameof(w));

    r = ND_rank(v);
    vi = ND_order(v);
    wi = ND_order(w);
    ND_order(v) = wi;
    GD_rank(Root)[r].v[wi] = v;
    ND_order(w) = vi;
    GD_rank(Root)[r].v[vi] = w;
}

static void balanceNodes(graph_t * g, int r, node_t * v, node_t * w)
{
    node_t *s;			/* separator node */
    int sepIndex = 0;
    int nullType;		/* type of null nodes */
    int cntDummy = 0, cntOri = 0;
    int k = 0, m = 0, k1 = 0, m1 = 0, i = 0;

    /* we only consider v and w of different types */
    if (ND_node_type(v) == ND_node_type(w))
	return;

    /* count the number of dummy and original nodes */
    for (i = 0; i < GD_rank(g)[r].n; i++) {
	if (ND_node_type(GD_rank(g)[r].v[i]) == NORMAL)
	    cntOri++;
	else
	    cntDummy++;
    }

    if (cntOri < cntDummy) {
	if (ND_node_type(v) == NORMAL)
	    s = v;
	else
	    s = w;
    } else {
	if (ND_node_type(v) == NORMAL)
	    s = w;
	else
	    s = v;
    }

    /* get the separator node index */
    for (i = 0; i < GD_rank(g)[r].n; i++) {
	if (GD_rank(g)[r].v[i] == s)
	    sepIndex = i;
    }

    nullType = ND_node_type(s) == NORMAL ? VIRTUAL : NORMAL;

    /* count the number of null nodes to the left and 
     * right of the separator node 
     */
    for (i = sepIndex - 1; i >= 0; i--) {
	if (ND_node_type(GD_rank(g)[r].v[i]) == nullType)
	    k++;
	else
	    break;
    }

    for (i = sepIndex + 1; i < GD_rank(g)[r].n; i++) {
	if (ND_node_type(GD_rank(g)[r].v[i]) == nullType)
	    m++;
	else
	    break;
    }

    /* now exchange v,w and calculate the same counts */

    exchange(v, w);

    /* get the separator node index */
    for (i = 0; i < GD_rank(g)[r].n; i++) {
	if (GD_rank(g)[r].v[i] == s)
	    sepIndex = i;
    }

    /* count the number of null nodes to the left and 
     * right of the separator node 
     */
    for (i = sepIndex - 1; i >= 0; i--) {
	if (ND_node_type(GD_rank(g)[r].v[i]) == nullType)
	    k1++;
	else
	    break;
    }

    for (i = sepIndex + 1; i < GD_rank(g)[r].n; i++) {
	if (ND_node_type(GD_rank(g)[r].v[i]) == nullType)
	    m1++;
	else
	    break;
    }

    if (abs(k1 - m1) > abs(k - m)) {
	exchange(v, w);		//revert to the original ordering
    }
}

static int balance(graph_t * g)
{
    int i, c0, c1, rv;
    node_t *v, *w;
    int r;

    rv = 0;

    for (r = GD_maxrank(g); r >= GD_minrank(g); r--) {

	GD_rank(g)[r].candidate = false;
	for (i = 0; i < GD_rank(g)[r].n - 1; i++) {
	    v = GD_rank(g)[r].v[i];
	    w = GD_rank(g)[r].v[i + 1];
	    assert(ND_order(v) < ND_order(w));
	    if (left2right(g, v, w))
		continue;
	    c0 = c1 = 0;
	    if (r > 0) {
		c0 += in_cross(v, w);
		c1 += in_cross(w, v);
	    }

	    if (GD_rank(g)[r + 1].n > 0) {
		c0 += out_cross(v, w);
		c1 += out_cross(w, v);
	    }

	    if (c1 <= c0) {
		balanceNodes(g, r, v, w);
	    }
	}
    }
    return rv;
}

static char *get_name(void *obj)
{
	// 获取节点或者边的自定义名字，按照预定义规则
	Agnode_t *n, *head, *tail;
	Agedge_t *e;

	if (AGTYPE(obj) == AGNODE) {
		// transform to type node n
		n = (node_t *)obj;
		if (ND_node_type(n) == NORMAL)
			return agnameof(n);
		
		static char buf[256];

		// 查找节点的输入边，追溯到非虚拟节点
		head = agtail(ND_in(n).list[0]);
		while (ND_node_type(head) == VIRTUAL) {
			assert(ND_in(head).size == 1);
			e = ND_in(head).list[0];
			head = agtail(e);
		}

		// 查找节点的输出边，追溯到非虚拟节点
		tail = aghead(ND_out(n).list[0]);
		while (ND_node_type(tail) == VIRTUAL) {
			assert(ND_out(tail).size == 1);
			e = ND_out(tail).list[0];
			tail = aghead(e);
		}

		sprintf(buf, "%d_%s_%s", ND_rank(n), agnameof(head), agnameof(tail));
		return buf;
	} else {
		e = (edge_t *)obj;
		head = agtail(e);
		tail = aghead(e);

		static char buf[256];
		sprintf(buf, "%s-%s", get_name(head), get_name(tail));
		return buf;
	}
}

static int get_type(void *obj)
{
	// 获取节点或者边的类型
	Agnode_t *n, *head, *tail;
	Agedge_t *e;

	if (AGTYPE(obj) == AGNODE) {
		// transform to type node n
		n = (node_t *)obj;
		if (ND_node_type(n) == NORMAL) {
			if (agnameof(n)[0] == 'l' || agnameof(n)[0] == 'r')
				return 2;
			else
				return 0;
		}

		// 查找节点的输入边，追溯到非虚拟节点
		head = agtail(ND_in(n).list[0]);
		while (ND_node_type(head) == VIRTUAL) {
			assert(ND_in(head).size == 1);
			e = ND_in(head).list[0];
			head = agtail(e);
		}

		// 查找节点的输出边，追溯到非虚拟节点
		tail = aghead(ND_out(n).list[0]);
		while (ND_node_type(tail) == VIRTUAL) {
			assert(ND_out(tail).size == 1);
			e = ND_out(tail).list[0];
			tail = aghead(e);
		}

		if (agnameof(head)[0] == 'l' || agnameof(tail)[0] == 'r')
			return 3;
		else
			return 1;
	} else {
		e = (edge_t *)obj;
		head = agtail(e);
		tail = aghead(e);
		int head_type = get_type(head);
		int tail_type = get_type(tail);

		if (head_type == 0 && tail_type == 0) return 0;
		if (head_type == 1 && tail_type == 0 ||\
		 	head_type == 0 && tail_type == 1) return 1;
		if (head_type == 1 && tail_type == 1) return 2;
		if (head_type == 0 && tail_type == 2 ||\
			head_type == 2 && tail_type == 0) return 3;
		if (head_type == 0 && tail_type == 3 ||\
			head_type == 3 && tail_type == 0) return 4;
		if (head_type == 3 && tail_type == 2 ||\
			head_type == 2 && tail_type == 3) return 5;
		if (head_type == 3 && tail_type == 3) return 6;
		if (head_type == 2 && tail_type == 2) return 7;
		
		return 8;
	}
}

static void examine_order(graph_t * g) {
	// for each rank, examine the node start with l and end with r
	// otherwise, exchange the node
	int r, i;
	node_t *v, *w;

	fprintf(stderr, "[sy]  - start examine_order\n");

	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
		fprintf(stderr, "[sy]    - examine_order(before) of rank %d, first %s, last %s\n", r, agnameof(GD_rank(g)[r].v[0]), agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]));
		if (agnameof(GD_rank(g)[r].v[0])[0] != 'l') {
			for (i = 1; i < GD_rank(g)[r].n; i++) {
				if (agnameof(GD_rank(g)[r].v[i])[0] == 'l') {
					// 不是和第一个节点交换，从0到i-1往后移动，再把i放到0
					// exchange(GD_rank(g)[r].v[0], GD_rank(g)[r].v[i]);
					node_t *tmp = GD_rank(g)[r].v[i];
					for (int j = i; j > 0; j--) {
						GD_rank(g)[r].v[j] = GD_rank(g)[r].v[j - 1];
						ND_order(GD_rank(g)[r].v[j]) = j;
						fprintf(stderr, "[sy]    - Move node %s from %d to %d\n", agnameof(GD_rank(g)[r].v[j]), j-1, j);
					}
					GD_rank(g)[r].v[0] = tmp;
					ND_order(GD_rank(g)[r].v[0]) = 0;
					fprintf(stderr, "[sy]    - Move node %s from %d to %d\n", agnameof(GD_rank(g)[r].v[0]), i, 0);
					break;
				}
			}
		}
		if (agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1])[0] != 'r') {
			for (i = GD_rank(g)[r].n - 2; i >= 0; i--) {
				if (agnameof(GD_rank(g)[r].v[i])[0] == 'r') {
					// 不是和最后一个节点交换，从i+1到n-1往前移动，再把i放到n-1
					// exchange(GD_rank(g)[r].v[GD_rank(g)[r].n - 1], GD_rank(g)[r].v[i]);
					node_t *tmp = GD_rank(g)[r].v[i];
					for (int j = i; j < GD_rank(g)[r].n - 1; j++) {
						GD_rank(g)[r].v[j] = GD_rank(g)[r].v[j + 1];
						ND_order(GD_rank(g)[r].v[j]) = j;
						fprintf(stderr, "[sy]    - Move node %s from %d to %d\n", agnameof(GD_rank(g)[r].v[j]), j+1, j);
					}
					GD_rank(g)[r].v[GD_rank(g)[r].n - 1] = tmp;
					ND_order(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]) = GD_rank(g)[r].n - 1;
					fprintf(stderr, "[sy]    - Move node %s from %d to %d\n", agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]), i, GD_rank(g)[r].n - 1);
					break;
				}
			}
		}
		fprintf(stderr, "[sy]    - examine_order(after) of rank %d, first %s, last %s\n", r, agnameof(GD_rank(g)[r].v[0]), agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]));
	}
}

static int transpose_step(graph_t * g, int r, bool reverse)
{
    int i, c0, c1, rv;
    node_t *v, *w;

    rv = 0;
    GD_rank(g)[r].candidate = false;

	// fprintf(stderr, "[sy]  - transpose_step(before) of rank %d, first %s, last %s\n", r, agnameof(GD_rank(g)[r].v[0]), agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]));

    for (i = 0; i < GD_rank(g)[r].n - 1; i++) {
	v = GD_rank(g)[r].v[i];
	w = GD_rank(g)[r].v[i + 1];
	assert(ND_order(v) < ND_order(w));
	if (left2right(g, v, w))
	    continue;
	if (is_context(v) || is_context(w))
	    continue;
	c0 = c1 = 0;
	if (r > 0) {
	    c0 += in_cross(v, w);
	    c1 += in_cross(w, v);
	}
	if (GD_rank(g)[r + 1].n > 0) {
	    c0 += out_cross(v, w);
	    c1 += out_cross(w, v);
	}
	if (c1 < c0 || (c0 > 0 && reverse && c1 == c0)) {
		exchange(v, w);
		rv += c0 - c1;
		GD_rank(Root)[r].valid = false;
		GD_rank(g)[r].candidate = true;

		if (r > GD_minrank(g)) {
		GD_rank(Root)[r - 1].valid = false;
		GD_rank(g)[r - 1].candidate = true;
		}
		if (r < GD_maxrank(g)) {
		GD_rank(Root)[r + 1].valid = false;
		GD_rank(g)[r + 1].candidate = true;
		}
	}
    }

	// fprintf(stderr, "[sy]  - transpose_step(after) of rank %d, first %s, last %s\n", r, agnameof(GD_rank(g)[r].v[0]), agnameof(GD_rank(g)[r].v[GD_rank(g)[r].n - 1]));
    return rv;
}

static void transpose(graph_t * g, bool reverse)
{
    int r, delta;

	fprintf(stderr, "[sy]- transpose of graph %s\n", agnameof(g));

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++)
	GD_rank(g)[r].candidate = true;
    do {
	delta = 0;
	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	    if (GD_rank(g)[r].candidate) {
		delta += transpose_step(g, r, reverse);
	    }
	}
    } while (delta >= 1);
}


// static void rename_virtual_nodes(graph_t * g) {
// 	int r, i;
// 	node_t *v;
// 	char buf[SMALLBUF];

// 	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
// 	for (i = 0; i < GD_rank(g)[r].n; i++) {
// 	    v = GD_rank(g)[r].v[i];
// 	    if (ND_node_type(v) == VIRTUAL) {
// 		sprintf(buf, "%s%d", ND_ranktype(v) == CLUSTER ? "c" : "l", i);
// 		agstrdup(&ND_node_name(v), buf);
// 	    }
// 	}
// 	}
// }

static void rename_virtual_nodes(Agraph_t *g) {
    Agnode_t *n, *head, *tail;
    Agedge_t *e;
    char new_name[1024];

	fprintf(stderr, "[sy]- rename_virtual_nodes\n");
    // 遍历所有节点
    // for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
	for (int r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	for (int i = 0; i < GD_rank(g)[r].n; i++) {
		n = GD_rank(g)[r].v[i];
		// 检查节点是否为虚拟节点
		if (ND_node_type(n) == VIRTUAL) {
			fprintf(stderr, "[sy]   - rename_virtual_nodes, node %s, type %d\n", agnameof(n), ND_node_type(n));
			// 断言虚拟节点有且只有一个入边和一个出边
			assert(ND_in(n).size == 1 && ND_out(n).size == 1);

			// 查找节点的输入边，追溯到非虚拟节点
			head = agtail(ND_in(n).list[0]);
			while (ND_node_type(head) == VIRTUAL) {
				fprintf(stderr, "[sy]   - rename_virtual_nodes, head %s is virtual\n", agnameof(head));
				assert(ND_in(head).size == 1);
				e = ND_in(head).list[0];
				head = agtail(e);
			}

			// 查找节点的输出边，追溯到非虚拟节点
			tail = aghead(ND_out(n).list[0]);
			while (ND_node_type(tail) == VIRTUAL) {
				fprintf(stderr, "[sy]   - rename_virtual_nodes, tail %s is virtual\n", agnameof(tail));
				assert(ND_out(tail).size == 1);
				e = ND_out(tail).list[0];
				tail = aghead(e);
			}

			// 生成新的虚拟节点名称
			sprintf(new_name, "%d_%s_%s", ND_rank(n), agnameof(tail), agnameof(head));
			fprintf(stderr, "[sy] - Setting attribute %s to %s on %s\n", "label", new_name, agnameof(n));
			// CGHDR_API Agattr_t *agattrrec(void *obj);
			// Segmentation fault when setting virtual node!!!!!!!!!!
			// virtual node has no data, so we cannot save attribute to it!!!
			agsafeset(n, "label", new_name, "");
			fprintf(stderr, "[sy]   - rename_virtual_nodes, set label %s to %s\n", agnameof(n), new_name);
		}
	}
    }
}

int get_crossing_type(graph_t *g) {
	const char *crossing_type = agget(g, "crossing_type");
	// convert string to int
	if (crossing_type) {
		if (strcmp(crossing_type, "0") == 0) return 0;
		if (strcmp(crossing_type, "1") == 0) return 1;
		if (strcmp(crossing_type, "2") == 0) return 2;
		if (strcmp(crossing_type, "3") == 0) return 3;
	}
	// raise error
	fprintf(stderr, "[sy] - get_crossing_type, crossing_type is not set or invalid\n");
}

static void print_ranks(graph_t * g) {
	int r, i;
	node_t *v;

	fprintf(stderr, "[sy]- print_ranks\n");
	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
		fprintf(stderr, "[sy]  - print_ranks on rank %d: ", r);
		for (i = 0; i < GD_rank(g)[r].n; i++) {
			v = GD_rank(g)[r].v[i];
			fprintf(stderr, "%s, ", get_name(v));
		}
		fprintf(stderr, "\n");
	}

	// const char *context_weight = agget(g, "context_weight");
	// if (context_weight) {
	// 	fprintf(stderr, "[sy]  - print_ranks, read context_weight: %s\n", context_weight);
	// } else {
	// 	agattr(g, AGRAPH, "context_weight", "0.5");
	// 	fprintf(stderr, "[sy]  - print_ranks, set context_weight: %s\n", "0.5");
	// }
}


static int mincross(graph_t * g, int startpass, int endpass, int doBalance)
{
	fprintf(stderr, "[sy]- mincross of graph %s (%d, %d, %d)\n", agnameof(g), startpass, endpass, doBalance);
    int maxthispass = 0, iter, trying, pass;
    int cur_cross, best_cross;

    if (startpass > 1) {
	cur_cross = best_cross = ncross(g);
	save_best(g);
    } else
	cur_cross = best_cross = INT_MAX;
    for (pass = startpass; pass <= endpass; pass++) {
		if (pass <= 1) {
			maxthispass = MIN(4, MaxIter);
			if (g == dot_root(g))
			build_ranks(g, pass);
			if (pass == 0)
			flat_breakcycles(g);
			flat_reorder(g);
			// 在flat_reorder中，打破了节点顺序约束，所以需要重新检查
			fprintf(stderr, "[sy]  - mincross, after flat_reorder: examine_order\n");
			examine_order(g);

			print_ranks(g);

			if ((cur_cross = ncross(g)) <= best_cross) {
			save_best(g);
			best_cross = cur_cross;
			}
		} else {
			maxthispass = MaxIter;
			if (cur_cross > best_cross)
			restore_best(g);
			cur_cross = best_cross;
		}
		trying = 0;
		for (iter = 0; iter < maxthispass; iter++) {
			if (Verbose)
			fprintf(stderr,
				"mincross: pass %d iter %d trying %d cur_cross %d best_cross %d\n",
				pass, iter, trying, cur_cross, best_cross);
			if (trying++ >= MinQuit)
			break;
			if (cur_cross == 0)
			break;
			mincross_step(g, iter);	
			
			if ((cur_cross = ncross(g)) <= best_cross) {
			save_best(g);
			if (cur_cross < Convergence * best_cross)
				trying = 0;
			best_cross = cur_cross;
			}
		}
		if (cur_cross == 0)
			break;
    }

    if (cur_cross > best_cross)
	restore_best(g);
    if (best_cross > 0) {
	transpose(g, FALSE);
	best_cross = ncross(g);
    }
	fprintf(stderr, "[sy]- mincross of graph %s done, best_cross: %d\n", agnameof(g), best_cross);
    if (doBalance) {
		fprintf(stderr, "[sy]- mincross of graph %s, doBalance for %d times\n", agnameof(g), maxthispass);
	for (iter = 0; iter < maxthispass; iter++)
	    balance(g);
    }

    return best_cross;
}

static void restore_best(graph_t * g)
{
    node_t *n;
    int i, r;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    n = GD_rank(g)[r].v[i];
	    ND_order(n) = saveorder(n);
	}
    }
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	GD_rank(Root)[r].valid = false;
	qsort(GD_rank(g)[r].v, GD_rank(g)[r].n, sizeof(GD_rank(g)[0].v[0]),
	      (qsort_cmpf) nodeposcmpf);
    }
}

static void save_best(graph_t * g)
{
    node_t *n;
    int i, r;
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    n = GD_rank(g)[r].v[i];
	    saveorder(n) = ND_order(n);
	}
    }
}

/* merges the connected components of g */
static void merge_components(graph_t * g)
{
    node_t *u, *v;

    if (GD_comp(g).size <= 1)
	return;
    u = NULL;
    for (size_t c = 0; c < GD_comp(g).size; c++) {
	v = GD_comp(g).list[c];
	if (u)
	    ND_next(u) = v;
	ND_prev(v) = u;
	while (ND_next(v)) {
	    v = ND_next(v);
	}
	u = v;
    }
    GD_comp(g).size = 1;
    GD_nlist(g) = GD_comp(g).list[0];
    GD_minrank(g) = GlobalMinRank;
    GD_maxrank(g) = GlobalMaxRank;
}

/* merge connected components, create globally consistent rank lists */
static void merge2(graph_t * g)
{
    int i, r;
    node_t *v;

    /* merge the components and rank limits */
    merge_components(g);

    /* install complete ranks */
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	GD_rank(g)[r].n = GD_rank(g)[r].an;
	GD_rank(g)[r].v = GD_rank(g)[r].av;
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    v = GD_rank(g)[r].v[i];
	    if (v == NULL) {
		if (Verbose)
		    fprintf(stderr,
			    "merge2: graph %s, rank %d has only %d < %d nodes\n",
			    agnameof(g), r, i, GD_rank(g)[r].n);
		GD_rank(g)[r].n = i;
		break;
	    }
	    ND_order(v) = i;
	}
    }
}

static void cleanup2(graph_t * g, int nc)
{
    int i, j, r, c;
    node_t *v;
    edge_t *e;

    if (TI_list) {
	free(TI_list);
	TI_list = NULL;
    }
    if (TE_list) {
	free(TE_list);
	TE_list = NULL;
    }
    /* fix vlists of clusters */
    for (c = 1; c <= GD_n_cluster(g); c++)
	rec_reset_vlists(GD_clust(g)[c]);

    /* remove node temporary edges for ordering nodes */
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    v = GD_rank(g)[r].v[i];
	    ND_order(v) = i;
	    if (ND_flat_out(v).list) {
		for (j = 0; (e = ND_flat_out(v).list[j]); j++)
		    if (ED_edge_type(e) == FLATORDER) {
			delete_flat_edge(e);
			free(e->base.data);
			free(e);
			j--;
		    }
	    }
	}
	free_matrix(GD_rank(g)[r].flat);
    }
    if (Verbose)
	fprintf(stderr, "mincross %s: %d crossings, %.2f secs.\n",
		agnameof(g), nc, elapsed_sec());
	
	fprintf(stderr, "[sy]- Metrics:\n");
	fprintf(stderr, "[sy]  - all cross: %d\n", ncross_with_type(g, 0));
	fprintf(stderr, "[sy]  - weigted cross: %d\n", ncross_with_type(g, 1));
	fprintf(stderr, "[sy]  - internal cross: %d\n", ncross_with_type(g, 2));
	fprintf(stderr, "[sy]  - external cross: %d\n", ncross_with_type(g, 3));
}

static node_t *neighbor(node_t * v, int dir)
{
    node_t *rv;

    rv = NULL;
assert(v);
    if (dir < 0) {
	if (ND_order(v) > 0)
	    rv = GD_rank(Root)[ND_rank(v)].v[ND_order(v) - 1];
    } else
	rv = GD_rank(Root)[ND_rank(v)].v[ND_order(v) + 1];
assert((rv == 0) || (ND_order(rv)-ND_order(v))*dir > 0);
    return rv;
}

static bool is_a_normal_node_of(graph_t *g, node_t *v) {
    return ND_node_type(v) == NORMAL && agcontains(g, v);
}

static bool is_a_vnode_of_an_edge_of(graph_t *g, node_t *v) {
    if (ND_node_type(v) == VIRTUAL
	&& ND_in(v).size == 1 && ND_out(v).size == 1) {
	edge_t *e = ND_out(v).list[0];
	while (ED_edge_type(e) != NORMAL)
	    e = ED_to_orig(e);
	if (agcontains(g, e))
	    return true;
    }
    return false;
}

static bool inside_cluster(graph_t *g, node_t *v) {
    return is_a_normal_node_of(g, v) | is_a_vnode_of_an_edge_of(g, v);
}

static node_t *furthestnode(graph_t * g, node_t * v, int dir)
{
    node_t *u, *rv;

    rv = u = v;
    while ((u = neighbor(u, dir))) {
	if (is_a_normal_node_of(g, u))
	    rv = u;
	else if (is_a_vnode_of_an_edge_of(g, u))
	    rv = u;
    }
    return rv;
}

void save_vlist(graph_t * g)
{
    int r;

    if (GD_rankleader(g))
	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	    GD_rankleader(g)[r] = GD_rank(g)[r].v[0];
	}
}

void rec_save_vlists(graph_t * g)
{
    int c;

    save_vlist(g);
    for (c = 1; c <= GD_n_cluster(g); c++)
	rec_save_vlists(GD_clust(g)[c]);
}


void rec_reset_vlists(graph_t * g)
{
    int r, c;
    node_t *u, *v, *w;

    /* fix vlists of sub-clusters */
    for (c = 1; c <= GD_n_cluster(g); c++)
	rec_reset_vlists(GD_clust(g)[c]);

    if (GD_rankleader(g))
	for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	    v = GD_rankleader(g)[r];
#ifdef DEBUG
	    node_in_root_vlist(v);
#endif
	    u = furthestnode(g, v, -1);
	    w = furthestnode(g, v, 1);
	    GD_rankleader(g)[r] = u;
#ifdef DEBUG
	    assert(GD_rank(dot_root(g))[r].v[ND_order(u)] == u);
#endif
	    GD_rank(g)[r].v = GD_rank(dot_root(g))[r].v + ND_order(u);
	    GD_rank(g)[r].n = ND_order(w) - ND_order(u) + 1;
	}
}

/* realFillRanks:
 * The structures in crossing minimization and positioning require
 * that clusters have some node on each rank. This function recursively
 * guarantees this property. It takes into account nodes and edges in
 * a cluster, the latter causing dummy nodes for intervening ranks.
 * For any rank without node, we create a real node of small size. This
 * is stored in the subgraph sg, for easy removal later.
 *
 * I believe it is not necessary to do this for the root graph, as these
 * are laid out one component at a time and these will necessarily have a
 * node on each rank from source to sink levels.
 */
static Agraph_t*
realFillRanks (Agraph_t* g, int rnks[], int rnks_sz, Agraph_t* sg)
{
    int i, c;
    Agedge_t* e;
    Agnode_t* n;

    for (c = 1; c <= GD_n_cluster(g); c++)
	sg = realFillRanks (GD_clust(g)[c], rnks, rnks_sz, sg);

    if (dot_root(g) == g)
	return sg;
    memset (rnks, 0, sizeof(int)*rnks_sz);
    for (n = agfstnode(g); n; n = agnxtnode(g,n)) {
	rnks[ND_rank(n)] = 1;
	for (e = agfstout(g,n); e; e = agnxtout(g,e)) {
	    for (i = ND_rank(n)+1; i <= ND_rank(aghead(e)); i++) 
		rnks[i] = 1;
	}
    }
    for (i = GD_minrank(g); i <= GD_maxrank(g); i++) {
	if (rnks[i] == 0) {
	    if (!sg) {
		sg = agsubg (dot_root(g), "_new_rank", 1);
	    }
	    n = agnode (sg, NULL, 1);
	    agbindrec(n, "Agnodeinfo_t", sizeof(Agnodeinfo_t), true);
	    ND_rank(n) = i;
	    ND_lw(n) = ND_rw(n) = 0.5;
	    ND_ht(n) = 1;
	    ND_UF_size(n) = 1;
	    alloc_elist(4, ND_in(n));
	    alloc_elist(4, ND_out(n));
	    agsubnode (g, n, 1);
	}
    }
    return sg;
}

static void
fillRanks (Agraph_t* g)
{
    int rnks_sz = GD_maxrank(g) + 2;
    int *rnks = gv_calloc(rnks_sz, sizeof(int));
    realFillRanks (g, rnks, rnks_sz, NULL);
    free (rnks);
}

static void init_mincross(graph_t * g)
{
    int size;

    if (Verbose)
	start_timer();

    ReMincross = false;
    Root = g;
    /* alloc +1 for the null terminator usage in do_ordering() */
    size = agnedges(dot_root(g)) + 1;
    TE_list = gv_calloc(size, sizeof(edge_t*));
    TI_list = gv_calloc(size, sizeof(int));
    mincross_options(g);
    if (GD_flags(g) & NEW_RANK)
	fillRanks (g);
    class2(g);
    decompose(g, 1);
    allocate_ranks(g);
    ordered_edges(g);
    GlobalMinRank = GD_minrank(g);
    GlobalMaxRank = GD_maxrank(g);
}

static void flat_rev(Agraph_t * g, Agedge_t * e)
{
    int j;
    Agedge_t *rev;

    if (!ND_flat_out(aghead(e)).list)
	rev = NULL;
    else
	for (j = 0; (rev = ND_flat_out(aghead(e)).list[j]); j++)
	    if (aghead(rev) == agtail(e))
		break;
    if (rev) {
	merge_oneway(e, rev);
	if (ED_edge_type(rev) == FLATORDER && ED_to_orig(rev) == 0)
	    ED_to_orig(rev) = e;
	elist_append(e, ND_other(agtail(e)));
    } else {
	rev = new_virtual_edge(aghead(e), agtail(e), e);
	if (ED_edge_type(e) == FLATORDER)
	    ED_edge_type(rev) = FLATORDER;
	else
	    ED_edge_type(rev) = REVERSED;
	ED_label(rev) = ED_label(e);
	flat_edge(g, rev);
    }
}

static void flat_search(graph_t * g, node_t * v)
{
    int i;
    bool hascl;
    edge_t *e;
    adjmatrix_t *M = GD_rank(g)[ND_rank(v)].flat;

    ND_mark(v) = TRUE;
    ND_onstack(v) = true;
    hascl = GD_n_cluster(dot_root(g)) > 0;
    if (ND_flat_out(v).list)
	for (i = 0; (e = ND_flat_out(v).list[i]); i++) {
	    if (hascl && !(agcontains(g, agtail(e)) && agcontains(g, aghead(e))))
		continue;
	    if (ED_weight(e) == 0)
		continue;
	    if (ND_onstack(aghead(e))) {
		assert(flatindex(aghead(e)) < M->nrows);
		assert(flatindex(agtail(e)) < M->ncols);
		ELT(M, flatindex(aghead(e)), flatindex(agtail(e))) = 1;
		delete_flat_edge(e);
		i--;
		if (ED_edge_type(e) == FLATORDER)
		    continue;
		flat_rev(g, e);
	    } else {
		assert(flatindex(aghead(e)) < M->nrows);
		assert(flatindex(agtail(e)) < M->ncols);
		ELT(M, flatindex(agtail(e)), flatindex(aghead(e))) = 1;
		if (!ND_mark(aghead(e)))
		    flat_search(g, aghead(e));
	    }
	}
    ND_onstack(v) = false;
}

static void flat_breakcycles(graph_t * g)
{
    int i, r, flat;
    node_t *v;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	flat = 0;
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    v = GD_rank(g)[r].v[i];
	    ND_mark(v) = FALSE;
	    ND_onstack(v) = false;
	    ND_low(v) = i;
	    if (ND_flat_out(v).size > 0 && flat == 0) {
		GD_rank(g)[r].flat =
		    new_matrix((size_t)GD_rank(g)[r].n, (size_t)GD_rank(g)[r].n);
		flat = 1;
	    }
	}
	if (flat) {
	    for (i = 0; i < GD_rank(g)[r].n; i++) {
		v = GD_rank(g)[r].v[i];
		if (!ND_mark(v))
		    flat_search(g, v);
	    }
	}
    }
}

/* allocate_ranks:
 * Allocate rank structure, determining number of nodes per rank.
 * Note that no nodes are put into the structure yet.
 */
void allocate_ranks(graph_t * g)
{
    int r, low, high;
    node_t *n;
    edge_t *e;

    int *cn = gv_calloc(GD_maxrank(g) + 2, sizeof(int)); // must be 0 based, not GD_minrank
    for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
	cn[ND_rank(n)]++;
	for (e = agfstout(g, n); e; e = agnxtout(g, e)) {
	    low = ND_rank(agtail(e));
	    high = ND_rank(aghead(e));
	    if (low > high) {
		int t = low;
		low = high;
		high = t;
	    }
	    for (r = low + 1; r < high; r++)
		cn[r]++;
	}
    }
    GD_rank(g) = gv_calloc(GD_maxrank(g) + 2, sizeof(rank_t));
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	GD_rank(g)[r].an = GD_rank(g)[r].n = cn[r];
	GD_rank(g)[r].av = GD_rank(g)[r].v = gv_calloc(cn[r] + 1, sizeof(node_t*));
    }
    free(cn);
}

/* install a node at the current right end of its rank */
void install_in_rank(graph_t * g, node_t * n)
{
    int i, r;

    r = ND_rank(n);
    i = GD_rank(g)[r].n;
    if (GD_rank(g)[r].an <= 0) {
	agerr(AGERR, "install_in_rank, line %d: %s %s rank %d i = %d an = 0\n",
	      __LINE__, agnameof(g), agnameof(n), r, i);
	return;
    }

    GD_rank(g)[r].v[i] = n;
    ND_order(n) = i;
    GD_rank(g)[r].n++;
    assert(GD_rank(g)[r].n <= GD_rank(g)[r].an);
#ifdef DEBUG
    {
	node_t *v;

	for (v = GD_nlist(g); v; v = ND_next(v))
	    if (v == n)
		break;
	assert(v != NULL);
    }
#endif
    if (ND_order(n) > GD_rank(Root)[r].an) {
	agerr(AGERR, "install_in_rank, line %d: ND_order(%s) [%d] > GD_rank(Root)[%d].an [%d]\n",
	      __LINE__, agnameof(n), ND_order(n), r, GD_rank(Root)[r].an);
	return;
    }
    if (r < GD_minrank(g) || r > GD_maxrank(g)) {
	agerr(AGERR, "install_in_rank, line %d: rank %d not in rank range [%d,%d]\n",
	      __LINE__, r, GD_minrank(g), GD_maxrank(g));
	return;
    }
    if (GD_rank(g)[r].v + ND_order(n) >
	GD_rank(g)[r].av + GD_rank(Root)[r].an) {
	agerr(AGERR, "install_in_rank, line %d: GD_rank(g)[%d].v + ND_order(%s) [%d] > GD_rank(g)[%d].av + GD_rank(Root)[%d].an [%d]\n",
	      __LINE__, r, agnameof(n),ND_order(n), r, r, GD_rank(Root)[r].an);
	return;
    }
}

/*	install nodes in ranks. the initial ordering ensure that series-parallel
 *	graphs such as trees are drawn with no crossings.  it tries searching
 *	in- and out-edges and takes the better of the two initial orderings.
 */
void build_ranks(graph_t * g, int pass)
{
	fprintf(stderr, "[sy]- Builindg ranks for graph %s\n", agnameof(g));
    int i, j;
    node_t *n, *n0;
    edge_t **otheredges;
    nodequeue *q;
	node_t **left_nodes = malloc(GD_n_nodes(g) * sizeof(node_t *));
    node_t **right_nodes = malloc(GD_n_nodes(g) * sizeof(node_t *));
	int left_count = 0, right_count = 0;

    q = new_queue(GD_n_nodes(g));
    for (n = GD_nlist(g); n; n = ND_next(n))
	MARK(n) = FALSE;

#ifdef DEBUG
    {
	edge_t *e;
	for (n = GD_nlist(g); n; n = ND_next(n)) {
	    for (i = 0; (e = ND_out(n).list[i]); i++)
		assert(!MARK(aghead(e)));
	    for (i = 0; (e = ND_in(n).list[i]); i++)
		assert(!MARK(agtail(e)));
	}
    }
#endif

    for (i = GD_minrank(g); i <= GD_maxrank(g); i++)
	GD_rank(g)[i].n = 0;

    for (n = GD_nlist(g); n; n = ND_next(n)) {
	otheredges = pass == 0 ? ND_in(n).list : ND_out(n).list;
	if (otheredges[0] != NULL)
	    continue;
	if (!MARK(n)) {
	    MARK(n) = TRUE;
	    enqueue(q, n);
	    while ((n0 = dequeue(q))) {

		fprintf(stderr, "[sy]- Builindg ranks: installing node %s in rank %d\n", get_name(n0), ND_rank(n0));
		if (ND_ranktype(n0) != CLUSTER) {
		    if (is_left_subgraph(n0)) {
				left_nodes[left_count++] = n0;  // 存储到左节点列表
			} else if (is_right_subgraph(n0)) {
				right_nodes[right_count++] = n0;  // 存储到右节点列表
			} else {
				install_in_rank(g, n0);
			}
		    enqueue_neighbors(q, n0, pass);
		} else {
		    install_cluster(g, n0, pass, q);
		}
	    }
	}
    }
	fprintf(stderr, "[sy]- Builindg ranks: left %d, right %d, maxrank %d\n", left_count, right_count, GD_maxrank(g));
	if (left_count == GD_maxrank(g) + 1 && right_count == GD_maxrank(g) + 1) {
		fprintf(stderr, "[sy]- Builindg ranks of context for graph %s: left=right=maxrank+1=%d\n", agnameof(g), left_count);
		for (int r = GD_minrank(g); r <= GD_maxrank(g); r++) {
			fprintf(stderr, "[sy]  - Builindg ranks of context on rank %d\n", r);
			node_t **vlist = GD_rank(g)[r].v;
			// 遍历left_nodes，寻找第一个ND_rank = r 的节点
			node_t *left_node = NULL, *right_node = NULL;
			for (int i = 0; i < left_count; i++) {
				if (ND_rank(left_nodes[i]) == r) {
					left_node = left_nodes[i];
					break;
				}
			}
			// 遍历right_nodes，寻找第一个ND_rank = r 的节点
			for (int i = 0; i < right_count; i++) {
				if (ND_rank(right_nodes[i]) == r) {
					right_node = right_nodes[i];
					break;
				}
			}
			int n = GD_rank(g)[r].n;
			for (int i = n - 1; i >= 0; i--) { // 为 left 节点腾出空间
				vlist[i + 1] = vlist[i];
				ND_order(vlist[i]) = i + 1;
				fprintf(stderr, "[sy]    - Move node %s from %d to %d\n", agnameof(vlist[i]), i, i + 1);
			}
			vlist[0] = left_node;  // 安装 left 节点
			ND_order(left_node) = 0;
			fprintf(stderr, "[sy]    - Install left node %s at %d\n", agnameof(left_node), 0);
			vlist[n + 1] = right_node;  // 安装 right 节点
			ND_order(right_node) = n + 1;
			fprintf(stderr, "[sy]    - Install right node %s at %d\n", agnameof(right_node), n + 1);
			GD_rank(g)[r].n += 2;  // 更新节点数量
		}

		free(left_nodes);
		free(right_nodes);
	}
	

    if (dequeue(q))
	agerr(AGERR, "surprise\n");
    for (i = GD_minrank(g); i <= GD_maxrank(g); i++) {
	GD_rank(Root)[i].valid = false;
	if (GD_flip(g) && GD_rank(g)[i].n > 0) {
	    node_t **vlist = GD_rank(g)[i].v;
	    int num_nodes_1 = GD_rank(g)[i].n - 1;
	    int half_num_nodes_1 = num_nodes_1 / 2;
	    for (j = 0; j <= half_num_nodes_1; j++)
		exchange(vlist[j], vlist[num_nodes_1 - j]);
	}
    }

    if (g == dot_root(g) && ncross(g) > 0)
	transpose(g, FALSE);
    free_queue(q);

	// rename_virtual_nodes(g);
}

void enqueue_neighbors(nodequeue * q, node_t * n0, int pass)
{
    edge_t *e;

    if (pass == 0) {
	for (size_t i = 0; i < ND_out(n0).size; i++) {
	    e = ND_out(n0).list[i];
	    if (!MARK(aghead(e))) {
		MARK(aghead(e)) = TRUE;
		enqueue(q, aghead(e));
	    }
	}
    } else {
	for (size_t i = 0; i < ND_in(n0).size; i++) {
	    e = ND_in(n0).list[i];
	    if (!MARK(agtail(e))) {
		MARK(agtail(e)) = TRUE;
		enqueue(q, agtail(e));
	    }
	}
    }
}

static bool constraining_flat_edge(Agraph_t *g, Agedge_t *e) {
  if (ED_weight(e) == 0)
    return false;
  if (!inside_cluster(g, agtail(e)))
    return false;
  if (!inside_cluster(g, aghead(e)))
    return false;
  return true;
}

/* construct nodes reachable from 'here' in post-order.
* This is the same as doing a topological sort in reverse order.
*/
static int postorder(graph_t * g, node_t * v, node_t ** list, int r)
{
    edge_t *e;
    int i, cnt = 0;

    MARK(v) = TRUE;
    if (ND_flat_out(v).size > 0) {
	for (i = 0; (e = ND_flat_out(v).list[i]); i++) {
	    if (!constraining_flat_edge(g, e)) continue;
	    if (!MARK(aghead(e)))
		cnt += postorder(g, aghead(e), list + cnt, r);
	}
    }
    assert(ND_rank(v) == r);
    list[cnt++] = v;
    return cnt;
}

static void flat_reorder(graph_t * g)
{
    int i, r, pos, n_search, local_in_cnt, local_out_cnt, base_order;
    node_t *v, **left, **right, *t;
    node_t **temprank = NULL;
    edge_t *flat_e, *e;

    if (!GD_has_flat_edges(g))
	return;
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	if (GD_rank(g)[r].n == 0) continue;
	base_order = ND_order(GD_rank(g)[r].v[0]);
	for (i = 0; i < GD_rank(g)[r].n; i++)
	    MARK(GD_rank(g)[r].v[i]) = FALSE;
	temprank = ALLOC(i + 1, temprank, node_t *);
	pos = 0;

	/* construct reverse topological sort order in temprank */
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    if (GD_flip(g)) v = GD_rank(g)[r].v[i];
	    else v = GD_rank(g)[r].v[GD_rank(g)[r].n - i - 1];

	    local_in_cnt = local_out_cnt = 0;
	    for (size_t j = 0; j < ND_flat_in(v).size; j++) {
		flat_e = ND_flat_in(v).list[j];
		if (constraining_flat_edge(g, flat_e)) local_in_cnt++;
	    }
	    for (size_t j = 0; j < ND_flat_out(v).size; j++) {
		flat_e = ND_flat_out(v).list[j];
		if (constraining_flat_edge(g, flat_e)) local_out_cnt++;
	    }
	    if ((local_in_cnt == 0) && (local_out_cnt == 0))
		temprank[pos++] = v;
	    else {
		if (!MARK(v) && local_in_cnt == 0) {
		    left = temprank + pos;
		    n_search = postorder(g, v, left, r);
		    pos += n_search;
		}
	    }
	}

	if (pos) {
	    if (!GD_flip(g)) {
		left = temprank;
		right = temprank + pos - 1;
		while (left < right) {
		    t = *left;
		    *left = *right;
		    *right = t;
		    left++;
		    right--;
		}
	    }
	    for (i = 0; i < GD_rank(g)[r].n; i++) {
		v = GD_rank(g)[r].v[i] = temprank[i];
		ND_order(v) = i + base_order;
	    }

	    /* nonconstraint flat edges must be made LR */
	    for (i = 0; i < GD_rank(g)[r].n; i++) {
		v = GD_rank(g)[r].v[i];
		if (ND_flat_out(v).list) {
		    for (size_t j = 0; (e = ND_flat_out(v).list[j]); j++) {
			if ( (!GD_flip(g) && ND_order(aghead(e)) < ND_order(agtail(e))) ||
				 ( (GD_flip(g)) && (ND_order(aghead(e)) > ND_order(agtail(e)) ))) {
			    assert(!constraining_flat_edge(g, e));
			    delete_flat_edge(e);
			    j--;
			    flat_rev(g, e);
			}
		    }
		}
	    }
	    /* postprocess to restore intended order */
	}
	/* else do no harm! */
	GD_rank(Root)[r].valid = false;
    }
    free(temprank);
}

static void reorder(graph_t * g, int r, bool reverse, bool hasfixed)
{
    int changed = 0, nelt;
    node_t **vlist = GD_rank(g)[r].v;
    node_t **lp, **rp, **ep = vlist + GD_rank(g)[r].n;

    for (nelt = GD_rank(g)[r].n - 1; nelt >= 0; nelt--) {
	lp = vlist;
	while (lp < ep) {
	    /* find leftmost node that can be compared */
	    while (lp < ep && ND_mval(*lp) < 0)
		lp++;
	    if (lp >= ep)
		break;
	    /* find the node that can be compared */
	    bool sawclust = false;
	    bool muststay = false;
	    for (rp = lp + 1; rp < ep; rp++) {
		if (sawclust && ND_clust(*rp))
		    continue;	/* ### */
		if (is_context(*lp) || is_context(*rp)) {
			muststay = true;
			break;
		}
		if (left2right(g, *lp, *rp)) {
		    muststay = true;
		    break;
		}
		if (ND_mval(*rp) >= 0)
		    break;
		if (ND_clust(*rp))
		    sawclust = true;	/* ### */
	    }
	    if (rp >= ep)
		break;
	    if (!muststay) {
		int p1 = ND_mval(*lp);
		int p2 = ND_mval(*rp);
		if (p1 > p2 || (p1 == p2 && reverse)) {
			exchange(*lp, *rp);
			changed++;
		}
	    }
	    lp = rp;
	}
	if (!hasfixed && !reverse)
	    ep--;
    }

    if (changed) {
	GD_rank(Root)[r].valid = false;
	if (r > 0)
	    GD_rank(Root)[r - 1].valid = false;
    }
}
static bool is_left_subgraph(node_t *n) {
    // graph_t* nodeGraph = agraphof(n);
	// const char* graphName = agnameof(nodeGraph);
	// if (Verbose)
	// fprintf(stderr,
	// 		"mincross: node %s is in graph %s\n",
	// 		agnameof(n), graphName);
	// if (strcmp(graphName, "left") == 0) {
	// 	return true;
	// }
    // return false;

	// char *label = agget(n, "label");  // 获取标签属性
	// if (Verbose)
	// fprintf(stderr, "[sy]	- Node: %s, Label: %s\n", agnameof(n), label ? label : "NULL");
    // if (label != NULL && strncmp(label, "left_", 5) == 0) {
    //     return true;
    // }
    // return false;

	const char *name = agnameof(n);
    return strncmp(name, "l", 1) == 0;
}

static bool is_right_subgraph(node_t *n) {
    const char *name = agnameof(n);
    return strncmp(name, "r", 1) == 0;
}

static bool is_context(node_t *n) {
    const char *name = agnameof(n);
    return strncmp(name, "l", 1) == 0 || strncmp(name, "r", 1) == 0;
}

static void mincross_step(graph_t * g, int pass)
{
    int r, other, first, last, dir;

    bool reverse = pass % 4 < 2;

    if (pass % 2 == 0) {	/* down pass */
	first = GD_minrank(g) + 1;
	if (GD_minrank(g) > GD_minrank(Root))
	    first--;
	last = GD_maxrank(g);
	dir = 1;
    } else {			/* up pass */
	first = GD_maxrank(g) - 1;
	last = GD_minrank(g);
	if (GD_maxrank(g) < GD_maxrank(Root))
	    first++;
	dir = -1;
    }

    for (r = first; r != last + dir; r += dir) {
	other = r - dir;
	bool hasfixed = medians(g, r, other);
	reorder(g, r, reverse, hasfixed);
    }
	fprintf(stderr, "[sy]    - mincross_step: examine_order\n");
	examine_order(g);
    transpose(g, !reverse);
}

static int local_cross(elist l, int dir, int type)
{
    int i, j;
    int cross = 0;
    edge_t *e, *f;
    bool is_out = dir > 0;
    for (i = 0; (e = l.list[i]); i++) {
	if (is_out)
	    for (j = i + 1; (f = l.list[j]); j++) {
		if ((ND_order(aghead(f)) - ND_order(aghead(e)))
			 * (ED_tail_port(f).p.x - ED_tail_port(e).p.x) < 0)
		    cross += get_edge_penalty(e) * get_edge_penalty(f);
	} else
	    for (j = i + 1; (f = l.list[j]); j++) {
		if ((ND_order(agtail(f)) - ND_order(agtail(e)))
			* (ED_head_port(f).p.x - ED_head_port(e).p.x) < 0)
		    cross += get_edge_penalty(e) * get_edge_penalty(f);
	    }
    }
    return cross;
}

static int get_edge_panelty_with_type(edge_t *e, int type) {
	if (type == 0) return ED_xpenalty(e);
	if (type == 1) return ED_weight(e);
	if (type == 2) return get_type(e) <= 2;
	if (type == 3) return get_type(e) > 2;
}

static int get_edge_penalty(edge_t *e) {
	graph_t *g = agraphof(agtail(e));
	return get_edge_panelty_with_type(e, get_crossing_type(g));
}

static int rcross(graph_t * g, int r, int type)
{
	// type=0: all crossing
	// type=1: weighted crossing
	// type=2: internal crossing
	// type=3: external crossing

    int top, bot, cross, max, i, k;
    node_t **rtop, *v;

    cross = 0;
    max = 0;
    rtop = GD_rank(g)[r].v;

    int *Count = gv_calloc(GD_rank(Root)[r + 1].n + 1, sizeof(int));

    for (top = 0; top < GD_rank(g)[r].n; top++) {
	edge_t *e;
	if (max > 0) {
	    for (i = 0; (e = ND_out(rtop[top]).list[i]); i++) {
		for (k = ND_order(aghead(e)) + 1; k <= max; k++)
		    cross += Count[k] * get_edge_panelty_with_type(e, type);
	    }
	}
	for (i = 0; (e = ND_out(rtop[top]).list[i]); i++) {
	    int inv = ND_order(aghead(e));
	    if (inv > max)
		max = inv;
	    Count[inv] += get_edge_panelty_with_type(e, type);

		// fprintf(stderr, "[sy]    - head %s, tail %s: %d, %d\n", get_name(agtail(e)), get_name(aghead(e)), get_type(agtail(e)), get_type(aghead(e)));
		// fprintf(stderr, "[sy]    - rcross: edge %s, type %d, weight: %d\n", get_name(e), get_type(e), ED_weight(e));
	}
    }
    for (top = 0; top < GD_rank(g)[r].n; top++) {
	v = GD_rank(g)[r].v[top];
	if (ND_has_port(v))
	    cross += local_cross(ND_out(v), 1, type);
    }
    for (bot = 0; bot < GD_rank(g)[r + 1].n; bot++) {
	v = GD_rank(g)[r + 1].v[bot];
	if (ND_has_port(v))
	    cross += local_cross(ND_in(v), -1, type);
    }
    free(Count);
    return cross;
}

int ncross(graph_t * g)
{
    int r, count, nc;

    g = Root;
    count = 0;
    for (r = GD_minrank(g); r < GD_maxrank(g); r++) {
	if (GD_rank(g)[r].valid)
	    count += GD_rank(g)[r].cache_nc;
	else {
	    nc = GD_rank(g)[r].cache_nc = rcross(g, r, get_crossing_type(g));
	    count += nc;
	    GD_rank(g)[r].valid = true;
	}
    }
    return count;
}

int ncross_with_type(graph_t * g, int type)
{
    int r, count, nc;

    g = Root;
    count = 0;
    for (r = GD_minrank(g); r < GD_maxrank(g); r++) {
	    nc = GD_rank(g)[r].cache_nc = rcross(g, r, type);
	    count += nc;
	    GD_rank(g)[r].valid = true;
    }
    return count;
}

static int ordercmpf(int *i0, int *i1)
{
  if (*i0 < *i1) {
    return -1;
  }
  if (*i0 > *i1) {
    return 1;
  }
  return 0;
}

/* flat_mval:
 * Calculate a mval for nodes with no in or out non-flat edges.
 * Assume (ND_out(n).size == 0) && (ND_in(n).size == 0)
 * Find flat edge a->n where a has the largest order and set
 * n.mval = a.mval+1, assuming a.mval is defined (>=0).
 * If there are no flat in edges, find flat edge n->a where a 
 * has the smallest order and set * n.mval = a.mval-1, assuming 
 * a.mval is > 0.
 * Return true if n.mval is left -1, indicating a fixed node for sorting.
 */
static bool flat_mval(node_t * n)
{
    int i;
    edge_t *e, **fl;
    node_t *nn;

    if (ND_flat_in(n).size > 0) {
	fl = ND_flat_in(n).list;
	nn = agtail(fl[0]);
	for (i = 1; (e = fl[i]); i++)
	    if (ND_order(agtail(e)) > ND_order(nn))
		nn = agtail(e);
	if (ND_mval(nn) >= 0) {
	    ND_mval(n) = ND_mval(nn) + 1;
	    return false;
	}
    } else if (ND_flat_out(n).size > 0) {
	fl = ND_flat_out(n).list;
	nn = aghead(fl[0]);
	for (i = 1; (e = fl[i]); i++)
	    if (ND_order(aghead(e)) < ND_order(nn))
		nn = aghead(e);
	if (ND_mval(nn) > 0) {
	    ND_mval(n) = ND_mval(nn) - 1;
	    return false;
	}
    }
    return true;
}

#define VAL(node,port) (MC_SCALE * ND_order(node) + (port).order)

static bool medians(graph_t * g, int r0, int r1)
{
    int i, j0, lspan, rspan, *list;
    node_t *n, **v;
    edge_t *e;
    bool hasfixed = false;

    list = TI_list;
    v = GD_rank(g)[r0].v;
    for (i = 0; i < GD_rank(g)[r0].n; i++) {
	n = v[i];
	size_t j = 0;
	if (r1 > r0)
	    for (j0 = 0; (e = ND_out(n).list[j0]); j0++) {
		if (get_edge_penalty(e) > 0)
		    list[j++] = VAL(aghead(e), ED_head_port(e));
	} else
	    for (j0 = 0; (e = ND_in(n).list[j0]); j0++) {
		if (get_edge_penalty(e) > 0)
		    list[j++] = VAL(agtail(e), ED_tail_port(e));
	    }
	switch (j) {
	case 0:
	    ND_mval(n) = -1;
	    break;
	case 1:
	    ND_mval(n) = list[0];
	    break;
	case 2:
	    ND_mval(n) = (list[0] + list[1]) / 2;
	    break;
	default:
	    qsort(list, j, sizeof(int), (qsort_cmpf) ordercmpf);
	    if (j % 2)
		ND_mval(n) = list[j / 2];
	    else {
		/* weighted median */
		size_t rm = j / 2;
		size_t lm = rm - 1;
		rspan = list[j - 1] - list[rm];
		lspan = list[lm] - list[0];
		if (lspan == rspan)
		    ND_mval(n) = (list[lm] + list[rm]) / 2;
		else {
		    double w = list[lm] * (double)rspan + list[rm] * (double)lspan;
		    ND_mval(n) = w / (lspan + rspan);
		}
	    }
	}
    }
    for (i = 0; i < GD_rank(g)[r0].n; i++) {
	n = v[i];
	if ((ND_out(n).size == 0) && (ND_in(n).size == 0))
	    hasfixed |= flat_mval(n);
    }
    return hasfixed;
}

static int nodeposcmpf(node_t ** n0, node_t ** n1)
{
  if (ND_order(*n0) < ND_order(*n1)) {
    return -1;
  }
  if (ND_order(*n0) > ND_order(*n1)) {
    return 1;
  }
  return 0;
}

static int edgeidcmpf(edge_t ** e0, edge_t ** e1)
{
  if (AGSEQ(*e0) < AGSEQ(*e1)) {
    return -1;
  }
  if (AGSEQ(*e0) > AGSEQ(*e1)) {
    return 1;
  }
  return 0;
}

/* following code deals with weights of edges of "virtual" nodes */
#define ORDINARY	0
#define SINGLETON	1
#define	VIRTUALNODE	2
#define NTYPES		3

#define C_EE		1
#define C_VS		2
#define C_SS		2
#define C_VV		4

static int table[NTYPES][NTYPES] = {
    /* ordinary */ {C_EE, C_EE, C_EE},
    /* singleton */ {C_EE, C_SS, C_VS},
    /* virtual */ {C_EE, C_VS, C_VV}
};

static int endpoint_class(node_t * n)
{
    if (ND_node_type(n) == VIRTUAL)
	return VIRTUALNODE;
    if (ND_weight_class(n) <= 1)
	return SINGLETON;
    return ORDINARY;
}

void virtual_weight(edge_t * e)
{
    int t;
    t = table[endpoint_class(agtail(e))][endpoint_class(aghead(e))];

    /* check whether the upcoming computation will overflow */
    assert(t >= 0);
    if (INT_MAX / t < ED_weight(e)) {
	agerr(AGERR, "overflow when calculating virtual weight of edge\n");
	graphviz_exit(EXIT_FAILURE);
    }

    ED_weight(e) *= t;
}

#ifdef DEBUG
void check_rs(graph_t * g, int null_ok)
{
    int i, r;
    node_t *v, *prev;

    fprintf(stderr, "\n\n%s:\n", agnameof(g));
    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	fprintf(stderr, "%d: ", r);
	prev = NULL;
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    v = GD_rank(g)[r].v[i];
	    if (v == NULL) {
		fprintf(stderr, "NULL\t");
		if (!null_ok)
		    abort();
	    } else {
		fprintf(stderr, "%s(%f)\t", agnameof(v), ND_mval(v));
		assert(ND_rank(v) == r);
		assert(v != prev);
		prev = v;
	    }
	}
	fprintf(stderr, "\n");
    }
}

void check_order(void)
{
    int i, r;
    node_t *v;
    graph_t *g = Root;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	assert(GD_rank(g)[r].v[GD_rank(g)[r].n] == NULL);
	for (i = 0; (v = GD_rank(g)[r].v[i]); i++) {
	    assert(ND_rank(v) == r);
	    assert(ND_order(v) == i);
	}
    }
}
#endif

static void mincross_options(graph_t * g)
{
    char *p;
    double f;

    /* set default values */
    MinQuit = 8;
    MaxIter = 24;
    Convergence = .995;

    p = agget(g, "mclimit");
    if (p && (f = atof(p)) > 0.0) {
	MinQuit = MAX(1, MinQuit * f);
	MaxIter = MAX(1, MaxIter * f);
    }
}

#ifdef DEBUG
void check_exchange(node_t * v, node_t * w)
{
    int i, r;
    node_t *u;

    if (ND_clust(v) == NULL && ND_clust(w) == NULL)
	return;
    assert(ND_clust(v) == NULL || ND_clust(w) == NULL);
    assert(ND_rank(v) == ND_rank(w));
    assert(ND_order(v) < ND_order(w));
    r = ND_rank(v);

    for (i = ND_order(v) + 1; i < ND_order(w); i++) {
	u = GD_rank(dot_root(v))[r].v[i];
	if (ND_clust(u))
	    abort();
    }
}

void check_vlists(graph_t * g)
{
    int c, i, j, r;
    node_t *u;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    u = GD_rank(g)[r].v[i];
	    j = ND_order(u);
	    assert(GD_rank(Root)[r].v[j] == u);
	}
	if (GD_rankleader(g)) {
	    u = GD_rankleader(g)[r];
	    j = ND_order(u);
	    assert(GD_rank(Root)[r].v[j] == u);
	}
    }
    for (c = 1; c <= GD_n_cluster(g); c++)
	check_vlists(GD_clust(g)[c]);
}

void node_in_root_vlist(node_t * n)
{
    node_t **vptr;

    for (vptr = GD_rank(Root)[ND_rank(n)].v; *vptr; vptr++)
	if (*vptr == n)
	    break;
    if (*vptr == 0)
	abort();
}
#endif				/* DEBUG code */
