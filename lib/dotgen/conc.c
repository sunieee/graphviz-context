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
 *	build edge_t concentrators for parallel edges with a common endpoint
 */

#include	<dotgen/dot.h>
#include	<stdbool.h>

#define		UP		0
#define		DOWN	1

static bool samedir(edge_t * e, edge_t * f)
{
    edge_t *e0, *f0;

    for (e0 = e; e0 != NULL && ED_edge_type(e0) != NORMAL; e0 = ED_to_orig(e0));
    if (e0 == NULL)
	return false;
    for (f0 = f; f0 != NULL && ED_edge_type(f0) != NORMAL; f0 = ED_to_orig(f0));
    if (f0 == NULL)
	return false;
    if (ED_conc_opp_flag(e0))
	return false;
    if (ED_conc_opp_flag(f0))
	return false;
    return ((ND_rank(agtail(f0)) - ND_rank(aghead(f0)))
	    * (ND_rank(agtail(e0)) - ND_rank(aghead(e0))) > 0);
}

static bool downcandidate(node_t * v)
{
    return ND_node_type(v) == VIRTUAL && ND_in(v).size == 1
	    && ND_out(v).size == 1 && ND_label(v) == NULL;
}

static bool bothdowncandidates(node_t * u, node_t * v)
{
    edge_t *e, *f;
    e = ND_in(u).list[0];
    f = ND_in(v).list[0];
    if (downcandidate(v) && agtail(e) == agtail(f)) {
	return samedir(e, f)
	    && portcmp(ED_tail_port(e), ED_tail_port(f)) == 0;
    }
    return false;
}

static bool upcandidate(node_t * v)
{
    return ND_node_type(v) == VIRTUAL && ND_out(v).size == 1
	    && ND_in(v).size == 1 && ND_label(v) == NULL;
}

static bool bothupcandidates(node_t * u, node_t * v)
{
    edge_t *e, *f;
    e = ND_out(u).list[0];
    f = ND_out(v).list[0];
    if (upcandidate(v) && aghead(e) == aghead(f)) {
	return samedir(e, f)
	    && portcmp(ED_head_port(e), ED_head_port(f)) == 0;
    }
    return false;
}

static void mergevirtual(graph_t * g, int r, int lpos, int rpos, int dir)
{
    int i, k;
    node_t *left, *right;
    edge_t *e, *f, *e0, *f0;

    left = GD_rank(g)[r].v[lpos];

	if (Verbose) fprintf(stderr, "[sy]  - mergevirtual: r=%d, lpos=%d, rpos=%d, dir=%d, left=%s\n", r, lpos, rpos, dir, get_name(left));
    /* merge all right nodes into the leftmost one */
    for (i = lpos + 1; i <= rpos; i++) {
	right = GD_rank(g)[r].v[i];
	if (Verbose) fprintf(stderr, "[sy]    - right=%s\n", get_name(right));
	if (dir == DOWN) {
	    while ((e = ND_out(right).list[0])) {
		for (k = 0; (f = ND_out(left).list[k]); k++)
		    if (aghead(f) == aghead(e))
			break;

		if (f == NULL)
		    f = virtual_edge(left, aghead(e), e);
		else 
			merge_with_penwidth(e, f, g);
		assert(ND_in(left).size == 1);
		f0 = ND_in(left).list[0];
		
		while ((e0 = ND_in(right).list[0])) {
			merge_with_penwidth(e0, f0, g);
		    /*ED_weight(f) += ED_weight(e0); */
		    delete_fast_edge(e0);
		}
		delete_fast_edge(e);
	    }
		// fprintf(stderr, "[sy]    - after:f0=%s, weight=%d\n", get_name(f0), ED_weight(f0));
	} else {
	    while ((e = ND_in(right).list[0])) {
		for (k = 0; (f = ND_in(left).list[k]); k++)
		    if (agtail(f) == agtail(e))
			break;
		if (f == NULL)
		    f = virtual_edge(agtail(e), left, e);
		else
			merge_with_penwidth(e, f, g);

		assert(ND_out(left).size == 1);
		f0 = ND_out(left).list[0];
		while ((e0 = ND_out(right).list[0])) {
		    merge_with_penwidth(e0, f0, g);
		    delete_fast_edge(e0);
		}
		delete_fast_edge(e);
	    }
	}
	assert(ND_in(right).size + ND_out(right).size == 0);
	delete_fast_node(g, right);
    }
    k = lpos + 1;
    i = rpos + 1;
    while (i < GD_rank(g)[r].n) {
	node_t *n;
	n = GD_rank(g)[r].v[k] = GD_rank(g)[r].v[i];
	ND_order(n) = k;
	k++;
	i++;
    }
    GD_rank(g)[r].n = k;
    GD_rank(g)[r].v[k] = NULL;
}

static void infuse(graph_t * g, node_t * n)
{
    node_t *lead;

    lead = GD_rankleader(g)[ND_rank(n)];
    if (lead == NULL || ND_order(lead) > ND_order(n))
	GD_rankleader(g)[ND_rank(n)] = n;
}

static int rebuild_vlists(graph_t * g)
{
    int c, i, r, maxi;
    node_t *n, *lead;
    edge_t *rep;

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++)
	GD_rankleader(g)[r] = NULL;
    dot_scan_ranks(g);
    for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
	infuse(g, n);
	for (edge_t *e = agfstout(g, n); e; e = agnxtout(g, e)) {
	    for (rep = e; ED_to_virt(rep); rep = ED_to_virt(rep));
	    while (rep != NULL && ND_rank(aghead(rep)) < ND_rank(aghead(e))) {
		infuse(g, aghead(rep));
		rep = ND_out(aghead(rep)).list[0];
	    }
	}
    }

    for (r = GD_minrank(g); r <= GD_maxrank(g); r++) {
	lead = GD_rankleader(g)[r];
	if (lead == NULL) {
		agerr(AGERR, "rebuild_vlists: lead is null for rank %d\n", r);
		return -1;
	}
	else if (GD_rank(dot_root(g))[r].v[ND_order(lead)] != lead) {
	    agerr(AGERR, "rebuild_vlists: rank lead %s not in order %d of rank %d\n", 
		agnameof(lead), ND_order(lead), r);
	    return -1;
	}
	GD_rank(g)[r].v =
	    GD_rank(dot_root(g))[r].v + ND_order((GD_rankleader(g)[r]));
	maxi = -1;
	for (i = 0; i < GD_rank(g)[r].n; i++) {
	    if ((n = GD_rank(g)[r].v[i]) == NULL)
		break;
	    if (ND_node_type(n) == NORMAL) {
		if (agcontains(g, n))
		    maxi = i;
		else
		    break;
	    } else {
		edge_t *e;
		for (e = ND_in(n).list[0]; e && ED_to_orig(e);
		     e = ED_to_orig(e));
		if (e && agcontains(g, agtail(e))
		    && agcontains(g, aghead(e)))
		    maxi = i;
	    }
	}
	if (maxi == -1)
	    agerr(AGWARN, "degenerate concentrated rank %s,%d\n", agnameof(g),
		  r);
	GD_rank(g)[r].n = maxi + 1;
    }

    for (c = 1; c <= GD_n_cluster(g); c++) {
	int ret = rebuild_vlists(GD_clust(g)[c]);
	if (ret != 0) {
	    return ret;
	}
    }
    return 0;
}

void dot_concentrate(graph_t * g)
{
    int c, r, leftpos, rightpos;
    node_t *left, *right;
	int concentrate_type = get_attr(g, "concentrate_type");
	fprintf(stderr, "[sy]- dot_concentrate: concentrate_type=%d\n", concentrate_type);

    if (GD_maxrank(g) - GD_minrank(g) <= 1)
	return;
    /* this is the downward looking pass. r is a candidate rank. */

	fprintf(stderr, "[sy]- dot_concentrate: start downward looking pass\n");
	if (Verbose) print_ranks(g);

    for (r = 1; GD_rank(g)[r + 1].n; r++) {
	for (leftpos = 0; leftpos < GD_rank(g)[r].n; leftpos++) {
	    left = GD_rank(g)[r].v[leftpos];
	    if (!downcandidate(left))
		continue;
		// 1: concentrate only context virtual node (type == 3)
		// 2: concentrate only focus virtual node (type == 1)
		if (concentrate_type == 1 && get_type(left) != 3)
		continue;
		if (concentrate_type == 2 && get_type(left) != 1)
		continue;
		int type = get_type(left);
	    for (rightpos = leftpos + 1; rightpos < GD_rank(g)[r].n;
		 rightpos++) {
		right = GD_rank(g)[r].v[rightpos];
		if (!bothdowncandidates(left, right))
		    break;
		// if concentrate >=1, type of right should be the same as left
		if (concentrate_type >= 1 && get_type(right) != type)
			break;
	    }
	    if (rightpos - leftpos > 1)
		mergevirtual(g, r, leftpos, rightpos - 1, DOWN);
	}
    }

	fprintf(stderr, "[sy]- dot_concentrate: start upward looking pass\n");
	if (Verbose) print_ranks(g);

    /* this is the corresponding upward pass */
    while (r > 0) {
	for (leftpos = 0; leftpos < GD_rank(g)[r].n; leftpos++) {
	    left = GD_rank(g)[r].v[leftpos];
	    if (!upcandidate(left))
		continue;
		if (concentrate_type == 1 && get_type(left) != 3)
		continue;
		if (concentrate_type == 2 && get_type(left) != 1)
		continue;
		int type = get_type(left);
	    for (rightpos = leftpos + 1; rightpos < GD_rank(g)[r].n;
		 rightpos++) {
		right = GD_rank(g)[r].v[rightpos];
		if (!bothupcandidates(left, right))
		    break;
	    }
		if (concentrate_type >= 1 && get_type(right) != type)
			break;
	    if (rightpos - leftpos > 1)
		mergevirtual(g, r, leftpos, rightpos - 1, UP);
	}
	r--;
    }

	fprintf(stderr, "dot_concentrate: start rebuild_vlists\n");
	if (Verbose) print_ranks(g);

    for (c = 1; c <= GD_n_cluster(g); c++) {
	if (rebuild_vlists(GD_clust(g)[c]) != 0) {
	    agerr(AGPREV, "concentrate=true may not work correctly.\n");
	    return;
	}
    }
	// set_penwidth(g);
}


void set_penwidth(graph_t *g) {
	/*
	在Graphviz中，直接控制一条边在不同分段上的粗细是比较复杂的，
	因为Graphviz本身并不直接支持在单条边的不同部分设置不同的penwidth。
	Graphviz中的边是被视为单一的对象，其属性（如penwidth）是均匀应用于整个边的。
	
	因此本函数不适用于直接控制一条边在不同分段上的粗细。需要D3.js对图形元素进行更精细的控制
	*/

	Agnode_t *n;
	Agedge_t *ep;
	int penwidth_type = get_attr(g, "penwidth_type");
	int focus_edge_weight = get_attr(g, "focus_edge_weight");

	if (penwidth_type == 0) {
		// 按照默认方式绘制，即penwitdh属性
		return;
	}

	fprintf(stderr, "[sy]- set_penwidth, penwidth_type=%d, focus_edge_weight=%d\n", penwidth_type, focus_edge_weight);
	// 遍历图中的所有节点
	for (n = agfstnode(g); n; n = agnxtnode(g, n)) {
		// 遍历与当前节点相关的所有边
		for (ep = agfstedge(g, n); ep; ep = agnxtedge(g, ep, n)) {
			if (penwidth_type == 4) {
				// 使用固定1的penwidth
				agsafeset(ep, "penwidth", "1", "");
				continue;
			}

			char *weight_str = agget(ep, "weight");  // 获取 weight 属性
			if (weight_str) {
				// 设置 penwidth 为 weight 的值
				int edge_type = get_type(ep);
				fprintf(stderr, "[sy]  - set_penwidth: ep=%s, type=%d, weight=%s\n", get_name(ep), edge_type, weight_str);
				if (edge_type <=2) {
					// focus edge: 0~2
					if (penwidth_type == 2 || penwidth_type == 3) {
						int weight = atoi(weight_str) / focus_edge_weight;
						char * weight_str = (char *)malloc(10);
						sprintf(weight_str, "%d", weight);
						agsafeset(ep, "penwidth", weight_str, "");
					}

				} else if (edge_type <= 6) {
					// context edge: 3~6
					if (penwidth_type == 1 || penwidth_type == 3) {
						agsafeset(ep, "penwidth", weight_str, "");
					}
				}
			}
		}
	}
}
