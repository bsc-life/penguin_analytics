import os
import re

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec
from matplotlib import colors as mcolors
import numpy as np
from random import random
import pandas as pd

import io
import base64

dpath = 'graph_tables'


def transform(x, y):
    """
    get coordinates further from circle center
    """
    d = (x**2 + y**2)**0.5
    a = np.arcsin(y / d)
    b = np.arccos(x / d)
    X = np.cos(b)
    Y = np.sin(a)
    k = 0.3 if d < 0.1 else 0.5 if d < 0.3 else 0.6 if d < 0.4 else 0.7 if d < 0.5 else 0.9
    k += random() * 0.1
    return X * k, Y * k

def annotate_axes(ax, text, fontsize=18):
    ax.text(0.5, 0.75, text, transform=ax.transAxes,
            ha="center", va="center", fontsize=fontsize, color="k", weight='normal')
    

def internal_norm(vals):
    mask = vals != np.inf
    vmax = np.nanmax(vals.loc[mask])
    vmin = np.nanmin(vals.loc[mask])
    vmax = max(vmax, abs(vmin))
    if vmax == True or vmax == False:
        vmax = 1
    if vmin < 0:  # some stats start at zero
        vmin = -vmax
    return np.array([v if np.isfinite(v) else vmax * 1.1 if v > 0 else vmin * 0.9
                    for v in vals]) / vmax


def single_promoter_graph(promoter, size_factor, top_genes, 
                          DBP_color, DBP_size,
                          intermediate_color, intermediate_size,
                          edge_color, edge_width, 
                          DBP_cmap, intermediate_cmap,
                          edge_cmap, edge_alpha,
                          Prom_displayed_perc, Enha_displayed_perc, int_displayed_perc, 
                          SNP_path, druggable, paintor_SNP, 
                         only_enriched_nodes,
                         enriched_edges,
                         paintor_SNP_edges):
    """
    :param 150 factor: scale factor (thehigher, the smaller willbe the balls)
    :param 40 top_genes: number of protein binding sites or intermediate protein names to be displayed
       in each ofthe 3 categories (promoter, middle or enhancer side)
    """

    edge_alpha = float(edge_alpha)
    node_druggability  = pd.read_csv("src/node_druggability.txt", sep='\t')
    paintor_snp_in_gene = pd.read_csv("src/nodes_with_paintorSNP.txt", sep='\t')
    paintor_snp_in_gene = paintor_snp_in_gene.groupby("protein").agg(
        n_rsid = pd.NamedAgg(column="rsid", aggfunc="nunique"))
    
    # read nodes
    nodes = pd.read_csv("".join([dpath,'/',re.sub('\W+','',str(promoter)),'_nodes.tsv']), sep='\t')
    nodes = nodes.merge(node_druggability, on='protein', how='left')
    nodes["Gene.name"] = nodes['protein'].str.replace("ENHA_", "").str.replace("PROM_", "")
    nodes = nodes.merge(paintor_snp_in_gene, left_on='Gene.name', right_on = "protein", how='left')
    nodes["paintor_SNP"] = np.where(nodes["n_rsid"] > 0  ,True , False)
    nodes["Node_druggability"] = nodes["Node_druggability"].fillna(False)
    #nodes = nodes.replace({'SNP-binding': {'True': True, 'False': False}})
    nodes['SNP-binding']=nodes['SNP-binding'].astype(object)
    nodes['SNP-binding'] = nodes['SNP-binding'].replace({'True': True, 'False': False})
    nodes["shape"] = np.where(
                nodes["Node_druggability"] == True, "^", "o" )

    paintor_SNP_nodes = list(set(nodes[(nodes["paintor_SNP"] == True) & (nodes["is_intermediate"] == True)]["protein"].values.tolist()))

    # read edges
    all_edges = pd.read_csv("".join([dpath,'/',re.sub('\W+','',str(promoter)),'_edges.tsv']), sep='\t')
    all_edges['XYs-edge'] = all_edges['XYs-edge'].apply(lambda x: eval(x.replace('array', '')))
    all_edges[['node1', 'node2']] = all_edges['edge'].str.split(';', expand=True)
      
    edges = pd.DataFrame()
    
    if only_enriched_nodes:
        #remove the non enriched nodes
        nodes =  nodes[(nodes["fisher_test_OR"] > 1) & (nodes["fisher_test_pval"] < 0.01)]
        #remove the edges attaching the non enriched nodes        
        all_edges = all_edges[(all_edges["node1"].isin(nodes["protein"])) & (all_edges["node2"].isin(nodes["protein"]))]
    
    if paintor_SNP_edges:
        edges = pd.concat([edges,
                           all_edges[(all_edges["node1"].isin(paintor_SNP_nodes)) | (all_edges["node2"].isin(paintor_SNP_nodes))] ])
    if SNP_path:
        edges = pd.concat([edges,
                           all_edges[all_edges['SNP-path'] == True]])
    if druggable:
        edges = pd.concat([edges,
                           all_edges[all_edges['druggable'] == True]])

    if druggable == False and SNP_path == False and enriched_edges == False and paintor_SNP_edges == False  :
        edges = all_edges

    # if i select the enricjed edges this will be a filter in the previous selection
    if enriched_edges:
        if edges.empty:
            edges = all_edges[(all_edges["fisher_test_edge_OR"] > 1) & (all_edges["fisher_test_edge_pvalue"] < 0.01)]
        else:
            edges = edges[(edges["fisher_test_edge_OR"] > 1) & (edges["fisher_test_edge_pvalue"] < 0.01)]
    
    edges = edges.loc[edges.astype(str).drop_duplicates().index]

    snp_path_nodes = list(set(all_edges[all_edges["SNP-path"] == True][["node1", "node2"]].stack().values.tolist()))
    if paintor_SNP == True:
        snp_path_nodes += nodes[nodes["paintor_SNP"] == True]["protein"].values.tolist()
        snp_path_nodes = list(set(snp_path_nodes))
        nodes["protein_display_name"] = np.where(
            nodes["paintor_SNP"] == True, nodes["protein"] + "*",  nodes["protein"]) 
    else:
        nodes["protein_display_name"] = nodes["protein"]


    ##################################
    # start plotting    
    fig = plt.figure(figsize=(14, 20), facecolor='white')

    gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[2, 1], hspace=0.1)
    gs00 = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=gs0[1], hspace=0.5)

    axe = fig.add_subplot(gs0[0])
    axes = [fig.add_subplot(gs00[i, j]) for i in range(4) for j in range(4)]

    #################
    # bottom histogram
    for nstat, stat in enumerate([
        'expression_log2FoldChange', 'expression_padj', 'fisher_test_OR', 'fisher_test_pval',
        'Betweenness_enrichment', 'Betweenness_enrichment_pval', 'Degree_enrichment', 'Degree_enrichment_pval',
        'Betweenness_global', 'Betweenness_in_cluster', 'Betweenness_out_cluster', 
        'Degree_global', 'Degree_in_cluster', 'Degree_out_cluster', ]):
        axes[nstat].hist([v for v in nodes[stat] if np.isfinite(v)], color='lightgrey', edgecolor='darkgrey')
        annotate_axes(axes[nstat], 'Node:\n' + stat.replace('_', ' '), fontsize=10)
        if stat in [DBP_color, DBP_size, intermediate_color, intermediate_size, edge_width, edge_color]:
            axes[nstat].spines['bottom'].set_color('tab:red')
            axes[nstat].spines['bottom'].set_linewidth(3)
            axes[nstat].spines['top'].set_color('tab:red') 
            axes[nstat].spines['top'].set_linewidth(3)
            axes[nstat].spines['right'].set_color('tab:red')
            axes[nstat].spines['right'].set_linewidth(3)
            axes[nstat].spines['left'].set_color('tab:red')
            axes[nstat].spines['left'].set_linewidth(3)
    for nstat, stat in enumerate(['fisher_test_edge_OR', 'fisher_test_edge_pvalue'], nstat + 1):
        axes[nstat].hist([v for v in edges[stat] if np.isfinite(v)], color='lightgrey', edgecolor='darkred')
        annotate_axes(axes[nstat], 'Edge:\n' + stat.replace('_', ' '), fontsize=10)
        if stat in [DBP_color, DBP_size, intermediate_color, intermediate_size, edge_width, edge_color]:
            axes[nstat].spines['bottom'].set_color('tab:red')
            axes[nstat].spines['bottom'].set_linewidth(3)
            axes[nstat].spines['top'].set_color('tab:red') 
            axes[nstat].spines['top'].set_linewidth(3)
            axes[nstat].spines['right'].set_color('tab:red')
            axes[nstat].spines['right'].set_linewidth(3)
            axes[nstat].spines['left'].set_color('tab:red')
            axes[nstat].spines['left'].set_linewidth(3)
    #################
    
    # X limits
    xprom = -1.1
    xenha = 1.1
    
    Prom_displayed_perc = 1 - Prom_displayed_perc / 100
    Enha_displayed_perc = 1 - Enha_displayed_perc / 100
    int_displayed_perc = 1 - int_displayed_perc / 100
        
    ## promoters
    sub_nodes = nodes[nodes['is_promoter']==True].sort_values(by='protein')
    
    # define size
    prom_sizes = internal_norm(sub_nodes[DBP_size])
    # plot promoter nodes
    for marker in ["o", "^"]:
        sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker]
        if len(sub_sub_nodes) == 0:
            continue
        # define size
        sub_sizes = internal_norm(sub_sub_nodes[DBP_size])
        # define colors
        prom_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_sub_nodes[DBP_color]))
        axe.scatter(sub_sub_nodes['X-graph-coord'], 
                    sub_sub_nodes['Y-graph-coord'],
                    s=(sub_sizes + 0.1) * size_factor,
                    marker = marker,
                    c=prom_colors)

    ## enhancers
    sub_nodes = nodes[nodes['is_enhancer']==True].sort_values(by='protein')
    # define size
    enha_sizes = internal_norm(sub_nodes[DBP_size])
    # plot enhancer nodes
    for marker in ["o", "^"]:
        sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker]
        if len(sub_sub_nodes) == 0:
            continue
        # define size    
        sub_sizes = internal_norm(sub_sub_nodes[DBP_size])
        # define colors
        enha_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_sub_nodes[DBP_color]))
        axe.scatter(sub_sub_nodes['X-graph-coord'], 
                    sub_sub_nodes['Y-graph-coord'], 
                    s=(sub_sizes + 0.1) * size_factor, 
                    marker = marker,
                    c=enha_colors)

        
        
    ## intermediates
    sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')
    

    # define size
    intr_sizes = internal_norm(sub_nodes[intermediate_size])
    sub_nodes["intr_sizes"] = intr_sizes
    # plot intermediate nodes
    
    if intermediate_color in ["fisher_test_OR" ] and only_enriched_nodes == False:
        sub_nodes["binary_sub_nodes"] = np.where((sub_nodes["fisher_test_OR"] > 1) & (sub_nodes["fisher_test_pval"] < 0.01 ), 1, 0)
 #       print(sub_nodes[sub_nodes["binary_sub_nodes"] == 1].shape)
 #       print(sub_nodes[sub_nodes["binary_sub_nodes"] == 0].shape)    
        intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_nodes["binary_sub_nodes"]))
 #       intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm_fisherOR(sub_nodes[intermediate_color]))

            
    else:
        intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_nodes[intermediate_color]))

    # define colors
    
    
    for marker in ["o", "^"]:
        #plot two markers, The normalisation has to be done before.
        # take the indeces
        sub_sub_nodes_index = (sub_nodes["shape"] == marker)
        
        #extract the colours
        sub_intr_colors = intr_colors[sub_sub_nodes_index]
        
        #extract the sizes
#         print(intr_sizes.size)        
#         print(intr_sizes)
#         print(sub_sub_nodes_index.size)
#         print(sub_sub_nodes_index)
        sub_sub_sizes = intr_sizes[sub_sub_nodes_index,]   
        
        #extract the nodes
        sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker] 
        
        if len(sub_sub_nodes) == 0:
            continue
        
        axe.scatter(sub_sub_nodes['X-graph-coord'],
                    sub_sub_nodes['Y-graph-coord'], 
                    s=(sub_sub_sizes + 0.1) * size_factor,
                    marker = marker,
                    c=sub_intr_colors)

    
    ### plot edges
    sub_edges = "empty"
    if edges.equals(all_edges) or edges.empty:
        filter_edges = False
    else:
        filter_edges = True
        sub_edges = [True] * len(edges.index)
#     if SNP_path and druggable:
#         sub_edges =  (edges['SNP-path']) & (edges['druggable'])
#         #sub_edges =  (edges['SNP-path']) | (edges['druggable'])
#     elif SNP_path:
#         sub_edges =  edges['SNP-path']
#     elif druggable:
#         sub_edges =  edges['druggable'] == druggable
#     else:
#         filter_edges = False
        
    
#     print("sub_edges", sub_edges)
#     print("edge_alpha", edge_alpha)
#     print("filter_edges", filter_edges)
#     print("edges", edges.shape)
    
    # fransuas initial line plotter
#     lc = LineCollection(edges['XYs-edge'], 
#                         alpha=(sub_edges * edge_alpha) if filter_edges else edge_alpha,
#                         #alpha= edge_alpha,
#                         lw=internal_norm(edges[edge_width]) * 0.8 + 0.5 if edge_width and filter_edges else 1, 
#                         color=(plt.cm.get_cmap(edge_cmap)(internal_norm(edges[edge_color]) * 0.9 + 0.1) if
#                                edge_cmap and edge_color else 'k'), zorder=-1)

    # alex  line plotter for two colours based on edge type

    conditions = [
        ((edges["node1"].isin(paintor_SNP_nodes)) | (edges["node2"].isin(paintor_SNP_nodes))) & (edges['SNP-path'] == True),
        (edges['SNP-path'] == True),
        ((edges["node1"].isin(paintor_SNP_nodes)) | (edges["node2"].isin(paintor_SNP_nodes)))
    ]
    #values = [mcolors.to_rgba("purple"), mcolors.to_rgba("red"), mcolors.to_rgba("blue")]
    choices = ['purple', 'orange', 'darkgreen']
    v = mcolors.to_rgba("black")
    edges["colour"] = np.select(conditions, choices, default = "black")
    lc = LineCollection(edges['XYs-edge'], 
                    alpha = 0.7,
                    lw = 0.7,
                    color = edges["colour"], 
                    zorder=-1)

    ### Annotation... this is the slowest part...
    ## promoter  
    sub_nodes = nodes[nodes['is_promoter']==True].sort_values(by='protein')
    sub_nodes["prom_sizes"] = prom_sizes 
    if SNP_path or len(snp_path_nodes) > 0:
        #sub_nodes = sub_nodes[(sub_nodes["prom_sizes"]  > Prom_displayed_perc ) | (sub_nodes["protein"].isin(snp_path_nodes))]
        sub_nodes = sub_nodes[(sub_nodes["protein"].isin(snp_path_nodes)) | 
                             (sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]
    else:
        sub_nodes = sub_nodes[(sub_nodes["prom_sizes"]  > Prom_displayed_perc )]
   
    #    sub_nodes = sub_nodes[sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist())))]
    ##print("prototer sub_nodes", sub_nodes.shape)
    wanted = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
    ypos = np.linspace(-1.15, 1.15, len(wanted)) 
    for n, (p, x, y) in enumerate(wanted):
        axe.annotate('', xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center',
                     arrowprops=dict(arrowstyle="->", facecolor='black',
                                     connectionstyle="arc,angleA=180,angleB=180,armA=0,armB=40,rad=15"))
        axe.annotate(p.replace('PROM_', ''), xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center')
        
        
    ## enhancer
    sub_nodes = nodes[nodes['is_enhancer']==True].sort_values(by='protein')
    sub_nodes["enha_sizes"] = enha_sizes 
#     print("all enha sub_nodes", sub_nodes.shape)
#     print()
#     print(sub_nodes[sub_nodes["SNP-binding"] != False].shape)
#     print(sub_nodes[sub_nodes["protein"].isin(snp_path_nodes)].shape)
#     print(sub_nodes[sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))) ].shape)
    
    if SNP_path or len(snp_path_nodes) > 0:
        sub_nodes = sub_nodes[(sub_nodes["SNP-binding"] != False) | (sub_nodes["protein"].isin(snp_path_nodes)) |
                              (sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]
#        print("sub_nodes", sub_nodes)                 
#        print("snp_path_nodes", snp_path_nodes)
        #sub_nodes = sub_nodes[( sub_nodes["enha_sizes"] > Enha_displayed_perc)]
        if len(snp_path_nodes) > 0:
#            print(sub_nodes.head(2))
            sub_nodes["protein_SNP"] = np.where(
                sub_nodes["SNP-binding"] != False,  sub_nodes["protein_display_name"].map(str) + "-\n" + sub_nodes["SNP-binding"].map(str)  , sub_nodes["protein_display_name"] ) 
            
            wanted = sub_nodes[['protein_SNP', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
            #print(sub_nodes[sub_nodes["SNP-binding"] != "False"])
        else:
            wanted = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord']].to_numpy()    
    else:
        sub_nodes = sub_nodes[( sub_nodes["enha_sizes"] > Enha_displayed_perc)]
        wanted = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
        
    ##print("enhancer sub_nodes", sub_nodes.shape)
    ypos = np.linspace(-1.15, 1.15, len(wanted)) 
    for n, (p, x, y) in enumerate(wanted):
        axe.annotate('', xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center',
                     arrowprops=dict(arrowstyle="->", facecolor='black',
                                     connectionstyle="arc,angleA=0,angleB=0,armA=0,armB=40,rad=15"))
        axe.annotate(p.replace('ENHA_', ''), xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center')

        
    ## intermediates
    sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')
    sub_nodes["intr_sizes"] = intr_sizes
    if SNP_path  or len(snp_path_nodes) > 0:
        #sub_nodes = sub_nodes[(sub_nodes["intr_sizes"] > int_displayed_perc) | (sub_nodes["protein"].isin(snp_path_nodes))]
        sub_nodes = sub_nodes[(sub_nodes["protein"].isin(snp_path_nodes))]
    else:
        sub_nodes = sub_nodes[(sub_nodes["intr_sizes"] > int_displayed_perc)]

    ##print("intermediates sub_nodes", sub_nodes.shape)
    wanted = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
    ypos = np.linspace(-1.15, 1.15, len(wanted)) 
    for n, (p, x, y) in enumerate(wanted):
        axe.annotate(p, xy=(x, y), xytext=[0.9 * v for v in transform(x, y)], ha="center", va='center',
                     arrowprops=dict(arrowstyle="->", facecolor='black'),
                     bbox=dict(boxstyle='round', fc='white', ec='none', alpha=0.4))

    # decorate plot
    axe.text(xprom - 0.03, 1.15, 'Promoter', size=14, ha='center')
    axe.text(xenha + 0.03, 1.15, 'Enhancer', size=14, ha='center')
    axe.text(0, 1.15, 'Intermediate', size=14, ha='center')
    axe.set_xlim(xprom - 0.2, xenha + 0.2)
    axe.set_ylim(-1.2, 1.2)
    axe.set_title(promoter, size=16, loc='left')
    axe = fig.axes[0]
    axe.add_collection(lc)
    axe.spines["right"].set_visible(False)
    axe.spines["bottom"].set_visible(False)
    axe.spines["left"].set_visible(False)
    axe.spines["top"].set_visible(False)
    axe.get_xaxis().set_visible(False)
    axe.get_yaxis().set_visible(False)
    #plt.show()
    buf = io.BytesIO() # in-memory files
    plt.savefig(buf, format = "png") # save to the above file object
    plt.close()
    data = base64.b64encode(buf.getbuffer()).decode("utf8") # encode to html elements
    return "data:image/png;base64,{}".format(data)
