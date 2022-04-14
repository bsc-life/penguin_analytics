import os

from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.gridspec as gridspec
import numpy as np
from random import random
import pandas as pd


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


def single_promoter_graph(prot, size_factor=150, top_genes=100, 
                          DBP_color='fisher_test_OR', DBP_size='degree',
                          intermediate_color='fisher_test_OR', intermediate_size='degree',
                          edge_width='druggable', edge_color='fisher_test_edge_OR',
                          DBP_cmap='coolwarm', 
                          intermediate_cmap='coolwarm',
                          edge_cmap='Greys',edge_alpha=0.5,
                          DBP_displayed_perc=75, 
                          int_displayed_perc=75, 
                          SNP_path=False, druggable=False):
    """
    :param 150 factor: scale factor (thehigher, the smaller willbe the balls)
    :param 40 top_genes: number of protein binding sites or intermediate protein names to be displayed
       in each ofthe 3 categories (promoter, middle or enhancer side)
    """
    #nodes = pd.read_csv(os.path.join(dpath, f'{prot}_nodes.tsv'), sep='\t')
    nodes = pd.read_csv("".join([dpath,'/',str(prot),'_nodes.tsv']), sep='\t')
    #edges = pd.read_csv(os.path.join(dpath, f'{prot}_edges.tsv'), sep='\t')
    edges = pd.read_csv("".join([dpath,'/',str(prot),'_edges.tsv']), sep='\t')
    edges['XYs-edge'] = edges['XYs-edge'].apply(lambda x: eval(x.replace('array', '')))
    
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
    
    DBP_displayed_perc = 1 - DBP_displayed_perc / 100
    int_displayed_perc = 1 - int_displayed_perc / 100
        
    ## promoters
    sub_nodes = nodes[nodes['is_promoter']==True].sort_values(by='protein')
    # define colors
    prom_colors = plt.cm.get_cmap(DBP_cmap)(internal_norm(sub_nodes[DBP_color]))
    # define size
    prom_sizes = internal_norm(sub_nodes[DBP_size])
    # plot promoter nodes
    axe.scatter(sub_nodes['X-graph-coord'], sub_nodes['Y-graph-coord'], s=(prom_sizes + 0.1) * size_factor,
                c=prom_colors)

    ## enhancers
    sub_nodes = nodes[nodes['is_enhancer']==True].sort_values(by='protein')
    # define colors
    enha_colors = plt.cm.get_cmap(DBP_cmap)(internal_norm(sub_nodes[DBP_color]))
    # define size
    enha_sizes = internal_norm(sub_nodes[DBP_size])
    # plot enhancer nodes
    axe.scatter(sub_nodes['X-graph-coord'], sub_nodes['Y-graph-coord'], s=(enha_sizes + 0.1) * size_factor, 
                c=enha_colors)

    ## intermediates
    sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')
    # define colors
    intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_nodes[intermediate_color]))
    # define size
    intr_sizes = internal_norm(sub_nodes[intermediate_size])
    # plot intermediate nodes
    axe.scatter(sub_nodes['X-graph-coord'], sub_nodes['Y-graph-coord'], s=(intr_sizes + 0.1) * size_factor,
                c=intr_colors)

    
    ### plot edges
    filter_edges = True
    if SNP_path and druggable:
        sub_edges =  (edges['SNP-path']) & (edges['druggable'])
    elif SNP_path:
        sub_edges =  edges['SNP-path']
    elif druggable:
        sub_edges =  edges['druggable'] == druggable
    else:
        filter_edges = False
    lc = LineCollection(edges['XYs-edge'], 
                        alpha=(sub_edges * edge_alpha) if filter_edges else edge_alpha,
                        lw=internal_norm(edges[edge_width]) * 0.8 + 0.5 if edge_width else 1, 
                        color=(plt.cm.get_cmap(edge_cmap)(internal_norm(edges[edge_color]) * 0.9 + 0.1) if
                               edge_cmap and edge_color else 'k'), zorder=-1)

    ### Annotation... this is the slowest part...
    ## promoter
    sub_nodes = nodes[nodes['is_promoter']==True].sort_values(by='protein')[prom_sizes > DBP_displayed_perc]
    wanted = sub_nodes[['protein', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
    ypos = np.linspace(-1.15, 1.15, len(wanted)) 
    for n, (p, x, y) in enumerate(wanted):
        axe.annotate('', xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center',
                     arrowprops=dict(arrowstyle="->", facecolor='black',
                                     connectionstyle="arc,angleA=180,angleB=180,armA=0,armB=40,rad=15"))
        axe.annotate(p.replace('PROM_', ''), xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center')
    ## enhancer
    sub_nodes = nodes[nodes['is_enhancer']==True].sort_values(by='protein')[enha_sizes > DBP_displayed_perc]
    wanted = sub_nodes[['protein', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
    ypos = np.linspace(-1.15, 1.15, len(wanted)) 
    for n, (p, x, y) in enumerate(wanted):
        axe.annotate('', xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center',
                     arrowprops=dict(arrowstyle="->", facecolor='black',
                                     connectionstyle="arc,angleA=0,angleB=0,armA=0,armB=40,rad=15"))
        axe.annotate(p.replace('ENHA_', ''), xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center')
    ## intermediates
    sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')[intr_sizes > int_displayed_perc]
    wanted = sub_nodes[['protein', 'X-graph-coord', 'Y-graph-coord']].to_numpy()
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
    axe.set_title(prot, size=16, loc='left')
    axe = fig.axes[0]
    axe.add_collection(lc)
    axe.spines["right"].set_visible(False)
    axe.spines["bottom"].set_visible(False)
    axe.spines["left"].set_visible(False)
    axe.spines["top"].set_visible(False)
    axe.get_xaxis().set_visible(False)
    axe.get_yaxis().set_visible(False)
    plt.show()

