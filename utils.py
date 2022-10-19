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
import csv

setupBaseDir = os.path.dirname(__file__)

dpath = os.path.join(setupBaseDir,'graph_tables')

expression_gepia  = pd.read_csv(os.path.join(setupBaseDir,"src/gepia2_logFC1_0.01_PRAD_degenes.txt"), sep='\t')[["Gene Symbol", "Log2(Fold Change)", "adjp"]]
expression_gepia.columns = [ 'Gepia_' + str(col) for col in expression_gepia.columns]

expression_FPKM  = pd.read_csv(os.path.join(setupBaseDir,"src/Cuff_Gene_Counts.txt"), sep='\t')
expression_FPKM.columns = [ 'FPKM_' + str(col) for col in expression_FPKM.columns]

expression_lhsar_lncap  = pd.read_csv(os.path.join(setupBaseDir,"src/DESeq_LNCaP_LHSAR.txt"), sep=' ')
expression_lhsar_lncap = expression_lhsar_lncap[["Gene.name","log2FoldChange", "padj"]]
expression_lhsar_lncap.columns = [ 'LNCAP_LHSAR_' + str(col) for col in expression_lhsar_lncap.columns]



node_druggability  = pd.read_csv(os.path.join(setupBaseDir,"src/node_druggability.txt"), sep='\t')
paintor_snp_in_gene = pd.read_csv(os.path.join(setupBaseDir,"src/nodes_with_paintorSNP.txt"), sep='\t')
paintor_snp_in_gene = paintor_snp_in_gene.groupby("protein").agg(
    n_rsid = pd.NamedAgg(column="rsid", aggfunc="nunique"))

proteins_22 = ["SMAD2", "KAT5", "NCOR2", "MAPK8", "SMAD4", "CREBBP", 
                "CTNNB1", "PGR", "HDAC3", "HDAC2", "GSK3B", "UBA52", "UBE2I", 
                "JUND", "PIAS1", "XRCC5", "CDK6", "XRCC6", "MAPK1", "FOS", "HIF1A", "MAPK3"]


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

def internal_norm_fisherOR(vals):
    vals = vals - 1 # i want the center to be zero
    mask = vals != np.inf
    vmax = np.nanmax(vals.loc[mask])
    vmin = np.nanmin(vals.loc[mask])
    vmax = max(vmax, abs(vmin))
    if vmax == True or vmax == False:
        vmax = 1
    vmin = -vmax
    #print(vmin, vmax)
    return np.array([v if np.isfinite(v) else vmax * 1.1 if v > 0 else vmin * 0.9
                    for v in vals]) / vmax






def single_promoter_graph(
    prot, 
    DBP_color,
    DBP_size,
    intermediate_color,
    intermediate_size,

    de_filter,
    DGE_tumorVSnormal_logFC_GEPIA,
    FPKM,
    DGE_lhsar_lncap,
    
    SNP_BindingSite_path,
    SNP_GenLoc_path,
    only_enriched_nodes,
    only_enriched_edges):

    druggable=False
    Prom_displayed_perc=100
    Enha_displayed_perc = 100
    int_displayed_perc=100
    size_factor=150
    top_genes=100
    edge_alpha = 0.8
    edge_color='fisher_test_edge_OR'
    edge_width='fisher_test_edge_pval'
    DBP_cmap = 'coolwarm'
    intermediate_cmap = 'coolwarm'
    edge_alpha = float(edge_alpha)
    """
    :param 150 factor: scale factor (thehigher, the smaller willbe the balls)
    :param 40 top_genes: number of protein binding sites or intermediate protein names to be displayed
        in each ofthe 3 categories (promoter, middle or enhancer side)
    """
    
    if prot == 'Select a promoter Gene...':
        # fig = plt.figure(figsize=(0,0), facecolor='white')
        # buf = io.BytesIO() # in-memory files
        # plt.savefig(buf, format = "png", bbox_inches='tight') # save to the above file object
        # plt.close()
        # data = base64.b64encode(buf.getbuffer()).decode("utf8") # encode to html elements
        # #return "data:image/png;base64,{}".format(data)
        return (None, 0)
    else:   
        #print(prot)
        # verbatim
        verbatim = []
        
        # read nodes
        nodes = pd.read_csv(os.path.join(dpath, f'{prot}_nodes.tsv'), sep='\t')
        nodes = nodes.merge(node_druggability, on='protein', how='left')
        nodes["Gene.name"] = nodes['protein'].str.replace("ENHA_", "").str.replace("PROM_", "")
        nodes = nodes.merge(paintor_snp_in_gene, left_on='Gene.name', right_on = "protein", how='left')
        nodes = nodes.merge(expression_gepia, left_on='Gene.name', right_on = "Gepia_Gene Symbol", how='left')
        nodes = nodes
        #nodes["Log2(Fold Change)"] = nodes["Log2(Fold Change)"].fillna(min(expression_gepia["Log2(Fold Change)"]))
        nodes = nodes.merge(expression_FPKM, left_on='Gene.name', right_on = "FPKM_Gene_ID", how='left')   
        nodes = nodes.merge(expression_lhsar_lncap, left_on='Gene.name', right_on = "LNCAP_LHSAR_Gene.name", how='left')   

        nodes = nodes.drop(['Gepia_Gene Symbol', 'FPKM_Gene_ID', 'LNCAP_LHSAR_Gene.name'], axis=1) 

        #nodes["LNCaP_1"] = nodes["LNCaP_1"].fillna(0)
        #nodes["LNCaP_2"] = nodes["LNCaP_2"].fillna(0)

        nodes["in_22"] = nodes["Gene.name"].isin(proteins_22)

        #paintor_SNP is if there are Paintor snps in the node that fall in its genomic location
        nodes["paintor_SNP"] = np.where(nodes["n_rsid"] > 0  ,True , False)
        nodes["Node_druggability"] = nodes["Node_druggability"].fillna(False)


        #SNP-binding is if there is a snp overlaping the binding site
        #nodes = nodes.replace({'SNP-binding': {'True': True, 'False': False}})
        

        nodes["shape"] = np.where(
                    nodes["Node_druggability"] == True, "^", "o" ) 

        # paintor_SNP_nodes are the nodes intermediates with paintor snp in their genomic loc
        paintor_SNP_nodes = list(set(nodes[(nodes["paintor_SNP"] == True) & (nodes["is_intermediate"] == True)]["protein"].values.tolist()))

        #print("paintor_SNP_nodes", paintor_SNP_nodes)


        # read edges
        all_edges = pd.read_csv(os.path.join(dpath, f'{prot}_edges.tsv'), sep='\t')
        ##print("total unfiltered edges: ", all_edges.shape)  
        
        all_edges['XYs-edge'] = all_edges['XYs-edge'].apply(lambda x: eval(x.replace('array', '')))
        all_edges[['node1', 'node2']] = all_edges['edge'].str.split(';', expand=True)
        
                
        all_edges["SNP.intermediate.path"] = np.where((all_edges["node1"].isin(paintor_SNP_nodes)) | (all_edges["node2"].isin(paintor_SNP_nodes)), True, False)
        
    ###### DOWNLOADS
        #prepare Downloadable table 
        downloadable = all_edges.drop(["XYs-edge"], axis=1).merge(
            nodes.drop(["X-graph-coord","Y-graph-coord", 
                                        "Degree_in_cluster","Degree_out_cluster","Degree_global"  ,
                                        "Betweenness_in_cluster", "Betweenness_out_cluster" ,"Betweenness_global",
                                        "in_22", "shape"], axis=1) , left_on = "node1", right_on = "protein").merge(
            nodes.drop(["X-graph-coord","Y-graph-coord", 
                                            "Degree_in_cluster","Degree_out_cluster","Degree_global"  ,
                                            "Betweenness_in_cluster", "Betweenness_out_cluster" ,"Betweenness_global",
                                            "in_22", "shape"], axis=1) , left_on = "node2", right_on = "protein").drop(
            ["protein_x", "protein_y", 'promoter_x', 'promoter_y'] , axis=1)
        
        downloadable.columns = [i.replace("-", "_") for i in downloadable.columns]
        downloadable["edge"] = downloadable["edge"].str.replace(';',':')
    #        print(downloadable.head(1))
        #file_name = f'{random.randint(0,5000)}.txt'
    #         file_name = f'{prot}_merge.txt'
    #         if not os.exist(file_name):
    #             downloadable.to_csv(file_name, sep="\t", header=True, index=False,quoting=csv.QUOTE_NONE, quotechar='')
        
        ###### ###### ###### 

        edges = pd.DataFrame(columns = all_edges.columns)
        
        nodes["grey"] = "normal"
        if de_filter == "Gepia":
            #print("DGE_tumorVSnormal_logFC_GEPIA", DGE_tumorVSnormal_logFC_GEPIA)
            nodes = nodes[(nodes["Gepia_Log2(Fold Change)"].abs() >= DGE_tumorVSnormal_logFC_GEPIA ) | (nodes["Gepia_Log2(Fold Change)"].isnull() )]
            nodes["grey"] = np.where(nodes["Gepia_Log2(Fold Change)"].isnull() , "grey", "normal")
            #print("nodes after GEPIA", nodes.shape)
        
        elif de_filter == "FPKM":            
        
#            nodes = nodes[(nodes["FPKM_LNCaP_1"] >= FPKM) & (nodes["FPKM_LNCaP_2"] >= FPKM)]
            nodes = nodes[((nodes["FPKM_LNCaP_1"] >= FPKM) & (nodes["FPKM_LNCaP_2"] >= FPKM)) | (nodes["FPKM_LNCaP_1"].isnull() )]
            nodes["grey"] = np.where(nodes["FPKM_LNCaP_1"].isnull() , "grey", "normal")
            #print("nodes after FPKM", nodes.shape)
            
        elif de_filter == 'LHSAR VS LNCaP':            
#            nodes = nodes[nodes["LNCAP_LHSAR_log2FoldChange"].abs() >= DGE_lhsar_lncap]
            nodes = nodes[(nodes["LNCAP_LHSAR_log2FoldChange"].abs() >= DGE_lhsar_lncap ) | (nodes["LNCAP_LHSAR_log2FoldChange"].isnull())]
            nodes["grey"] = np.where(nodes["LNCAP_LHSAR_log2FoldChange"].isnull() , "grey", "normal")
            
            #print("nodes after LHSAR VS LNCaP", nodes.shape)
    

        if only_enriched_nodes:
            #remove the non enriched nodes
            nodes =  nodes[(nodes["fisher_test_OR"] > 1) & (nodes["fisher_test_pval"] < 0.01)]
            #remove the edges attaching the non enriched nodes        
    #        print("nodes after only_enriched_nodes", nodes.shape[1])

        #print("nodes after only_enriched_nodes", nodes.shape)
        #print("all_edges before", all_edges.shape)    

        all_edges = all_edges[(all_edges["node1"].isin(nodes["protein"])) & (all_edges["node2"].isin(nodes["protein"]))]
        #print("all_edges list", edges.drop_duplicates())
        #print("all_edges after node filtering", all_edges.shape)    

        # SNP_GenLoc_path are the edges that have in one or in 
        # the other end an intermediate with Paintor snp in its genomic location
        if SNP_GenLoc_path:
            edges = pd.concat([edges,
                                all_edges[all_edges["SNP.intermediate.path"] == True ]])

        # SNP-path is if the edge is part of the subnetwork that start from SPN in binding site
        if SNP_BindingSite_path:
            edges = pd.concat([edges,
                                all_edges[all_edges['SNP-path'] == True]])

        # if edge is druggable
        if druggable:
            edges = pd.concat([edges,
                                all_edges[all_edges['druggable'] == True]])


        if druggable == False and SNP_BindingSite_path == False and only_enriched_edges == False and SNP_GenLoc_path == False  :
            edges = all_edges

        # if i select the enricjed edges this will be a filter in the previous selection
        if only_enriched_edges:
            if edges.empty:
                edges = edges
                #print("edges.empty:", edges.empty)
                #edges = all_edges[(all_edges["fisher_test_edge_OR"] > 1) & (all_edges["fisher_test_edge_pvalue"] < 0.01)]
            else:
                edges = edges[(edges["fisher_test_edge_OR"] > 1) & (edges["fisher_test_edge_pvalue"] < 0.01)]

        #   print("total filtered edges: ", edges.shape)
        edges = edges.loc[edges.astype(str).drop_duplicates().index]
        #print("all_edges after  filtering", edges.shape) 
        #print("edge list", edges["edge"].drop_duplicates())

        #the nodes in the SNP path (nodes that are in the subnetwork that start from SPN in binding site
        # PLUS the paintor SNP nodes (in hte genomic loc) if i select them
        #snp_path_nodes = list(set(all_edges[all_edges["SNP-path"] == True][["node1", "node2"]].stack().values.tolist()))
        #SNP_intermediate_path_nodes = list(set(all_edges[all_edges["SNP.intermediate.path"] == True][["node1", "node2"]].stack().values.tolist()))
        #print("SNP_intermediate_path_nodes", SNP_intermediate_path_nodes)
        nodes["protein_display_name"] = np.where(
                nodes["paintor_SNP"] == True, nodes["protein"] + "*",  nodes["protein"]) 

            
        ##################################
        # start plotting    
        fig = plt.figure(figsize=(14, 20), facecolor='white')
        gs0 = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[2, 1], hspace=0.1)
        axe = fig.add_subplot(gs0[0])

        #################
        verbatim.append(["Number of edges: ", edges.shape[0]])
        # X limits
        xprom = -1.1
        xenha = 1.1
        
        Prom_displayed_perc = 1 - Prom_displayed_perc / 100
        Enha_displayed_perc = 1 - Enha_displayed_perc / 100
        int_displayed_perc = 1 - int_displayed_perc / 100
            
        ######    ######
        ## promoters
        sub_nodes = nodes[nodes['is_promoter']==True].sort_values(by='protein')
        # define size
        prom_sizes = internal_norm(sub_nodes[DBP_size]) if sub_nodes.size != 0  else None    

        sub_nodes["prom_sizes"] = prom_sizes 
        if SNP_BindingSite_path or SNP_GenLoc_path:
            #sub_nodes = sub_nodes[(sub_nodes["prom_sizes"]  > Prom_displayed_perc ) | (sub_nodes["protein"].isin(snp_path_nodes))]
            sub_nodes = sub_nodes[(sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]
        else:
            sub_nodes = sub_nodes[(sub_nodes["prom_sizes"]  > Prom_displayed_perc )]
        #print("promoters sub_nodes", sub_nodes)
        #    sub_nodes = sub_nodes[sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist())))]
        ##print("prototer sub_nodes", sub_nodes.shape)
        wanted_promoters = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord', 'in_22']].to_numpy()

        #print("Number of Promoter Binders: ", wanted_promoters.shape[0])

        # plot promoter nodes
        for marker in ["o", "^"]:
            sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker]
            if len(sub_sub_nodes) == 0:
                continue
            # define size
            sub_sizes = internal_norm(sub_sub_nodes[DBP_size])
            # define colors
            
            #the normal coloured nodes            
            prom_colors = plt.cm.get_cmap(DBP_cmap)(internal_norm(sub_sub_nodes[DBP_color]))            
            prom_colors[np.where(sub_sub_nodes["grey"] == "grey") ] = mcolors.to_rgba("grey")

            axe.scatter(sub_sub_nodes['X-graph-coord'], 
                        sub_sub_nodes['Y-graph-coord'],
                        s=(sub_sizes + 0.1) * size_factor,
                        marker = marker,
                        c=prom_colors)

            


        ######    ######
        ###### enhancers
        sub_nodes = nodes[nodes['is_enhancer']==True].sort_values(by='protein')
        # define size
        #print("enhancer sub_nodes", sub_nodes["protein_display_name"])
        enha_sizes = internal_norm(sub_nodes[DBP_size]) if sub_nodes.size != 0  else None
    #    enha_sizes = internal_norm(sub_nodes[DBP_size])

        sub_nodes["enha_sizes"] = enha_sizes  
        if SNP_BindingSite_path or SNP_GenLoc_path:
            sub_nodes = sub_nodes[(sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]
    #         sub_nodes = sub_nodes[(sub_nodes["SNP-binding"] != False) | (sub_nodes["protein"].isin(snp_path_nodes)) |
    #                               (sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]

            if SNP_BindingSite_path:
    #            print(sub_nodes.head(2))
                sub_nodes["protein_SNP"] = np.where(
                    sub_nodes["SNP-binding"] != False,  sub_nodes["protein_display_name"].map(str) + "-\n" + sub_nodes["SNP-binding"].map(str)  , sub_nodes["protein_display_name"] ) 

                wanted_enhancers = sub_nodes[['protein_SNP', 'X-graph-coord', 'Y-graph-coord','in_22']].to_numpy()
                #print(sub_nodes[sub_nodes["SNP-binding"] != "False"])
            else:
                wanted_enhancers = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord','in_22']].to_numpy()    
        else:
            sub_nodes = sub_nodes[( sub_nodes["enha_sizes"] > Enha_displayed_perc)]
            wanted_enhancers = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord','in_22']].to_numpy()


        #print("Number of Enhancer Binders: ", wanted_enhancers.shape[0])
        # plot enhancer nodes
        for marker in ["o", "^"]:
            sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker]
            if len(sub_sub_nodes) == 0:
                continue
            # define size    
            sub_sizes = internal_norm(sub_sub_nodes[DBP_size])

            #the normal coloured nodes
            enha_colors = plt.cm.get_cmap(DBP_cmap)(internal_norm(sub_sub_nodes[DBP_color]))            
            enha_colors[np.where(sub_sub_nodes["grey"] == "grey") ] = mcolors.to_rgba("grey")

            axe.scatter(sub_sub_nodes['X-graph-coord'], 
                        sub_sub_nodes['Y-graph-coord'],
                        s=(sub_sizes + 0.1) * size_factor,
                        marker = marker,
                        c=enha_colors)



        ##############
        ######## intermediates
        sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')
        
        # define size
        if sub_nodes.size != 0:
            intr_sizes = internal_norm(sub_nodes[DBP_size])
            sub_nodes["intr_sizes"] = intr_sizes

            if intermediate_color in ["fisher_test_OR" ] and only_enriched_nodes == False:
                sub_nodes["binary_sub_nodes"] = np.where((sub_nodes["fisher_test_OR"] > 1) & (sub_nodes["fisher_test_pval"] < 0.01 ), 1, 0)
                intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_nodes["binary_sub_nodes"]))
            else:
                intr_colors = plt.cm.get_cmap(intermediate_cmap)(internal_norm(sub_nodes[intermediate_color]))
            intr_colors[np.where(sub_nodes["grey"] == "grey") ] = mcolors.to_rgba("grey")

        else:
            intr_colors = None
            intr_sizes = None



        for marker in ["o", "^"]:
            sub_sub_nodes = sub_nodes[sub_nodes["shape"] == marker] 

            if len(sub_sub_nodes) == 0:
                continue

            #plot two markers, The normalisation has to be done before.
            # take the indeces
            sub_sub_nodes_index = (sub_nodes["shape"] == marker)

            #extract the colours
            sub_intr_colors = intr_colors[sub_sub_nodes_index]

            sub_sub_sizes = intr_sizes[sub_sub_nodes_index,]   


            axe.scatter(sub_sub_nodes['X-graph-coord'],
                        sub_sub_nodes['Y-graph-coord'], 
                        s=(sub_sub_sizes + 0.1) * size_factor,
                        marker = marker,
                        c=sub_intr_colors)

        ### plot edges
        #sub_edges = "empty"
        #if edges.equals(all_edges) or edges.empty:
        #    filter_edges = False
        #else:
        #    filter_edges = True
            #   sub_edges = [True] * len(edges.index)
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

        if SNP_GenLoc_path and SNP_BindingSite_path:
            conditions = [
                ((edges["node1"].isin(paintor_SNP_nodes)) | (edges["node2"].isin(paintor_SNP_nodes))) & (edges['SNP-path'] == True),
                (edges['SNP-path'] == True),
                ((edges["node1"].isin(paintor_SNP_nodes)) | (edges["node2"].isin(paintor_SNP_nodes)))
            ]
            choices = ['purple', 'orange', 'darkgreen']
            v = mcolors.to_rgba("black")
            edges["colour"] = np.select(conditions, choices, default = "black")

        elif SNP_GenLoc_path:
            conditions = [
                ((edges["node1"].isin(paintor_SNP_nodes)) | (edges["node2"].isin(paintor_SNP_nodes)))
            ]
            choices = ['darkgreen']
            v = mcolors.to_rgba("black")
            edges["colour"] = np.select(conditions, choices, default = "black")        

        elif SNP_BindingSite_path:
            conditions = [(edges['SNP-path'] == True)]
            choices = [ 'orange']
            v = mcolors.to_rgba("black")
            edges["colour"] = np.select(conditions, choices, default = "black")

        else:
            edges["colour"] = "black"

        lc = LineCollection(edges['XYs-edge'], 
                        alpha = 0.7,
                        lw = 0.7,
                        color = edges["colour"], 
                        zorder=-1)


        ### Annotation... this is the slowest part...
        ## promoter  
        ypos = np.linspace(-1.15, 1.15, len(wanted_promoters)) 
        for n, (p, x, y, in_22) in enumerate(wanted_promoters):
            axe.annotate('', xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center',
                         arrowprops=dict(arrowstyle="->", facecolor='black',
                                         connectionstyle="arc,angleA=180,angleB=180,armA=0,armB=40,rad=15"))
            if in_22:
                axe.annotate(p.replace('PROM_', ''), xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center',
                             weight='bold', fontsize = 12)  
            else:
                axe.annotate(p.replace('PROM_', ''), xy=(x, y), xytext=(x - 0.2, ypos[n]), ha="right", va='center')



        ## enhancer
        ##print("enhancer sub_nodes", sub_nodes.shape)
        ypos = np.linspace(-1.15, 1.15, len(wanted_enhancers)) 
        for n, (p, x, y, in_22) in enumerate(wanted_enhancers):
            axe.annotate('', xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center',
                         arrowprops=dict(arrowstyle="->", facecolor='black',
                                         connectionstyle="arc,angleA=0,angleB=0,armA=0,armB=40,rad=15"))
            if in_22:
                axe.annotate(p.replace('ENHA_', ''), xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center',
                             weight='bold', fontsize = 12)
            else:
                axe.annotate(p.replace('ENHA_', ''), xy=(x, y), xytext=(x + 0.2, ypos[n]), ha="left", va='center')



        ## intermediates
        sub_nodes = nodes[nodes['is_intermediate']==True].sort_values(by='protein')
        sub_nodes["intr_sizes"] = intr_sizes
        if SNP_BindingSite_path or SNP_GenLoc_path:
            sub_nodes = sub_nodes[(sub_nodes["protein"].isin(list(set(edges[["node1", "node2"]].stack().values.tolist()))))]

            #sub_nodes = sub_nodes[(sub_nodes["intr_sizes"] > int_displayed_perc) | (sub_nodes["protein"].isin(snp_path_nodes))]
            #sub_nodes = sub_nodes[(sub_nodes["protein"].isin(snp_path_nodes))]
        else:
            sub_nodes = sub_nodes[(sub_nodes["intr_sizes"] > int_displayed_perc)]

        #print("intermediates sub_nodes", sub_nodes["protein_display_name"])

        ##print("intermediates sub_nodes", sub_nodes.shape)
        wanted = sub_nodes[['protein_display_name', 'X-graph-coord', 'Y-graph-coord','in_22']].to_numpy()
        #print("Number of Intermediate nodes: ", wanted.shape[0])
        ypos = np.linspace(-1.15, 1.15, len(wanted)) 
        for n, (p, x, y,in_22) in enumerate(wanted):
            if in_22:
                axe.annotate(p, xy=(x, y), xytext=[0.9 * v for v in transform(x, y)], ha="center", va='center',
                             arrowprops=dict(arrowstyle="->", facecolor='black'),
                             weight='bold', fontsize = 12,
                             bbox=dict(boxstyle='round', fc='white', ec='none', alpha=0.4))
            else:
                axe.annotate(p, xy=(x, y), xytext=[0.9 * v for v in transform(x, y)], ha="center", va='center',
                             arrowprops=dict(arrowstyle="->", facecolor='black'),
                             bbox=dict(boxstyle='round', fc='white', ec='none', alpha=0.4))

        # decorate plot
        axe.text(xprom - 0.03, 1.15, 'Promoter', size=14, ha='center')
        axe.text(xenha + 0.03, 1.15, 'Enhancer(s)', size=14, ha='center')
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
            
        buf = io.BytesIO() # in-memory files
        plt.savefig(buf, format = "png", bbox_inches='tight') # save to the above file object
        plt.close()
        data = base64.b64encode(buf.getbuffer()).decode("utf8") # encode to html elements



        
        with open(os.path.join(setupBaseDir,"outfiles/network_stats.txt"), "w") as f:
            wr = csv.writer(f, delimiter =' ')
            wr.writerows(verbatim)

        downloadable.to_csv(os.path.join(setupBaseDir,f"outfiles/{prot}_network.tsv"),sep='\t',header=True, index=False)
        #downloadable.to_csv("outfiles/network.tsv",sep='\t',header=True, index=False)
        #button for all merged data
        # html_buttons = '''<html>
        #     <head>
        #     <meta name="viewport" content="width=device-width, initial-scale=1">
        #     </head>
        #     <body>
        #     <a download="{filename}" href="data:text/csv;base64,{payload}" download>
        #     <button style="width: 250px;" class="p-Widget jupyter-widgets jupyter-button widget-button mod-warning">{label}</button>
        #     </a>
        #     </body>
        #     </html>
        #     '''
        # b64 = base64.b64encode(downloadable.to_csv(index=False).encode())
        # payload = b64.decode()
        # html_download_data = html_buttons.format(payload=payload, filename=f'{prot}.csv', label = "Download Data")
        #display(HTML(html_download_data))


        return ("data:image/png;base64,{}".format(data), 1)
