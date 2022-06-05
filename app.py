from dash import Dash, html, dcc, html, Input, Output
import dash_daq as daq
import os
import utils

app = Dash(__name__)
server = app.server
app.title = "PENGUIN"

dpath = 'graph_tables'
prots = sorted([f[:-10] for f in os.listdir(dpath) if f.endswith('_nodes.tsv')])
node_stats_short = ['degree', 'expression_log2FoldChange', 'expression_padj', 'SNP-binding', 'fisher_test_OR', 'fisher_test_pval', 'Node_druggability']

node_stats = [
    'degree', 
    'expression_log2FoldChange', 'expression_padj', 
    'SNP-binding',
    'Betweenness_global', 'Degree_global', 'Betweenness_in_cluster',
    'Betweenness_out_cluster', 'Betweenness_enrichment',
    'Degree_in_cluster', 'Degree_out_cluster', 'Degree_enrichment',
    'Betweenness_enrichment_pval', 'Degree_enrichment_pval', 'fisher_test_OR', 'fisher_test_pval',
    'Node_druggability'
]
edge_stats = [
    'None',
    'fisher_test_edge_OR', 
    'fisher_test_edge_pvalue',
    'SNP-path', 
    'druggable'
]

DBP_cmap=['None', 'coolwarm', 'Greys', 'Reds', 'viridis']

app.layout = html.Div(
            className="content",
            children=[

                html.Div(
                    className="left_menu",
                    children=[
                        html.Div(
                            [
                                
                                html.Label("Promoter"),
                                dcc.Dropdown(prots,'FOXA1',
                                    id='promoter'
                                    ),

                                html.Label("Size factor"),
                                dcc.Slider(0, 1000, 5,
                                    value=150,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='size_factor'
                                    ),

                                html.Label("Top genes"),
                                dcc.Slider(0, 300, 1,
                                    value=100,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='top_genes'
                                    ),

                                html.Label("DBPs color"),
                                dcc.Dropdown(node_stats_short,'fisher_test_OR',
                                    id='DBP_color'
                                    ),

                                html.Label("DBPs size"),
                                dcc.Dropdown(node_stats_short,'degree',
                                    id='DBP_size'
                                    ),

                                html.Label("Intermediates color"),
                                dcc.Dropdown(node_stats,'fisher_test_OR',
                                    id='intermediate_color'
                                    ),
                                html.Label("Intermediates size"),
                                dcc.Dropdown(node_stats,'Betweenness_enrichment',
                                    id='intermediate_size'
                                    ),

                                html.Label("Edges color"),
                                dcc.Dropdown(edge_stats,'fisher_test_edge_OR',
                                    id='edge_color'
                                    ),

                                html.Label("Edges width"),
                                dcc.Dropdown(edge_stats,'druggable',
                                    id='edge_width'
                                    ),

                                html.Label("DBPs CMAP"),
                                dcc.Dropdown(DBP_cmap,'coolwarm',
                                    id='DBP_cmap'
                                    ),

                                html.Label("Intermediates CMAP"),
                                dcc.Dropdown(DBP_cmap,'coolwarm',
                                    id='intermediate_cmap'
                                    ),

                                html.Label("Edges CMAP"),
                                dcc.Dropdown(DBP_cmap, 'None',
                                    id='edge_cmap'
                                    ),

                                html.Label("Edges alpha"),
                                dcc.Slider(0, 1, 0.1,
                                    value=0.5,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='edge_alpha'
                                    ),

                                html.Label("Percentage of displayed promoter"),
                                dcc.Slider(0, 100, 1,
                                    value=75,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='Prom_displayed_perc'
                                    ),

                                html.Label("Percentage of displayed enhancer"),
                                dcc.Slider(0, 100, 1,
                                    value=75,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='Enha_displayed_perc'
                                    ),

                                html.Label("Percentage of displayed intermediates"),
                                dcc.Slider(0, 100, 1,
                                    value=75,
                                    marks=None,
                                    tooltip={"placement": "bottom", "always_visible": True},
                                    id='int_displayed_perc'
                                    ),

                                html.Label("SNP_path"),
                                daq.ToggleSwitch(value=False,
                                    id='SNP_path'),

                                html.Label("druggable"),
                                daq.ToggleSwitch(value=True,
                                    id='druggable'),

                                html.Label("paintor_SNP"),
                                daq.ToggleSwitch(value=False,
                                    id='paintor_SNP'),

                                html.Label("only_enriched_nodes"),
                                daq.ToggleSwitch(value=False,
                                    id='only_enriched_nodes'),

                                html.Label("enriched_edges"),
                                daq.ToggleSwitch(value=False,
                                    id='enriched_edges'),

                                html.Label("paintor_SNP_edges"),
                                daq.ToggleSwitch(value=False,
                                    id='paintor_SNP_edges'),
                            ]
                        ),
                ]),

                html.Div(
                    className="right_content",
                    children=[
                        html.Div(
                            className="top_metrics",
                            children=[
                                html.H2('PENGUIN (Promoter-ENancher-GUided Interaction Networks)')
                            ]
                        ),
                        
                        html.Img(id='display-value', style={'height':'65%', 'width':'65%'})
                    
                    ], style={'textAlign': 'center'}

                )
            ])

@app.callback(
    Output('display-value', 'src'),
    Input('promoter', 'value'),
    Input('size_factor', 'value'),
    Input('top_genes', 'value'),
    Input('DBP_color', 'value'),
    Input('DBP_size', 'value'),
    Input('intermediate_color', 'value'),
    Input('intermediate_size', 'value'),
    Input('edge_color', 'value'),
    Input('edge_width', 'value'),
    Input('DBP_cmap', 'value'),
    Input('intermediate_cmap', 'value'),
    Input('edge_cmap', 'value'),
    Input('edge_alpha', 'value'),
    Input('Prom_displayed_perc', 'value'),
    Input('Enha_displayed_perc', 'value'),
    Input('int_displayed_perc', 'value'),
    Input('SNP_path', 'value'),
    Input('druggable', 'value'),
    Input('paintor_SNP', 'value'),
    Input('only_enriched_nodes', 'value'),
    Input('enriched_edges', 'value'),
    Input('paintor_SNP_edges', 'value')
)

def display_value(promoter, size_factor, top_genes, DBP_color, DBP_size,
                    intermediate_color, intermediate_size, edge_color, edge_width,
                    DBP_cmap, intermediate_cmap, edge_cmap, edge_alpha,
                    Prom_displayed_perc, Enha_displayed_perc, int_displayed_perc,
                    SNP_path, druggable, paintor_SNP, only_enriched_nodes, enriched_edges, paintor_SNP_edges):
    p = utils.single_promoter_graph(promoter, size_factor, top_genes, DBP_color, DBP_size,
                    intermediate_color, intermediate_size, edge_color, edge_width,
                    DBP_cmap, intermediate_cmap, edge_cmap, edge_alpha,
                    Prom_displayed_perc, Enha_displayed_perc, int_displayed_perc,
                    SNP_path, druggable, paintor_SNP, only_enriched_nodes, enriched_edges, paintor_SNP_edges)
    return(p)

if __name__ == '__main__':
    app.run_server(debug=True)
