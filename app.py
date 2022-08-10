
# setup
import dash
from dash import Dash, html, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_daq as daq
import os
import utils
import pandas as pd
from settings import config

# App Instance
app = dash.Dash(name=config.app_name, assets_folder="src", external_stylesheets=[dbc.themes.LUX, config.fontawesome])
app.title = config.app_name
server = app.server

########################## Navbar ##########################

navbar = dbc.Nav(className="nav nav-pills", children=[
    ## logo/home
    dbc.NavItem(html.Img(src=app.get_asset_url("penguin.png"), height="40px")),
    ## about
    dbc.NavItem(html.Div([
        dbc.NavLink("About", href="/", id="about-popover", active=False),
        dbc.Popover(id="about", is_open=False, target="about-popover", children=[
            #dbc.PopoverHeader("How it works"), 
            dbc.PopoverBody(config.about)
        ])
    ])),
    ## links
    dbc.DropdownMenu(label="Links", nav=True, children=[
    	dbc.DropdownMenuItem([html.I(className="fa fa-github"), "  Code"], href=config.code, target="_blank"),
        dbc.DropdownMenuItem([html.I(className="fa fa-file-text-o"), "  Tutorial"], href=config.tutorial, target="_blank")
    ])
])

@app.callback(output=[Output(component_id="about", component_property="is_open"), 
                      Output(component_id="about-popover", component_property="active")], 
              inputs=[Input(component_id="about-popover", component_property="n_clicks")], 
              state=[State("about","is_open"), State("about-popover","active")])
def about_popover(n, is_open, active):
    if n:
        return not is_open, active
    return is_open, active


########################## Body ##########################

dpath = 'graph_tables'
prots = sorted([f[:-10] for f in os.listdir(dpath) if f.endswith('_nodes.tsv')])
                               
node_stats = [
    'degree', 
    'expression_log2FoldChange', 'expression_padj', 
    'Betweenness_global', 'Degree_global', 'Betweenness_in_cluster',
    'Betweenness_out_cluster', 'Betweenness_enrichment',
    'Degree_in_cluster', 'Degree_out_cluster', 'Degree_enrichment',
    'Betweenness_enrichment_pval', 'Degree_enrichment_pval', 'fisher_test_OR', 'fisher_test_pval',
    'Node_druggability'
]
edge_stats = [
    None,
    'fisher_test_edge_OR', 
    'fisher_test_edge_pvalue',
    'SNP-path', 
    'druggable'
]

options = [s for s in node_stats if not 'enrich' in s and not 'cluster' in s and not 'global' in s]

# Output
body = dbc.Row([

	dbc.Col(width=2, children=[

		
			dbc.Label("Promoter"),
		    dcc.Dropdown(prots,'MYC',
		        id='promoter'
		        ),

		    dbc.Label("DBPs displayed size"),
		    dcc.Dropdown(options,'fisher_test_pval',
		        id='DBP_size'
		        ),

		    dbc.Label("DBPs displayed color"),
		    dcc.Dropdown(options,'fisher_test_OR',
		        id='DBP_color'
		        ),

		    dbc.Label("Intermediates displayed size"),
		    dcc.Dropdown(node_stats,'fisher_test_pval',
		        id='intermediate_size'
		        ),

		    dbc.Label("Intermediates displayed color"),
		    dcc.Dropdown(node_stats,'fisher_test_OR',
		        id='intermediate_color'
		        ),
		    
		    html.Br(),
		    
		    dbc.Label("Differential expression cutoff (GEPIA) [log2FC]"),
		    dcc.Slider(-5, 8, 0.5,
		        value=-5,
		        marks=None,
		        tooltip={"placement": "bottom", "always_visible": True},
		        id='DGE_tumorVSnormal_logFC_GEPIA_input'
		        ),

		    dbc.Label("Differential expression cutoff (LNCaP/LHSAR) [FPKM]"),
		    dcc.Slider(0, 100, 2,
		        value=0,
		        marks=None,
		        tooltip={"placement": "bottom", "always_visible": True},
		        id='FPKM_input'
		        ),
            
		    html.Br(),
		    
		    dbc.Label("SNP paths (PrCa SNPs in enhancer binding motifs)"),
		    daq.ToggleSwitch(value=True,
		        id='SNP_BindingSite_path_input'),

		    dbc.Label("SNP paths (PrCa SNPs in intermediate proteins)"),
		    daq.ToggleSwitch(value=True,
		        id='SNP_GenLoc_path_input'),

		    dbc.Label("Nodes enriched in cluster GWAS+"),
		    daq.ToggleSwitch(value=True,
		        id='only_enriched_nodes_input'),

		    dbc.Label("Edges enriched in cluster GWAS+"),
		    daq.ToggleSwitch(value=True,
		        id='only_enriched_edges_input')

	], style={"fontSize": "10px"}),

	dbc.Col(width=10, children=[
			dbc.Spinner([
				html.Img(id='display-value', style={'height':'83%', 'width':'83%'})
			], color="primary", type="grow"),
	])
])


@app.callback(
	Output('display-value', 'src'),
	Input('promoter', 'value'),
	Input('DBP_size', 'value'),
	Input('DBP_color', 'value'),
	Input('intermediate_size', 'value'),
	Input('intermediate_color', 'value'),
	Input('DGE_tumorVSnormal_logFC_GEPIA_input', 'value'),
	Input('FPKM_input', 'value'),
	Input('SNP_BindingSite_path_input', 'value'),
	Input('SNP_GenLoc_path_input', 'value'),
	Input('only_enriched_nodes_input', 'value'),
	Input('only_enriched_edges_input', 'value')
)
def display_value(promoter, 
                          DBP_size, DBP_color,
                          intermediate_size, intermediate_color,
                          DGE_tumorVSnormal_logFC_GEPIA_input,
                          FPKM_input,
                          SNP_BindingSite_path_input,
                          SNP_GenLoc_path_input,
                          only_enriched_nodes_input,
                          only_enriched_edges_input):
    p = utils.single_promoter_graph(promoter, 
                          DBP_size, DBP_color,
                          intermediate_size, intermediate_color,
                          DGE_tumorVSnormal_logFC_GEPIA_input,
                          FPKM_input,
                          SNP_BindingSite_path_input,
                          SNP_GenLoc_path_input,
                          only_enriched_nodes_input,
                          only_enriched_edges_input)
    return(p)

@app.callback(
    Output("download_network_tsv", "data"),
    Input("network_tsv", "n_clicks"),
    prevent_initial_call=True
)
def func(n_clicks):
    return dcc.send_file(
        "outfiles/network.tsv"
    )

########################## App Layout ##########################
app.layout = dbc.Container(fluid=True, children=[
    html.Br(),
    html.H1(config.app_name, id="nav-pills"),
    html.H4('Promoter-ENancher-GUided Interaction Networks'),
    navbar,
    html.Br(),
    html.Div(
    [
        html.Button("Download the network", id="network_tsv"),
        dcc.Download(id="download_network_tsv"),
    ]),
    html.Br(),
    body
])



########################## Run ##########################
if __name__ == "__main__":
    debug = True if config.ENV == "DEV" else False
    app.run_server(debug=debug, host=config.host, port=config.port)
