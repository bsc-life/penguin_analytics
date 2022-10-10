
# setup
import dash
from dash import Dash, html, dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import dash_daq as daq

from dash.dependencies import State

from dash_extensions import Download
from dash_extensions.snippets import send_bytes

import os
import sys
setupBaseDir = os.path.dirname(__file__)
sys.path.insert(0, setupBaseDir)
import utils
import pandas as pd
from settings import config

# App Instance
app = dash.Dash(name=config.app_name, assets_folder=os.path.join(setupBaseDir, "src"), external_stylesheets=[dbc.themes.LUX, config.fontawesome])
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
    dbc.DropdownMenu(label="Code & Tutorial", nav=True, children=[
    	dbc.DropdownMenuItem([html.I(className="fa fa-github"), "  Code"], href=config.code, target="_blank"),
        #dbc.DropdownMenuItem([html.I(className="fa fa-file-text-o"), "  Tutorial"], href=app.get_asset_url("documentation.html"), target="_blank")
        dbc.DropdownMenuItem(html.A("  Tutorial",className="fa fa-file-text-o",href=app.get_asset_url("documentation.html"), target="_blank"))
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

dpath = os.path.join(setupBaseDir,'graph_tables')
prots = sorted([f[:-10] for f in os.listdir(dpath) if f.endswith('_nodes.tsv')])
#prots.append('Select a promoter Gene...')

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
            html.Br(),
		    dbc.Label("DBPs displayed size"),
		    dcc.Dropdown(options,'fisher_test_pval',
		        id='DBP_size'
		        ),
            html.Br(),
		    dbc.Label("DBPs displayed color"),
		    dcc.Dropdown(options,'fisher_test_OR',
		        id='DBP_color'
		        ),
            html.Br(),
		    dbc.Label("Intermediates displayed size"),
		    dcc.Dropdown(node_stats,'fisher_test_pval',
		        id='intermediate_size'
		        ),
            html.Br(),
		    dbc.Label("Intermediates displayed color"),
		    dcc.Dropdown(node_stats,'fisher_test_OR',
		        id='intermediate_color'
		        ),
            html.Br(),		    
		    
		    
		    dbc.Label("Apply expression filtering on nodes:"),            
            dcc.Dropdown(
                id = 'de_filter',
                options=[                    
                    {'label': 'No filter', 'value': 'None'},
                    {'label': 'DEG (Gepia)', 'value': 'Gepia'},
                    {'label': 'DEG (LHSAR vs LNCaP)', 'value': 'LHSAR VS LNCaP'},
                    {'label': 'Expression cutoff (LNCaP)', 'value': 'FPKM'}                    
                ],
                value = 'None'
            ),
            html.Br(),
            # Create Div to place a conditionally visible element inside
            html.Div([                
                dbc.Label("Differential expression cutoff (GEPIA) [|log2FC|]"),
                dcc.Slider(0, 8, 0.5,
                    value=0,
                    marks=None,
                    tooltip={"placement": "bottom", "always_visible": True},
                    id='DGE_tumorVSnormal_logFC_GEPIA_input'
                    )], 
                id = 'element-to-hide1',
                style= {'display': 'block'} ),

            html.Div([
                dbc.Label("Expression cutoff (LNCaP) [FPKM]"),
                dcc.Slider(0, 100, 2,
                    value=0,
                    marks=None,
                    tooltip={"placement": "bottom", "always_visible": True},
                    id='FPKM_input'
                    )], 
                id = 'element-to-hide2',
                style= {'display': 'block'} ),

            html.Div([
                dbc.Label("Differential expression cutoff (LHSAR vs LNCaP) [log2FoldChange]"),
                dcc.Slider(0, 20, 0.5,
                    value=0,
                    marks=None,
                    tooltip={"placement": "bottom", "always_visible": True},
                    id='DGE_lhsar_lncap_input'
                    )], 
                id = 'element-to-hide3',
                style= {'display': 'block'} ),
            
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
                html.Div([
                        html.Button("Download the network", id="network_tsv"),
                        dcc.Download(id="download_network_tsv"),
                    ],
                    id="download_div",
                    style= {'display': 'none'}
                ),

				html.Img(id='display-value', style={'height':'83%', 'width':'83%'})
			], color="primary", type="grow"),
	])
])


#####
@app.callback(
   [Output(component_id='element-to-hide1', component_property='style'),
   Output(component_id='element-to-hide2', component_property='style'),
   Output(component_id='element-to-hide3', component_property='style')],
   [Input(component_id='de_filter', component_property='value')])


def show_hide_element(visibility_state):
    if visibility_state == 'Gepia':
        return ({'display': 'block'}, {'display': 'none'}, {'display': 'none'})
    elif visibility_state == 'FPKM':
        return ({'display': 'none'}, {'display': 'block'}, {'display': 'none'})
    elif visibility_state == 'LHSAR VS LNCaP':
        return ({'display': 'none'}, {'display': 'none'}, {'display': 'block'})
    elif visibility_state == 'None':
        return ({'display': 'none'}, {'display': 'none'}, {'display': 'none'})


#####
@app.callback(
    Output("download_network_tsv", "data"),
    Input("network_tsv", "n_clicks"),
    State('promoter', 'value'),
    prevent_initial_call=True
)
def func(n_clicks,promoter):
    downloadfile = f"outfiles/{promoter}_network.tsv"
    if os.path.exists(downloadfile):
        return dcc.send_file(downloadfile)
    else:    
        return None

# @app.callback(
#     [Output("download_network_tsv", "data")],
#     [Input("network_tsv", "n_clicks"), 
#     State('promoter', 'value')],
    
#     prevent_initial_call=True
# )
# def func(network_tsv, promoter):
#     if promoter != 'Select a promoter Gene...':
#         print(f"outfiles/{promoter}_network.tsv")
#         return dcc.send_file(
#             "outfiles/network.tsv"
#         )
#     else:
#         return None


# def de_filter_observer(filter_state):
#     if filter_state == 'FPKM':        
#         vb.children = [ widgets.VBox(child_part1 + child_fpkm + child_part2, layout=box_layout),  outt]
#     elif filter_state == 'Gepia':
#         vb.children =  [ widgets.VBox( child_part1 + child_gepia + child_part2, layout=box_layout),  outt] 
#     elif button['new'] == 'None':
#         vb.children =  [ widgets.VBox( child_part1 + child_part2, layout=box_layout),  outt]
#     elif button['new'] == 'LHSAR VS LNCaP':
#         vb.children =  [ widgets.VBox( child_part1 + child_lhsar_lncap + child_part2, layout=box_layout),  outt] 
# #    vb.layout = outt


#####
# @app.callback(	
#     [Output('download_div', component_property = 'style')],
# 	[Input("submit_analysis", "n_clicks")])
# def remove_download_button(submit_analysis):
#     return ({'display': 'none'})

@app.callback(	
    [Output('display-value', 'src'),
    Output('download_div', component_property = 'style'),
    Output ('network_tsv','children')],    
	[Input("submit_analysis", "n_clicks"), 
    State('promoter', 'value'),
	State('DBP_size', 'value'),
	State('DBP_color', 'value'),
	State('intermediate_size', 'value'),
	State('intermediate_color', 'value'),
    State('de_filter', 'value'),
	State('DGE_tumorVSnormal_logFC_GEPIA_input', 'value'),
	State('FPKM_input', 'value'),
    State('DGE_lhsar_lncap_input', 'value'),
	State('SNP_BindingSite_path_input', 'value'),
	State('SNP_GenLoc_path_input', 'value'),
	State('only_enriched_nodes_input', 'value'),
	State('only_enriched_edges_input', 'value')],
    prevent_initial_call=True
)
def display_value(
        n_clicks,
        promoter, 
        DBP_size, 
        DBP_color,
        intermediate_size,
        intermediate_color,

        de_filter,
        DGE_tumorVSnormal_logFC_GEPIA_input,
        FPKM_input,
        DGE_lhsar_lncap_input,

        SNP_BindingSite_path_input,
        SNP_GenLoc_path_input,
        only_enriched_nodes_input,
        only_enriched_edges_input):

    if promoter == 'Select a promoter Gene...':
        return (None, {'display': 'none'})
    else:
        p = utils.single_promoter_graph(
            promoter, 
            DBP_size, 
            DBP_color,
            intermediate_size, 
            intermediate_color,

            de_filter,
            DGE_tumorVSnormal_logFC_GEPIA_input,
            FPKM_input,
            DGE_lhsar_lncap_input,

        
            SNP_BindingSite_path_input,
            SNP_GenLoc_path_input,
            only_enriched_nodes_input,
            only_enriched_edges_input)
        
        new_label_button = f'Download {promoter} network'
        if p[1] == 1:
            return (p[0], {'display': 'block'}, new_label_button)
        else:
            return (p[0], {'display': 'none'}, new_label_button)





########################## App Layout ##########################
app.layout = dbc.Container(fluid=True, children=[
    html.Br(),
    html.H1(config.app_name, id="nav-pills"),
    html.H4('Promoter-ENancher-GUided Interaction Networks'),
    navbar,
    html.Br(),
    html.Div(
    [
        html.Button("Submit", id="submit_analysis"),
    ]),
    
    html.Br(),        
    
    # html.Div(
    # [
    #     html.Button("Download the network", id="network_tsv"),
    #     dcc.Download(id="download_network_tsv"),
    # ],
    #     id="download_div",
    #     style= {'display': 'none'}),
    body,
    html.Br(),
    html.Br(),
    html.Br()
])


########################## Run ##########################

if __name__ == "__main__":
    debug = True if config.ENV == "DEV" else False
    app.run_server(debug=debug, host=config.host, port=config.port)

#debug = True if config.ENV == "DEV" else False
#app.run_server(debug=debug, host=config.host, port=config.port)
