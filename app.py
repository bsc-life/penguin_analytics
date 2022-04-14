from dash import Dash, dcc, html, Input, Output
import os
import utils
import ipywidgets as widgets

external_stylesheets = [
    {
        "href": "https://fonts.googleapis.com/css2?family=Lato:wght@400;700&display=swap",
        "rel": "stylesheet",
    },
]

app = Dash(__name__, external_stylesheets=external_stylesheets)

app.title = "Avocado Analytics: Understand Your Avocados!"

server = app.server

dpath = 'graph_tables'
node_stats = [
    'degree', 
    'expression_log2FoldChange', 'expression_padj', 
    'SNP-binding',
    'Betweenness_global', 'Degree_global', 'Betweenness_in_cluster',
    'Betweenness_out_cluster', 'Betweenness_enrichment',
    'Degree_in_cluster', 'Degree_out_cluster', 'Degree_enrichment',
    'Betweenness_enrichment_pval', 'Degree_enrichment_pval', 'fisher_test_OR', 'fisher_test_pval'
]
edge_stats = [
    None,
    'fisher_test_edge_OR', 
    'fisher_test_edge_pvalue',
    'SNP-path', 
    'druggable'
]

prots=sorted([f[:-10] for f in os.listdir(dpath) if f.endswith('_nodes.tsv')])
app.layout = html.Div([
    html.H2('Hello World'),
    dcc.Dropdown(prots,
        'FOXA1',
        id='dropdown'
    ),
    html.Div(id='display-value')
])

@app.callback(Output('display-value', 'children'),
                [Input('dropdown', 'value')])
#def display_value(value):
#    return f'You have selected {value}'
def display_value(prot):
    p = widgets.interact(
    utils.single_promoter_graph({prot})
    )
    return(p)

if __name__ == '__main__':
    app.run_server(debug=True)

