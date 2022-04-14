from dash import Dash, dcc, html, Input, Output
import os
import utils

external_stylesheets = [
    {
        "href": "https://fonts.googleapis.com/css2?family=Lato:wght@400;700&display=swap",
        "rel": "stylesheet",
    },
]

app = Dash(__name__, external_stylesheets=external_stylesheets)

app.title = "PENGUIN"

server = app.server

dpath = 'graph_tables'
prots=sorted([f[:-10] for f in os.listdir(dpath) if f.endswith('_nodes.tsv')])

app.layout = html.Div([
    dcc.Dropdown(prots,
        'FOXA1',
        id='dropdown'
    ),
    html.Img(id='display-value', style={'height':'10%', 'width':'10%'})
])

@app.callback(
	Output('display-value', 'src'),
	[Input('dropdown', 'value')]
)

def display_value(prot):
    p = utils.single_promoter_graph({prot})
    return(p)

if __name__ == '__main__':
    app.run_server(debug=True)

