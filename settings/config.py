
import os

ENV = "DEV"
#ENV = "PROD"

## fs
#root = os.path.dirname(os.path.dirname(__file__)) + "/"

## server
host = "0.0.0.0"
port = int(os.environ.get("PORT", 5000))


## info
app_name = "PENGUIN"
contacts = ""
code = "https://github.com/bsc-life/penguin_analytics"
tutorial = "https://raw.githubusercontent.com/bsc-life/penguin_analytics/main/src/documentation.html?token=GHSAT0AAAAAABXOBGV44IJWCVIDWGWJOO3UYXT3ICQ"
fontawesome = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"

about = "PENGUIN (Promoter-ENancher-GUided Interaction Networks) is a tool to prioritize protein-protein interactions in enhancer-promoter loops involved in prostate cancer by integrating several sources of information (high-coverage H3K27ac-HiChIP data, DNA-binding motifs, tissue-specific physical protein interactions, gene expression, and genome-wide genetic variants)."
