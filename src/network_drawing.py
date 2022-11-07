# %%
import networkx as nx
import pandas as pd
import matplotlib as mpl
import numpy as np
import pygraphviz
from sys import argv
import os
import re

# %%
script, infile, outfile = argv
font = "Times"
fontsize = 150
shade_color = "#c6a1d1"
make_label = False
# %%
A = pd.read_csv('data/ref_files/unordered_patterns_adjacency_matrix.csv', index_col=0, header=0)

# %%
G = nx.from_pandas_adjacency(A)

# %%
attrs = pd.read_csv('data/ref_files/unordered_patterns_base_attributes.tsv', sep='\t', header=0)
#pos_dict = {row['node']: (row['x'], row['y']*-1) for i, row in attrs.iterrows()}
layer = {row['node']: row['x'] for i, row in attrs.iterrows()}
#nx.set_node_attributes(G, pos_dict, "pos")
nx.set_node_attributes(G, layer, 'layer')

# %%
sx = 200
sy = 150

nx.set_node_attributes(G, nx.multipartite_layout(G, subset_key="layer", scale=3), "pos")
for n in G:
    G.nodes[n]['pos'] = "{},{}!".format(G.nodes[n]['pos'][0]*sx, G.nodes[n]['pos'][1]*sy)

# %%
def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def colorFunc(c1,c2,prop):
    if prop >= 0.1:
        return(colorFader(c1,c2,1/2))
    elif prop < 0.1 and prop >= 0.01:
        return(colorFader(c1,c2,1/4))
    elif prop < 0.01 and prop >= 0.001:
        return(colorFader(c1,c2,1/16))
    else:
        return(colorFader(c1,c2,0))

def colorFunc2(c1,c2,prop,cap=0.1,scale=2/3):
    if prop > cap:
        return(colorFader(c1,c2,1*scale))
    else:
        prop = prop/cap
        return(colorFader(c1,c2,prop*scale))

# %%
# add proportion data to graph
proportion_file = pd.read_csv(infile, sep='\t', header=0)

# %%
prop_labels = {}
fontscale = 1
for i, row in proportion_file.iterrows():
    if row["prop"] > 0 and row["prop"] < 0.0001:
        decim, expon = f"{row['prop']:.2e}".split("e")
        sci_val = decim + '&#x00B7;10<sup>' + str(int(expon.replace('0',''))+2) + "</sup>%"
        prop_labels[row["node"]] = (
            f'<TR><TD><FONT FACE="{font}" POINT-SIZE="{fontsize*0.65}">{sci_val}</FONT>'
            f'<FONT FACE="{font}" POINT-SIZE="{fontsize*0.85}"> </FONT>'
            f'</TD></TR>'
        )                    
    else:
        sci_val = f"{row['prop']:.2%}"
        fontscale = 0.85
        prop_labels[row["node"]] = (
            f'<TR><TD><FONT FACE="{font}" POINT-SIZE="{fontsize*0.85}">'
            f'{sci_val}</FONT></TD></TR>'
        )

# %%
node_labels = {row['node']: f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="2">'
                            f'<TR><TD><FONT FACE="{font}" POINT-SIZE="{fontsize}">'
                            f'<B>({row["node"].replace(",",", ")})</B>'
                            f'</FONT></TD></TR><TR><TD></TD></TR>'
                            f'{prop_labels[row["node"]]}</TABLE>>'
                            for i, row in proportion_file.iterrows()}

nx.set_node_attributes(G, node_labels, 'label')

# %%
prop_colors = {row['node']: colorFunc2("white", shade_color, row['prop'], scale=1) for i, row in proportion_file.iterrows()}
nx.set_node_attributes(G, prop_colors, 'fillcolor')

# %%
if not os.path.exists('figures'):
    os.mkdir('figures')

# %%
g_label = re.match(r'\d+_g-(\d+)_.*', os.path.basename(infile)).group(1)
plot_label = f"g = {g_label:4}"

# %%
A = nx.nx_agraph.to_agraph(G)  # convert to a graphviz graph
A.layout(prog="neato", args="-n")
A.node_attr["shape"] = "box"
A.node_attr["style"] = "rounded,filled,solid"
A.graph_attr["overlap"] = "scale"
A.node_attr["penwidth"] = "18"
A.node_attr["margin"] = "0"
A.edge_attr["penwidth"] = "20"
A.graph_attr["margin"] = "0.25,0.25"
#A.graph_attr["nodesep"] = "2"
A.graph_attr["size"] = "3.5,2.3"
A.graph_attr["ratio"] = "fill"
if make_label:
    A.graph_attr["label"] = plot_label
    A.graph_attr["labelloc"] = "bottom"
    A.graph_attr["labeljust"] = "right"
    A.graph_attr["fontsize"] = fontsize*2.5
A.graph_attr["fontname"] = font
A.graph_attr["dpi"] = "300"
A.draw(outfile, prog="neato", args="-n")  # Draw with pygraphviz
#A.write("figures/out.dot")

# %%
