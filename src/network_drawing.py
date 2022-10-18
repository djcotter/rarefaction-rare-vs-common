import networkx as nx
import matplotlib.pyplot as plt
from grave import plot_network

graph = nx.powerlaw_cluster_graph(50, 1, .2)


