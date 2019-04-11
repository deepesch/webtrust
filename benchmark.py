################################################################################
# CS 224W (Fall 2017) - Project
# Weighted Sign Networks
# Author: akshayk@stanford.edu, vhying@stanford.edu
# Last Updated: Dec 06, 2017
################################################################################


import graph_tool.all as gt
import networkx as nx
import tidaltrust as tt


# Data Load
print "Loading Weighted Graphs:"
print "Loading bitcoin alpha graph..."
fh = open("Datasets/soc-sign-bitcoinalpha.csv", "r")
btcAlphaGr = nx.read_weighted_edgelist(fh, delimiter=',', nodetype=int)
fh.close()


# Tidaltrust of single edge.
u = 1
v = 888
if (btcAlphaGr.has_edge(u, v)):
  print tt.compute_trust(btcAlphaGr, 1, 888)["trust"]


# Eigentrust
# Trust values must lie between [0,1]

g = gt.Graph()
# Add vertices
for node in btcAlphaGr.nodes():
  g.add_vertex(node)
# Add edges
for u, v in btcAlphaGr.edges():
  g.add_edge(u, v)
# Add edge weight???????

#gt.centrality.eigentrust(g, trust_map, vprop=None, norm=False, epsilon=1e-06, max_iter=0, ret_iter=False)
