#!/usr/bin/python
import networkx as nx
from matplotlib import pyplot
from scipy import io
from pylab import *
import sys, getopt


if len(sys.argv) == 1:
 print "sys.argv[1] file.mtx"
 exit(2)

"... iniciando grafo" 
A = io.mmread(sys.argv[1])
G = nx.from_scipy_sparse_matrix(A)
"......................................................................"

"..."
adjacency_matrix = nx.to_numpy_matrix(G)

fig = pyplot.figure(figsize=(10,10))
pyplot.imshow(adjacency_matrix,cmap="Greys",interpolation="none")
show()
