import numpy as np
from scipy import spatial
import networkx as nx
import matplotlib.pyplot as plt
from treelib import Node, Tree
import sys

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

def plot_tree(infilename,outfilename):
    Np=1
    fin = open(infilename, "r")
    positions = dict()
    G=list()
    line = fin.readline()
    numbers = line.split()
    if numbers[0]=='A)':
        Np=int(numbers[2])
    for i in range(Np):
        G.append(nx.Graph())
    for line in fin:
        numbers = line.split()
        if str(numbers[0])=="C)":
            positions[str(numbers[3])] = (float(numbers[4]), float(numbers[5]))
        if numbers[0]=="D)":
            G[int(numbers[1])].add_edge(str(numbers[3]), str(numbers[4]))    
    fin.close()
    n = len(colors)
    for i in range(Np):
        nx.draw(G[i], positions, node_color=colors[i%n], node_shape='.', node_size=100, with_labels=False)
    plt.axis('equal')
    plt.savefig(outfilename)
    plt.clf()
    


def print_tree(filename) :
    Np=1
    fin = open(filename, "r")
    positions = dict()
    tree=list()
    line = fin.readline()
    positions = dict()
    numbers = line.split()
    if numbers[0]=='A)':
        Np=int(numbers[2])
    for i in range(Np):
        tree.append(Tree())    
    for line in fin:
        numbers = line.split()
        if str(numbers[0])=="B)":
            Harry =  "("+str(numbers[4]) +", " +str(numbers[5]) + ")"
            harry = str(numbers[3])
            tree[int(numbers[1])].create_node(Harry, harry) 
        if str(numbers[0])=="C)":
            positions[ numbers[3]] = (float(numbers[4]), float(numbers[5])) 
        if numbers[0]=="D)": 
            Jane = str(positions[numbers[4]])
            jane = str(numbers[4])
            harry=str(numbers[3])
            tree[int(numbers[1])].create_node(Jane,jane, parent=harry)
    fin.close()         
    for i  in range(Np):        
        tree[i].show()

        
if  len(sys.argv) < 2 or  len(sys.argv) > 3  :
    print("Usage: python %s inputfile [outputfile]" % sys.argv[0] )
    exit(1)
infile = sys.argv[1]
if len(sys.argv) == 3:
    outfile = sys.argv[2]
else:
    outfile = "kd_tree.png"
print("input file: %s" % infile )
print("output file: %s" % outfile )
print("Tree:\n")
plot_tree(infile,outfile)
print_tree(infile)
