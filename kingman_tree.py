__author__ = "Manuel Aguero"

import random
import numpy as np
#This function simulates a tree according to Kingsman coalescent model with
#n = number of leaves and Ne = effective population size and returns a the simulated tree
def kingman_tree(n,ne):
    k = n
    t = 0
    nodes = list()
    #Make n leaf nodes with time t = 0 and labeled from 1 to n. This is the set of available nodes.
    for i in range(n):
        nodes.append(Node(str(i+1)))
        nodes[i].set_height(t)
    while k > 1:
        kc2 = k*(k-1)/2
        tk = np.random.exponential(ne/kc2) #A time sample is generated using the total rate of coalescence.
        t = t + tk
        indices = random.sample(range(len(nodes)),2) #randomly choosing two discting nodes index
        indices.sort()    
        label = nodes[indices[0]].get_label()+nodes[indices[1]].get_label()
        nodeM = Node(label) #create new parent node
        nodeM.set_height(t)
        nodeM.add_child(nodes[indices[0]])  #assign sampled nodes as childs for the new node
        nodeM.add_child(nodes[indices[1]])
        nodes[indices[0]] = nodeM #remove samples from available nodes
        del nodes[indices[1]]
        k = k-1 #Adjust the size of the set of available nodes
    return Tree(nodeM)
