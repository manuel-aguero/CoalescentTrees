__author__ = "Manuel Aguero"

import numpy as np

#This function generates a random sequence
def randseq(l):
    seq = ""
    for i in range(l):
        seq += np.random.choice(["A","C","G","T"],p=[0.25,0.25,0.25,0.25])
    return seq

#This function mutates given sequence according to the Jukes-Cantor model of mutation over a given length.
#mutations from a base to itself are allowed
def mutate(x,t,mu):
    l = len(x)
    numMutation = np.random.poisson(l*mu*t) #the number of mutations from poisson distribution with total rate l*mu*t
    for i in range(numMutation):
        site = np.random.randint(0,l-1) #pick a site to mutate
        x = x[:site]+np.random.choice(["A","C","G","T"],p=[0.25,0.25,0.25,0.25])+x[site+1:]
    return x

#This function simulates a sequence down a given tree according to Jukes-Cantor substitution model by iterating 
#through each nodes and setting the new generated sequence
def jc_sequence_simulation(tree,l,mu):
    tree.get_root().set_sequence(randseq(l))
    jc_tree_traversal(tree.get_root(),l,mu)

#This function mutates and sets sequence to the child nodes of a given input node
def jc_tree_traversal(node,l,mu):
    if not node.is_leaf():
        children = node.get_children()
        for child in children:
            branch_length = node.get_height() - child.get_height()
            mutated_sequence = mutate(node.get_sequence(),branch_length,mu)
            child.set_sequence(mutated_sequence)
            jc_tree_traversal(child,l,mu)
