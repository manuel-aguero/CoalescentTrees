__author__  = "Manuel Aguero"
import numpy
from scipy.stats import lognorm

def internal_nodes_heights(node,nodes_heights):
    if not node.is_leaf():
        nodes_heights.append(node.get_height())
        for child in node.get_children():
            internal_nodes_heights(child,nodes_heights)

def mcmc(tree,max_steps,chain):
    ne = 1
    #generate list of times with k linages
    time_linage = list()
    internal_nodes_heights(tree.get_root(),time_linage)
    time_linage.append(0)
    time_linage.sort()
    for i in range(1,len(time_linage)):
        time_linage[i] = (time_linage[i] - time_linage[i-1])
    del time_linage[0]
    time_linage.reverse()
    for i in range(max_steps):
        u = np.random.uniform(-10,10) #random step
        ne2 = ne + u
        if(ne2 > 0): #handle by rejectioon
            likelihood1 = 0
            likelihood2 = 0
            #generate log of the prior probability densities
            prior1 = lognorm.logpdf(ne, s = 1.25,scale=np.exp(3)) 
            prior2 = lognorm.logpdf(ne2, s = 1.25,scale=np.exp(3)) 
            for k in range(2,len(tree.get_leaves())):  #time interval with k linages
                tk = time_linage[k-2]
                b = (k*(k-1)*tk)/(2*ne)
                c = (k*(k-1)*tk)/(2*ne2)
                #compute the coalescent likelihood
                likelihood1 += (-np.log(ne) - b) 
                likelihood2 += (-np.log(ne2) - c)
            #obtain posterior distributions
            post1 = likelihood1 + (prior1) 
            post2 = likelihood2 + (prior2)
            post = post2 - post1
            #compute acceptance probability
            if post > -800 and post < 700: #handle overflow by reject large values
                a = min(1,np.exp(post))
                v = np.random.uniform(0,1)
                if(a > v):
                    ne = ne2
        chain.append(ne) #add to MCMC chain
