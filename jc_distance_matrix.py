__author__ = "Manuel Aguero"
import numpy as np

#This function calculates the Jukes-Cantor distance between two sequences
def jc_distance(s1,s2):
    count = sum(s1[i]!= s2[i] for i in range(len(s1))) #count the number of different sites between s1 and s2
    if count == 0:
        return 0  
    p = min((count/len(s1)),(0.75-1/len(s1))) #p-distance between two sequences so that the distances are well defined
    d = -(3/4)*np.log(1-(4/3)*p) #Jukes-Cantor distance between s1 and s2
    return d

#This function creates the Jukes-Cantor distance matrix
def jc_distance_matrix(sequences):
    distance_matrix = np.zeros(shape=[len(sequences), len(sequences)])
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            distance_matrix[i][j] = (jc_distance(sequences[i],sequences[j]))
    return(distance_matrix)
