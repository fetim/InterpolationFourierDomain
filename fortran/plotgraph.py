#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
plotgraph.py:

Plot a Graph from a ascii file with two columns.
1st column is x direction
2nd comlun is y direction

@author: felipe
@email : felipetimoteo@id.uff.br
"""

def main():
    
    plot_from_file('gaussfunction.dat')

    plot_from_file('gaussfunction_fourierdomain.dat')

    plot_from_file('gaussfunction_interpolated.dat')

def plot_from_file(filename):
    import numpy as np
    import matplotlib.pyplot as plt
    
    #filename="gaussfunction.dat"
    
    data = np.loadtxt(filename,float)
    
    x = data[:,0]
    y = data[:,1]


    plt.plot(x,y,marker='o')
    plt.show()
    
if __name__=="__main__":
    main()
