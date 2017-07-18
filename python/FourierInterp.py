#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
FourierInterp:

Interpolation in Fourier Domain.
The number of samples in Frequency domain
is proportional to the sampling rate in time 
domain.

@author: felipe
@email : felipetimoteo@id.uff.br
"""

def main():
    #%%
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.fftpack 
     
    incr_sampling = 1
    L  = 2.0
    dt = 0.001
    Nt = L/dt + 1
    
     
    t = np.linspace(0.0,dt*Nt,Nt)
    #y = np.exp(-100*(t-0.3)*(t-0.3))
    y = np.sin(2*np.pi*50*t)
     
    plt.figure(1)
    plt.plot(t,y,color='k',linestyle='dashed',marker='o')
    plt.show()
     
    
    Y = scipy.fftpack.rfft(y)
     
    size_Y = np.size(Y)    
    X = np.linspace(0.0,1.0/(2.0*dt),Nt/2)
     
    
    Y_ex = np.insert(Y,size_Y/2,np.zeros(incr_sampling*Nt))
     
    plt.figure(2)
     
    plt.plot(X,2.0/Nt*np.abs(Y[:Nt/2]))
    plt.show()
     
    y_interp = (incr_sampling+1)*scipy.fftpack.irfft(Y_ex)
    t_interp = np.linspace(0.0,dt*Nt,(incr_sampling+1)*Nt)
     
    plt.figure(1)
    plt.plot(t,y,color='k',linestyle='dashed',marker='o')
    plt.plot(t_interp,y_interp,color='r',marker='o')
    plt.show()

    #%%
if __name__=="__main__":
    main()
