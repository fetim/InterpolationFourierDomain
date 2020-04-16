def readbinaryfile(dim1,dim2,filename):
    """
    readbinaryfile - Functions that read a binary file.
    Usage
    Input:
    dim1     = Number of sample of 1st Dimension
    dim2     = Number of sample of 2nd Dimension
    filename = path of binary file     
    """      
    with open(filename, 'rb') as f:    
        data   = np.fromfile(filename, dtype= np.float32, count= dim1*dim2)
        matrix = np.reshape(data, [dim1,dim2], order='F')
    return matrix
      
# -*- coding: utf-8 -*-
"""
FourierTransform

This program test Fourier Transform in a simple example.
"""

import numpy as np
import matplotlib.pyplot as pl
import scipy.fftpack

# parameters
dt = 0.001
Nt = 1501
Nx = 250
time = np.arange(0,dt*Nt,dt)
f_nqust = 1/(2*dt)

# filename="/home/truta/Desktop/ViscoelasticCongressPaper/python/sismo_tpp_e.bin"
# matrix = readbinaryfile(Nt,Nx,filename)
# matrix = np.transpose(matrix)
# signal = matrix[:,100]

signal = (np.sin(2*np.pi*time*(20/(dt*Nt))) + np.sin(2*np.pi*time*(50/(dt*Nt))))

pl.plot(time,signal)
pl.show(block=False)

# Fourier parameters
freq_ini=0.1
freq_fin=40
df= 0.1
frequency = np.arange(freq_ini,freq_fin+df,df)

Cos_signal = np.zeros(len(frequency))
Sin_signal = np.zeros(len(frequency))

# Fourier transform
for m in range(len(frequency)):
    print(m,frequency[m])    
    for n in range(0,Nt):        
        cos = signal[n]*np.cos(2*np.pi*frequency[m]*n*dt)
        Cos_signal[m] = Cos_signal[m] + cos
        
        sin = signal[n]*np.sin(2*np.pi*frequency[m]*n*dt)
        Sin_signal[m] = Sin_signal[m] + sin

# spectrum plot
Pot_signal = np.sqrt(Cos_signal*Cos_signal + Sin_signal*Sin_signal)

pl.figure()
pl.plot(frequency,Pot_signal)
xmin, xmax, ymin, ymax = pl.axis()
xmin=0
xmax=40
ymin=0
# ymax=14.785355758666991
pl.axis([xmin, xmax, ymin, ymax])
pl.show(block=False)

# compare with scipy
spectrum = np.abs(scipy.fftpack.fft(signal))
freq = scipy.fftpack.fftfreq(Nt,dt)
pl.figure()
pl.plot(freq,spectrum,color='black',label='sample')
xmin, xmax, ymin, ymax = pl.axis()
xmin=0
xmax=40
ymin=0
# ymax=14.785355758666991
pl.axis([xmin, xmax, ymin, ymax])

pl.show(block=False)

# gradiente frequencia
# 1 - criar conjunto de frequencia
# frq:0.1:frq_max
# 2*pi*(frq_ini:frq_fin)
# 2 - transformada de fourier numerida
# criar seno e coseno
# sinF = sin(Freq(f)*n*dt)*20.0f*dt
# usen =  u(x) * sen(w*t*dt)
# ucos =  u(x) * cos(w*t*dt)