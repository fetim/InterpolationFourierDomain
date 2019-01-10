# InterpolationFourierDomain
This project implements Fourier Domain Interpolation in Fortran and Python

Prerequisite:
-------
* FFTW library
* python

### Installing FFTW
[FFTW]( http://www.fftw.org/) - downloading from the FFTW home page.
unpacking with
```bash
$ tar zxvf fftw-x.x.x.tar.gz
```
Enter in the folder and execute: FFTW3 (default: gfortran and gcc)

```bash
$ ./configure --enable-threads --enable-openmp --enable-avx
$ ./configure 
```
（For gromacs ）
```bash
$ ./configure --enable-threads --enable-float
```
and then:
```bash
$ make
$ sudo make install
```