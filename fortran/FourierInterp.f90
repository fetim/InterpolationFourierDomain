PROGRAM FourierInterpolation
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ### Compiling ### - $ gfortran FourierInterp.f90 -Ldirectory -lfftw3 -o run_Fourier
    !  ### Execute ###  - $ ./run_Fourier
    !
    ! FourierInterp.f90 - Fourier Interpolation /  Zero Padding in the Fourier Domain.
    !   This program interpolates a function y(t).
    !   The function y(t) is defined inside the program.
    !   
    ! INPUT:  
    !       t0 = initial t  
    !       tf = final t
    !       Nt = number of sample
    !       factor = interpolation factor
    !       t0 = initial t
    !       frequency1 = 10 Hz
    !       frequency2 = 20 Hz
    !       frequency3 =  2 Hz
    !       y(t) = 0.1*sin(2*pi*frequency1*t)  + sin(2*pi*frequency2*t) + sin(2*pi*frequency3*t)
    !
    ! OUTPUT: 
    !       function.dat - t,y(t)
    !       function_fourierdomain.dat - f,Y(f)
    !       function_interpolated.dat - t,y_interp(t)
    ! 
    ! Code Written by Felipe Timoteo
    !                 email : felipetimoteo@id.uff.br
    !                 Last update: July 22th, 2017
    !
    ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
    !                    Departamento de Geologia e Geofísica
    !                    Universidade Federal Fluminense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE

  INCLUDE "fftw3.f90"

  INTEGER    ,PARAMETER                :: sp=4, dp=8, qp=16
  INTEGER(sp),PARAMETER                :: seed = 8198272  
  INTEGER                              :: k,Nt,factor
  REAL(dp)  ,PARAMETER                 :: pi = 4.0*atan(1.0)
  REAL(dp)                             :: dt,dt_interp,t0,tf
  REAL(dp)                             :: df,f0,ff,Fs
  REAL(dp)                             :: frequency
  REAL(dp)                             :: plan_backward,plan_forward
  

  REAL(dp)   ,DIMENSION(:),ALLOCATABLE    :: t,t_interp,y,y_interp,f
  COMPLEX(dp),DIMENSION(:),ALLOCATABLE :: out,out_expanded

  
  ! Parameters  
  t0 = 0.0
  tf = 1.0
  Nt = 100
  factor = 100


  ! Alocating
  ALLOCATE(out(Nt))
  ALLOCATE(out_expanded(factor*Nt))
  ALLOCATE(t(Nt))
  ALLOCATE(t_interp(factor*Nt))
  ALLOCATE(y(Nt))
  ALLOCATE(y_interp(factor*Nt))
  ALLOCATE(f(Nt))
  
  out_expanded = 0.0

  !Sampling
  dt = (tf - t0)/(Nt-1)
  !Frequency Sampling
  Fs = 1/dt
  print*,''
  print*,'Sample rate = ',dt, 'Frequency sampling = ', Fs
  print*,'Nyquist Frequency = ', Fs/4
  print*,''

  ! vector
  t = (/(t0 + (k-1)*dt,k=1,Nt,1)/)
  
  ! calculation
  !  y = exp(-1000*(t-0.3)*(t-0.3))
  frequency = 10 !Hz
  y = 0.1*sin(2*pi*frequency*t)  + sin(2*pi*2*t) + sin(2*pi*20*t)


  ! save function
  open(25,file='function.dat',status='unknown',form='formatted')
  do k=1,Nt     
     !     y(k) = rand()
     write(25,*) t(k),y(k)
  end do
  close(25)

  !prepare fft
  CALL dfftw_plan_dft_r2c_1d_ (plan_forward,Nt,y,out,FFTW_ESTIMATE)
  
  !execute fft
  CALL dfftw_execute_ (plan_forward)


  ! estimating frequencies
  f0 = 0.0
  df = Fs/size(y)
  f  = (/(f0 + (k-1)*df ,k=1,Nt)/)
 
  ! save fourier transform of function
  open(26,file='function_fourierdomain.dat',status='unknown',form='formatted')
  do k=1,Nt
     out_expanded(k) = out(k)
     write(26,*) f(k),abs(out(k)/Nt)
  end do
  close(26)

  ! prepare inverse fft
  CALL dfftw_plan_dft_c2r_1d_(plan_backward,factor*Nt,out_expanded,y_interp,FFTW_ESTIMATE)

  ! execute inverse fft
  CALL dfftw_execute(plan_backward)

  !New Sampling
  dt_interp = dt/factor
!  dt_interp = (tf - t0)/(2*Nt-1)
  !New vector
  t_interp = (/(t0 + (k-1)*dt_interp,k=1,factor*Nt,1)/)

  ! save interpolated function
  open(27,file='function_interpolated.dat',status='unknown',form='formatted')
  do k=1,factor*Nt     
     write(27,*) t_interp(k),y_interp(k)/Nt
  end do
  close(25)

  !end fft
  CALL dfftw_destroy_plan_ ( plan_forward )
  CALL dfftw_destroy_plan_ ( plan_backward )



END PROGRAM FourierInterpolation
