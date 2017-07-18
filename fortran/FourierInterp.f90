PROGRAM FourierInterpolation
!Compile - $ gfortran FourierInterp.f90 -Ldirectory -lfftw3 -o run_Fourier
  IMPLICIT NONE

  INCLUDE "fftw3.f90"

  INTEGER    ,PARAMETER                :: sp=4, dp=8, qp=16
  INTEGER(sp),PARAMETER                :: seed = 8198272  
  INTEGER                              :: k,Nt,Nt_c
  REAL(dp)  ,PARAMETER                 :: pi = 4.0*atan(1.0)
  REAL(dp)                             :: dt,t0,tf
  REAL(dp)                             :: df,f0,ff,Fs
  REAL(dp)                             :: frequency
  REAL(dp)                             :: plan_backward,plan_forward
  
  REAL(dp)   ,DIMENSION(:),ALLOCATABLE    :: in
  REAL(dp)   ,DIMENSION(:),ALLOCATABLE    :: t,y,f
  COMPLEX(dp),DIMENSION(:),ALLOCATABLE :: out 

  ! seed random number mechanism
  call srand(seed)
  
  ! Parameters
  t0 = 0.0
  tf = 1.0
  Nt = 1001
  Nt_c = Nt/2 + 1

  ! Alocating
  ALLOCATE(in(Nt))
  ALLOCATE(out(Nt))
  ALLOCATE(t(Nt))
  ALLOCATE(y(Nt))
  ALLOCATE(f(Nt))

  !Sampling
  dt = (tf - t0)/(Nt-1)
  !Frequency Sampling
  Fs = 1/dt
  print*,''
  print*,'Sample rate = ',dt, 'Frequency sampling = ', Fs
  print*,'Nyquist Frequency = ', Fs/2
  print*,''

  ! vector
  t = (/(t0 + (k-1)*dt,k=1,Nt,1)/)
  
  ! calculation
  !  y = exp(-1000*(t-0.3)*(t-0.3))
  frequency = 10 !Hz
  y = sin(2*pi*frequency*t)  +sin(2*pi*2*t) + sin(2*pi*20*t)


  ! save function
  open(25,file='gaussfunction.dat',status='unknown',form='formatted')
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
  open(26,file='gaussfunction_fourierdomain.dat',status='unknown',form='formatted')
  do k=1,Nt
     print*,'k,f(k)=>',k,f(k),abs(out(k)/Nt)
     write(26,*) f(k),abs(out(k)/Nt)
  end do
  close(26)

  !end fft
  call dfftw_destroy_plan_ ( plan_forward )

END PROGRAM FourierInterpolation
