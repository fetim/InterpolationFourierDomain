PROGRAM FourierInterpolation

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
  Nt = 100
  Nt_c = Nt/2 + 1

  ! Alocating
  ALLOCATE(in(Nt))
  ALLOCATE(out(Nt))
  ALLOCATE(t(Nt))
  ALLOCATE(y(Nt))
  ALLOCATE(f(Nt_c))

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
  y = sin(2*pi*frequency*t)


  ! save function
  open(25,file='gaussfunction.dat',status='unknown',form='formatted')
  do k=1,Nt     
     !     y(k) = rand()
     write(25,*) t(k),y(k)
  end do
  close(25)

  !prepare fft
  CALL dfftw_plan_dft_r2c_1d_ (plan_forward,Nt,y,out,FFTW_ESTIMATE)
  
  !execute
  CALL dfftw_execute_ (plan_forward)
  
  print*, size(out)
  f0 = 0.0
  tf =2.0/(2*dt) 
  df = (tf - f0)/(Nt_c-1)
  
  f  = (/(f0 + (k-1)*df ,k=1,Nt_c,1)/)
  open(26,file='gaussfunction_fourierdomain.dat',status='unknown',form='formatted')

  ! do k=1,Nt_c
  !    write(26,*) f(k),2.0/Nt*abs(out(k))
  ! end do


  do k=1,Nt
     write(26,*) k,abs(out(k))
  end do

  close(26)

  call dfftw_destroy_plan_ ( plan_forward )

END PROGRAM FourierInterpolation
