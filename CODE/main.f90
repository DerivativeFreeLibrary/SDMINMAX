program minmax
	IMPLICIT NONE
	INTEGER N,m

	INTEGER NUM_FUNCT,NUM_ITER,IPRINT
	real*8 f,mu,alfamax
	character*40 nomefun

	REAL*8,allocatable :: X(:),Z(:),Z1(:),D(:),D1(:),XOLD(:)    
	REAL*8,allocatable :: DOLDALFA(:),DCONV(:),FINIT(:,:) 
	REAL*8,allocatable :: FSTOP(:),XFSTOP(:,:),fs(:)

	COMMON /NUM/F
	COMMON /NUMNEW/NUM_FUNCT
	common/calfamax/alfamax
	common/cfun/m

	write(*,*)
	write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'


	call set_dimension(n,m,nomefun)

	allocate ( X(N),Z(N),Z1(N),D(N),D1(N),XOLD(N) )    
	allocate ( DOLDALFA(N),DCONV(N),FINIT(N,2)    ) 
	allocate ( FSTOP(N+1),XFSTOP(N,N+1),fs(m)     )


	NUM_FUNCT=0

	IPRINT=0

	CALL STARTP(N,X)

	call functs(n,m,x,fs)

	CALL MAINBOX(N,X,D,D1,Z,Z1,XOLD,NUM_ITER,DOLDALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,mu)
	
	call functs(n,m,x,fs)
	call functmu(n,x,mu,f)

	write(*,*) ' fmu = ',f
	write(*,*) ' max = ',maxval(fs)
	write(*,*) '  mu = ',mu

	write(*,1000) nomefun, n, m, num_funct, maxval(fs), f, mu

	deallocate ( X,Z,Z1,D,D1,XOLD     )    
	deallocate ( DOLDALFA,DCONV,FINIT ) 
	deallocate ( FSTOP,XFSTOP,fs      )
	
	write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
	write(*,*)

1000 format(1x,A40,    ' & ', i3,' & ', i3,' & ',    i8, ' & ',   d9.3, ' & ',     d9.3,' & ',    d9.3,' \\ \hline'  )

END program minmax

!===================================================================
!===================================================================

subroutine functmu(n,x,mu,fmu)
	implicit none
	integer n,m,i
	real*8  x(n),mu,fmu
	real*8  f(m),sum,alfamax
	logical lflag,lflag1

	common/calfamax/alfamax
	common/cflag/lflag,lflag1
	common/cfun/m

	call functs(n,m,x,f)
	
	lflag  = .true.

	sum = 0.0d0
	do i = 1,m
		!write(*,*) 'fmu: fi/mu = ',f(i)/mu
		if(f(i)/mu > 500.0d0) then
			sum = exp(500.0d0)
			exit
		else
			sum = sum + dexp(f(i)/mu)
		endif
	enddo

	!write(*,*) 'fmu: sum = ',sum
	if(sum < 1.0d-300 ) then
		fmu = mu*log(1.0d-300)
	else
		fmu = mu*dlog(sum)
	endif

	return
end subroutine functmu

