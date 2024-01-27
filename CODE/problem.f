C-----------------------------------------
C	Hald Madsen 2 test problem
C	n = 5		m = 42
C	f(x*) = 0.000122
C-----------------------------------------
	subroutine set_dimension(n,m,nomefun)
		implicit none
		integer n,m
		character*40 nomefun
		real*8	eps

		common/cvinceps/eps

		eps = 1.0d-1

		n = 5
		m = 42
		nomefun = 'hald-madsen 2'

		return
	end subroutine set_dimension

        subroutine STARTP(N,X)
	integer n
	real*8 x(n)
	x(1) =  0.5d0
	x(2) =  0.0d0
	x(3) =  0.0d0
	x(4) =  0.0d0
	x(5) =  0.0d0

	return
	end subroutine startp

	subroutine functs(n,m,x,f)
	implicit none
	integer n,m,i
	real*8 x(n),f(m),y(21),num,den

	do i = 1,21
		y(i) = -1.0d0 + 0.1d0*( dble(i) - 1.0d0)
		num  = ( x(1) + x(2)*y(i) )
		den  = ( 1.0d0 + x(3)*y(i) + x(4)*y(i)**2 + x(5)*y(i)**3)
		f(i) = num / den - exp(y(i))
	enddo

	do i = 22,42
		f(i) = -f(i-21)
	enddo

	return
	end subroutine functs
