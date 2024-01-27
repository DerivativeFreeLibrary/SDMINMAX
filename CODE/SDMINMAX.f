C============================================================================================
C    SDMINMAX - A Derivative-free Method for Nonlinear Finite Minimax Optimization 
C    Copyright (C) 2011  G.Liuzzi, S.Lucidi, M.Sciandrone
C
C    This program is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program.  If not, see <http://www.gnu.org/licenses/>.
C
C    G. Liuzzi, S. Lucidi, M. Sciandrone. A Derivative-free Algorithm for Linearly
C    Constrained Finite Minimax Problems, SIAM J. on Optimization, 16(4): 1054-1075 (2006)
C
C============================================================================================
      SUBROUTINE MAINBOX(N,X,D,D1,Z,Z1,XOLD,NUM_ITER,
     *DOLDALFA,IPRINT,DCONV,FSTOP,XFSTOP,FINIT,mu)

      IMPLICIT NONE
      INTEGER N,I,J,I_CORR,NUM_FUNCT,NUM_ITER
      INTEGER NUM_FAL,NUM_SUCC,ISTOP
      INTEGER IPRINT,I_CORR_FALL
	INTEGER IMIN,IMAX,IMINALFA,IMAXALFA

      REAL*8 DCONV(N),DNR,mu,maxalfa
      REAL*8 X(N),Z(N),Z1(N),D(N),D1(N),XOLD(N)
      REAL*8 DOLDALFA(N),GAMMA,RHO,ALFA,ALFAEX,ALFAOLD
      REAL*8 F,FZ,FZDELTA,FZOLD,ftemp 
      REAL*8 DOLDALFAMEDIO,DALFAMAX
      REAL*8 FMIN,FMAX,ALFA0,DOLDALFAMIN,DOLDALFAMAX,SUM
      REAL*8 RAPALFA,FINIT(N,2)

      REAL*8 FSTOP(N+1),XFSTOP(N,N+1),WEIGHT(N+1)
      logical lflag, lflag1

      COMMON /NUM/F
      COMMON /NUMNEW/NUM_FUNCT
      common/calfamax/dalfamax
      common/cflag/lflag,lflag1

      NUM_SUCC=0

	mu = 1.0d1

      lflag = .false.
      lflag1 = .false.

      GAMMA=1.D-6
      RHO=1.D-1

      NUM_FUNCT = 0
      NUM_ITER = 0 
      NUM_FAL=0
      ISTOP = 0

      I_CORR=1

      ALFA0=1.D0

      DO I=1,N

        DOLDALFA(I)=1.D0

        DOLDALFA(I)=1.0D0*DMAX1(1.D-3,DMIN1(1.D+0,DABS(X(I))))
      
        IF(IPRINT.GE.2) THEN
          WRITE(*,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
          WRITE(1,*) ' ALFAiniz(',I,')=',DOLDALFA(I)
        ENDIF

      END DO
	dalfamax = maxval(doldalfa)

        IF(IPRINT.GE.0) THEN
          WRITE(*,*) ' ALFAmax =',Dalfamax
          WRITE(1,*) ' ALFAmax =',Dalfamax
        ENDIF
	
      DO I=1,N      
        D(I)=1.D0 
      END DO
       
      CALL FUNCTMU(N,X,mu,F)
	NUM_FUNCT=NUM_FUNCT+1

      FSTOP(I_CORR)=F

      DO I=1,N
        XFSTOP(I,I_CORR)=X(I)
	  Z(I)=X(I)
      END DO

      IF(IPRINT.GE.2) THEN
        WRITE(*,*) ' ----------------------------------'
        WRITE(1,*) ' ----------------------------------'
        WRITE(*,*) ' Finiz =',F
        WRITE(1,*) ' Finiz =',F
        DO I=1,N
          WRITE(*,*) ' Xiniz(',I,')=',X(I)
          WRITE(1,*) ' Xiniz(',I,')=',X(I)
        ENDDO
      ENDIF
 
      NUM_FAL=0

  1   CONTINUE

	maxalfa = maxval(doldalfa)

      IF(I_CORR.EQ.1) THEN
           DO I=1,N
                DCONV(I)=-D(I)
           END DO
      ENDIF

      IF(IPRINT.GE.1) THEN
        WRITE(*,*) '----------------------------------------------'
        WRITE(1,*) '----------------------------------------------'
        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
      ENDIF

      IF (ISTOP.EQ.1) THEN
   50   FORMAT(1X,a3,i5,a3,i5,a3,d13.5,a6)
        WRITE(*,*) '----------------------------------------------'
        WRITE(1,*) '----------------------------------------------'
        WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
        return
      END IF

      CALL LINESEARCH(N,X,F,D,ALFA,ALFAOLD,DOLDALFA,Z,FZ,FZOLD,GAMMA,
     *RHO,I_CORR,NUM_FAL,NUM_ITER,DALFAMAX,I_CORR_FALL,IPRINT,mu)
      
	maxalfa = max(maxalfa,alfa)
	  
      IF(DABS(ALFA).GE.1.D-24) THEN

	    FSTOP(I_CORR)=FZOLD

	    DO J=1,N
	       IF(J.EQ.I_CORR) THEN
	          XFSTOP(I_CORR,I_CORR)= X(I_CORR)+ALFAOLD*D(I_CORR)
             ELSE
                XFSTOP(J,I_CORR)=X(J)
	       ENDIF
          END DO 
		           
          
          X(I_CORR) = X(I_CORR)+ALFA*D(I_CORR)
          Z(I_CORR) = X(I_CORR)

          F=FZ
      
      
          NUM_FAL=0
          NUM_ITER=NUM_ITER+1
          NUM_SUCC=NUM_SUCC+1
    
          IF(IPRINT.GE.1) THEN
             WRITE(*,*) ' F =',F
             WRITE(1,*) ' F =',F
          ENDIF

          IF(IPRINT.GE.2) THEN
	       DO I=1,N
                WRITE(*,*) ' X(',I,')=',X(I)
                WRITE(1,*) ' X(',I,')=',X(I)
             ENDDO
          ENDIF
      
	ELSE
      
	    IF(I_CORR_FALL.EQ.0) THEN 

		   FSTOP(I_CORR)=FZ
	       DO J=1,N
                XFSTOP(J,I_CORR)=Z(J)
             END DO            

             NUM_FAL=NUM_FAL+1
             NUM_ITER=NUM_ITER+1
             Z(I_CORR) = X(I_CORR)
	    ENDIF
		      
	END IF

      IF(I_CORR.LT.N) THEN

          I_CORR=I_CORR+1

      ELSE

	    IF (ISTOP.EQ.1) THEN
             WRITE(*,*) '----------------------------------------------'
             WRITE(1,*) '----------------------------------------------'
             WRITE(*,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
             WRITE(1,*) 'NF=',NUM_FUNCT,'   F=',F,'   ALFAmax=',DALFAMAX
             return
          END IF


          FMIN=F
          FMAX=F
          IMIN=1
          IMAX=1
          DOLDALFAMIN=DOLDALFA(1)
          DOLDALFAMAX=DOLDALFA(1)
          IMINALFA=1
          IMAXALFA=1
          do i=2,n
             if(doldalfa(i).lt.doldalfamin) then
                doldalfamin=doldalfa(i)
                iminalfa=i
             end if 
             if(doldalfa(i).gt.doldalfamax) then
                doldalfamax=doldalfa(i)
                imaxalfa=i
             end if 
          end do
   
          rapalfa=1.d+6

          if(doldalfamax/doldalfamin.gt.rapalfa) then
             do i=1,n
                 D1(i)=dconv(i)
             end do
             dnr=dsqrt(DFLOAT(N))
	    else
             do i=1,n
                if(fstop(i).lt.fmin) then
                   fmin=fstop(i)
                   imin=i
                end if
                if(fstop(i).gt.fmax) then
                   fmax=fstop(i)
                   imax=i
                end if
             end do


             DNR=0.D0
             DOLDALFAMEDIO=(DOLDALFAMAX+DOLDALFAMIN)/2.D0
             DO I=1,N
	          IF(IMIN.GT.0) THEN
                   D1(I)=(XFSTOP(I,IMIN)-XFSTOP(I,IMAX))/DOLDALFAMAX
                ELSE
                   D1(I)=(X(I)-XFSTOP(I,IMAX))/DOLDALFAMAX
	          ENDIF
             DNR=DNR+D1(I)*D1(I)
             END DO
             DNR=DSQRT(DNR)


             DNR=DSQRT(DNR)

             IF(DNR.LE.1.D-24) THEN
	           do i=1,n
                   D1(i)=dconv(i)
                 end do
                 dnr=dsqrt(DFLOAT(N))
             ENDIF
          

          endif

          CALL LINESEARCH1(n,x,f,z1,d1,dnr,alfa0,gamma,alfa,fz,
     *                   iprint,NUM_FUNCT,mu)

		maxalfa = max(maxalfa,alfa)

          if(dabs(alfa).ge.1.d-24) then
             do i=1,n
                x(i) = x(i)+alfa*d1(i)
	          Z(I)=X(I)
             end do
             f=fz
 
             num_fal=0


             num_iter=num_iter+1
    
             if(iprint.ge.1) then
                write(*,*) ' f =',f
                write(1,*) ' f =',f
             endif

             if(iprint.ge.2) then
                do i=1,n
                   write(*,*) ' x(',i,')=',x(i)
                   write(1,*) ' x(',i,')=',x(i)
                enddo
             endif

          else

             alfa0=0.5d0*alfa0

          end if

		CALL STOP(N,DOLDALFA,ISTOP,DALFAMAX,NUM_FUNCT,FSTOP,F)

		if(mu > 1.0d+0*sqrt(maxalfa)) then
			mu = 1.0d+0*sqrt(maxalfa)
			CALL FUNCTMU(N,X,mu,F)
			NUM_FUNCT=NUM_FUNCT+1
		endif

          I_CORR=1

      END IF 

	if(((.NOT.lflag).AND.(dalfamax <  1.0d-3)).or.
     *   (lflag .AND.(dalfamax >= 1.0d-3))) then
		CALL FUNCTMU(N,X,mu,F)
		NUM_FUNCT=NUM_FUNCT+1
	endif

	lflag1 = lflag
	CALL FUNCTMU(N,X,mu,Ftemp)
	lflag = lflag1

      GO TO 1

      END
        

      SUBROUTINE STOP(N,DOLDALFA,ISTOP,DALFAMAX,NUM_FUNCT,FSTOP,F)

      IMPLICIT NONE
      
      INTEGER N,ISTOP,I,NUM_FUNCT
      REAL*8 DOLDALFA(N),DALFAMAX,FSTOP(N+1),FFSTOP,FFM,F


      ISTOP=0

      DALFAMAX=DOLDALFA(1)
      DO I=1,N
        IF(DOLDALFA(I).GT.DALFAMAX) THEN
          DALFAMAX=DOLDALFA(I)
        END IF
      END DO
     
      FFM=F
      DO I=1,N
        FFM=FFM+FSTOP(I)
      ENDDO
      FFM=FFM/DFLOAT((N+1))

      FFSTOP=(F-FFM)*(F-FFM)
      DO I=1,N
        FFSTOP=FFSTOP+(FSTOP(I)-FFM)*(FSTOP(I)-FFM)
      ENDDO
      FFSTOP=DSQRT(FFSTOP/DFLOAT(N+1))
      IF(DALFAMAX.LE.1.D-4) THEN
        ISTOP = 1
      END IF
      IF(FFSTOP.LE.1.D-4) THEN
      END IF
      IF(NUM_FUNCT.GT.500000) THEN
        ISTOP = 1
      END IF

      RETURN

      END

      SUBROUTINE LINESEARCH(N,X,F,D,ALFA,ALFAOLD,DOLDALFA,Z,FZ,
     *FZOLD,GAMMA,RHO,I_CORR,NUM_FAL,NUM_ITER,
     *DALFAMAX,I_CORR_FALL,IPRINT,mu)
     
      IMPLICIT NONE

      INTEGER N,I_CORR,NUM_FUNCT
      INTEGER I,J,L,LL
      INTEGER NUM_ITER,NUM_FAL
      INTEGER IPRINT,I_CORR_FALL
      REAL*8 X(N),D(N),DOLDALFA(N),Z(N),Z1(N)
      REAL*8 F,ALFA,FZ,GAMMA,RHO,DNR,mu
      REAL*8 DELTA,DELTA1,FPAR,FZDELTA
      REAL*8 ALFAEX,FMIN 
      REAL*8 DALFAMAX,FCOMMON
      REAL*8 ALFAOLD,FZOLD

      COMMON /NUM/FCOMMON
	COMMON /NUMNEW/NUM_FUNCT

      DELTA =0.5D0
      DELTA1 =0.5D0
      DNR=1.D0
      I_CORR_FALL=0

      J=I_CORR

      IF(IPRINT.GE.1) THEN
         WRITE(*,*) ' J=',J,'  D(J)=',D(J),'  DOLDALFA=',DOLDALFA(J)
         WRITE(1,*) ' J=',J,'  D(J)=',D(J),'  DOLDALFA=',DOLDALFA(J)
      ENDIF

 10   CONTINUE

      ALFA=DOLDALFA(J)

      IF(DABS(ALFA).LE.1.D-3*DALFAMAX) THEN

          I_CORR_FALL=1
          D(J)=-D(J)
	   
	   ALFA=0.D0
         RETURN
      
      END IF

      Z(J) = X(J)+ALFA*D(J)

      ALFAEX=ALFA
       
      CALL FUNCTMU(N,Z,mu,FZ)
      NUM_FUNCT=NUM_FUNCT+1

      IF(IPRINT.GE.1) THEN
         WRITE(*,*) ' FZ =',FZ,'   ALFA =',ALFA
         WRITE(1,*) ' FZ =',FZ,'   ALFA =',ALFA
      ENDIF
      IF(IPRINT.GE.2) THEN
         DO I=1,N
            WRITE(*,*) ' Z(',I,')=',Z(I)
            WRITE(1,*) ' Z(',I,')=',Z(I)
         ENDDO
      ENDIF

        
      FPAR= F-GAMMA*ALFA*ALFA*DNR*DNR

      IF(FZ.LE.FPAR) THEN
         FMIN=FZ 
         FZOLD=F
	   ALFAOLD=0.D0

   11    CONTINUE

         ALFAEX=ALFA/DELTA1

         Z(J) = X(J)+ALFAEX*D(J)
               
         CALL FUNCTMU(N,Z,mu,FZDELTA)
         NUM_FUNCT=NUM_FUNCT+1

         IF(IPRINT.GE.1) THEN
            WRITE(*,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX  
            WRITE(1,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX
         ENDIF
         IF(IPRINT.GE.2) THEN
             DO I=1,N
                WRITE(*,*) ' Z(',I,')=',Z(I)
                WRITE(1,*) ' Z(',I,')=',Z(I)
             ENDDO
         ENDIF

         FPAR= F-GAMMA*ALFAEX*ALFAEX*DNR*DNR
         IF(FZDELTA.LE.FPAR) THEN
                FZOLD=FZ
	          ALFAOLD=ALFA
                FZ=FZDELTA
                ALFA=ALFAEX
                GO TO 11
         ELSE
                  
              DOLDALFA(J)=DELTA*ALFA

              RETURN
         END IF

      ELSE 

          DOLDALFA(J)=DELTA*DOLDALFA(J)

          D(J)=-D(J)

          IF(IPRINT.GE.1) THEN
              WRITE(*,*) ' direzione opposta'
              WRITE(1,*) ' direzione opposta'
          ENDIF

          ALFA=0.D0

          RETURN

      END IF

      END

      SUBROUTINE LINESEARCH1(N,X,F,Z,D,DNR,ALFA0,GAMMA,
     *                       ALFA,FZ,IPRINT,NUM_FUNCT,mu)
      IMPLICIT NONE
      INTEGER N,IMIN,IMAX,I,J,NUM_FUNCT,IPRINT
      REAL*8 X(N),Z(N),F,ALFA0,FZ,GAMMA,DNR,mu
      REAL*8 FPAR,ALFA,ALFAEX,DELTA,DELTA1,FZDELTA,D(N),FMIN



      DELTA =0.5D0
      DELTA1 =0.5D0

      ALFA=ALFA0

      DO I=1,N
        Z(I)=X(I)+ALFA*D(I)
      END DO

      CALL FUNCTMU(N,Z,mu,FZ)
      NUM_FUNCT=NUM_FUNCT+1

      IF(IPRINT.GE.1) THEN
	   WRITE(*,*) ' DNR =',DNR
         WRITE(1,*) ' DNR =',DNR
         WRITE(*,*) ' FZ =',FZ,'   ALFA =',ALFA,'   ALFA0 =',ALFA0
         WRITE(1,*) ' FZ =',FZ,'   ALFA =',ALFA,'   ALFA0 =',ALFA0
      ENDIF
      IF(IPRINT.GE.2) THEN
         DO I=1,N
            WRITE(*,*) ' Z(',I,')=',Z(I)
            WRITE(1,*) ' Z(',I,')=',Z(I)
         ENDDO
       ENDIF

        FPAR= F-GAMMA*ALFA*ALFA*DNR*DNR
 11     CONTINUE
        IF(FZ.LE.FPAR) THEN
              ALFAEX=ALFA/DELTA1
              DO I=1,N
                Z(I) = X(I)+ALFAEX*D(I)
              END DO      
              CALL FUNCTMU(N,Z,mu,FZDELTA)
              NUM_FUNCT=NUM_FUNCT+1
              FPAR= F-GAMMA*ALFAEX*ALFAEX*DNR*DNR
              IF(IPRINT.GE.1) THEN
                WRITE(*,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX  
                WRITE(1,*) ' FZex=',FZDELTA,'  ALFAEX=',ALFAEX
              ENDIF
              IF(IPRINT.GE.2) THEN
                DO I=1,N
                  WRITE(*,*) ' Z(',I,')=',Z(I)
                  WRITE(1,*) ' Z(',I,')=',Z(I)
                ENDDO
              ENDIF
               IF(FZDELTA.LE.FPAR) THEN
                   FZ=FZDELTA
                   ALFA=ALFAEX
                   GO TO 11
                ELSE
                   ALFA0=ALFA
                   RETURN
                END IF
        ELSE
         ALFA=0.D0
        END IF

              
       END

