MODULE VARIABLES
IMPLICIT NONE
INTEGER, PARAMETER :: N=121,M=101,L=481		!181 141 661 !121 101 421
INTEGER, PARAMETER :: Imax1=20,Imax2=40,      Imax3=40,Imax4=20        ! total= all this +1!Imax1=35,Imax2=55,      Imax3=55,Imax4=35
INTEGER, PARAMETER :: Jmax1=10,Jmax2=40,      Jmax3=40,Jmax4=10
INTEGER, PARAMETER :: Kmax1=10,Kmax2=50,      Kmax3=50,Kmax4=10
INTEGER :: TTT,OT,OTF,FT,RERUN,REIT,ITER
DOUBLE PRECISION :: T1,T2
DOUBLE PRECISION :: X(-1:N+1),Y(-1:M+1),Z(-1:L+1)
DOUBLE PRECISION :: U(-1:N+1,-1:M+1,-1:L+1) 
DOUBLE PRECISION :: V(-1:N+1,-1:M+1,-1:L+1)
DOUBLE PRECISION :: W(-1:N+1,-1:M+1,-1:L+1)
DOUBLE PRECISION :: UN1(0:N,0:M,0:L)
DOUBLE PRECISION :: VN1(0:N,0:M,0:L)
DOUBLE PRECISION :: WN1(0:N,0:M,0:L)
DOUBLE PRECISION :: PRES(-1:N+1,-1:M+1,-1:L+1)
DOUBLE PRECISION :: FUU(0:N,0:M,0:L),FVV(0:N,0:M,0:L),FWW(0:N,0:M,0:L)
DOUBLE PRECISION :: FUUN1(0:N,0:M,0:L),FVVN1(0:N,0:M,0:L),FWWN1(0:N,0:M,0:L)
DOUBLE PRECISION :: FUUN2(0:N,0:M,0:L),FVVN2(0:N,0:M,0:L),FWWN2(0:N,0:M,0:L)
DOUBLE PRECISION :: DP(-1:N+1,-1:M+1,-1:L+1)
DOUBLE PRECISION :: DX(-1:N),DX1(-1:N)
DOUBLE PRECISION :: DY(-1:M),DY1(-1:M)
DOUBLE PRECISION :: DZ(-1:L),DZ1(-1:L)
DOUBLE PRECISION :: A(-1:M+1),B(-1:M+1),C(-1:L+1)
DOUBLE PRECISION :: CD,CDP,CDV,CL,CLP,CLV
DOUBLE PRECISION :: LAMDA=0.78,PI=3.141592654
DOUBLE PRECISION :: DT,OMEGA,RE,KC,PERIOD,TIME,PHI
DOUBLE PRECISION :: LENGTH=50.0,HEIGHT=40.0,WIDTH=12.0
INTEGER :: N_THREADS,THREAD_NUM
INTEGER :: I_WORK_CHUNK,J_WORK_CHUNK,K_WORK_CHUNK
INTEGER :: I_LBOUND,I_UBOUND,J_LBOUND,J_UBOUND,K_LBOUND,K_UBOUND
END MODULE VARIABLES

MODULE IBM
USE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION :: Fx(0:N,0:M,0:L),Fy(0:N,0:M,0:L),Fz(0:N,0:M,0:L)                           
DOUBLE PRECISION :: USOLID=0.D0
DOUBLE PRECISION :: VSOLID=0.D0
DOUBLE PRECISION :: WSOLID=0.D0
DOUBLE PRECISION :: DISTANCE,XC=25.0,YC=20.0,ZC=6.0,R=0.5
INTEGER :: ETA(0:N,0:M,0:L)                                                
END MODULE IBM

MODULE COEFFICIENT
USE VARIABLES
IMPLICIT NONE 
DOUBLE PRECISION ::totalFx,totalFy,totalFz,CENTRAL_LW
DOUBLE PRECISION ::CDx,CDy,CDZ
DOUBLE PRECISION ::XL=0.D0    
DOUBLE PRECISION ::YL=0.D0
DOUBLE PRECISION ::ZL=0.D0                                                   
END MODULE COEFFICIENT
!
!       THIS PROGRAM IS FOR SOLVING 3D OSCILLATING FLOW BY 
!       PROJECTION METHOD
!       
!       AUTHOR: M. J. CHERN   DATE:17 NOVEMBER,2016
!              
!       STAGGERED GRID (MAC)
!
PROGRAM FLOW3D
USE VARIABLES
IMPLICIT NONE

OPEN(43,FILE='CLOCK.dat')
CALL cpu_time(T1)
!
!       INPUT DATA OF THE FUNCTION U(X,Y) ON COLLOCATION POINTS
!
CALL INPUT
PRINT *,'INPUT PASS'
!
CALL SOLVER 
CALL cpu_time(T2)
WRITE(43,*) 'Time gaken by code was',T2-T1,'seconds'
PRINT *,'SOLVER PASS'
!
STOP
END PROGRAM FLOW3D
!-------------------------------------------------------------
!
   SUBROUTINE INPUT
   USE VARIABLES
   USE IBM
   IMPLICIT NONE
   INTEGER :: I,J,K

	OPEN(UNIT=15,FILE='XX.dat')
        DO I=-1,N+1,1
          READ(15,*)X(I)
        END DO
	CLOSE(15)
!
	OPEN(UNIT=16,FILE='YY.dat')
        DO J=-1,M+1,1
           READ(16,*)Y(J)
        END DO
	CLOSE(16)
!
        OPEN(UNIT=17,FILE='ZZ.dat')
        DO K=-1,L+1,1
           READ(17,*)Z(K)
        END DO
        CLOSE(17)
!
        DO I=-1,N
           DX(I)=X(I+1)-X(I)
        END DO
!        
        DO I=-1,N-1
           DX1(I)=0.5*(X(I+2)-X(I))
        END DO
!        DX1(N)=0.5*(X(N)-X(N-1))
!        DX1(0)=0.5*(X(2)-X(1))
!        DX1(N-1)=0.5*(X(N)-X(N-1))
!        
        DO J=-1,M
           DY(J)=Y(J+1)-Y(J)
        END DO
!
        DO J=-1,M-1
           DY1(J)=0.5*(Y(J+2)-Y(J))
        END DO
!        DY1(N)=0.5*(Y(M)-Y(M-1))
!        DY1(0)=0.5*(Y(2)-Y(1))
!        DY1(M-1)=0.5*(Y(M)-Y(M-1))        
!
        DO K=-1,L
           DZ(K)=Z(K+1)-Z(K)
        END DO
!
        DO K=-1,L-1
           DZ1(K)=0.5*(Z(K+2)-Z(K))
        END DO
!        DZ1(L)=0.5*(Z(L)-Z(L-1))
!        DZ1(0)=0.5*(Z(2)-Z(1))
!        DZ1(L-1)=0.5*(Z(L)-Z(L-1))
!
	OPEN(UNIT=19,FILE='PARA.F')
        READ(19,*)TTT,OT,OTF,DT,FT,ITER,RE,KC,OMEGA,RERUN
	CLOSE(19)
        TIME=0.D0
     DO I=-1,N+1,1
        DO J=-1,M+1,1
           DO K=-1,L+1,1
              U(I,J,K) = 0.D0
              V(I,J,K) = 0.D0
              W(I,J,K) = 0.D0
              PRES(I,J,K) = 0.D0
              END DO
           END DO
        END DO
!
    DO I=-1,N+1,1
       DO J=-1,M+1,1  
          DO K=-1,L+1,1               
	     USOLID=0.D0
             VSOLID=0.D0
             WSOLID=0.D0
          END DO
       END DO
    END DO
CALL VOS
RETURN
END SUBROUTINE INPUT
!-------------------------------------------------------------
! NM means Neumann boundary condition
!-----------------------------------------------------------
SUBROUTINE BOUNDARY
USE VARIABLES
IMPLICIT NONE
INTEGER I,J,K
PERIOD=KC 
!
! Y PLANE
!
DO K=0,L
  DO I=1,N-1
  
    ! the first y-plzne ( Oscillatoryy flow: U = sin(2*pi*Time/period) )       
    U(I,0,K)=U(I,1,K)	!+(2.D0*1.D0*DSIN(2.D0*PI*TIME/PERIOD))
    V(I,0,K)=0.D0
    W(I,0,K)=W(I,1,K)

    U(I,-1,K)=U(I,0,K)
    V(I,-1,K)=V(I,0,K)
    W(I,-1,K)=W(I,0,K)

    ! the final y-plzne ( Oscillatoryy flow: U = sin(2*pi*Time/period) )      
    U(I,M,K)=U(I,M-1,K)	!+(2.D0*1.D0*DSIN(2.D0*PI*TIME/PERIOD))
    V(I,M,K)=0.D0
    W(I,M,K)=W(I,M-1,K)

    U(I,M+1,K)=U(I,M,K)
    V(I,M+1,K)=V(I,M,K)
    W(I,M+1,K)=W(I,M,K)
  END DO      
END DO
!
! Z PLANE
!
DO J=0,M
  DO I=1,N-1

    ! the first z-plzne (NM)   
    U(I,J,0)=U(I,J,1)		!+2.D0*DSIN(2.D0*PI*TIME/PERIOD)
    V(I,J,0)=V(I,J,1)
    W(I,J,0)=0.D0
    
    U(I,J,-1)=U(I,J,0)   
    V(I,J,-1)=V(I,J,0)
    W(I,J,-1)=W(I,J,0)

    ! the final z-plzne (NM)
    U(I,J,L)=U(I,J,L-1)		!+2.D0*DSIN(2.D0*PI*TIME/PERIOD)
    V(I,J,L)=V(I,J,L-1)
    W(I,J,L)=0.D0
    
    U(I,J,L+1)=U(I,J,L)
    V(I,J,L+1)=V(I,J,L)
    W(I,J,L+1)=W(I,J,L)
  END DO
END DO	
!
! X PLANE
!
DO K=-1,L+1
  DO J=-1,M+1
  
    ! the first X-plzne ( Oscillatoryy flow: U = sin(2*pi*Time/period) )    
    U(0,J,K)=1.D0*DSIN(2.D0*PI*TIME/PERIOD)
    V(0,J,K)=-V(1,J,K)
    W(0,J,K)=-W(1,J,K)

    U(-1,J,K)=U(0,J,K)
    V(-1,J,K)=V(0,J,K)
    W(-1,J,K)=W(0,J,K)

    ! the final X-plzne ( Oscillatoryy flow: U = sin(2*pi*Time/period) )     
    U(N,J,K)=1.D0*DSIN(2.D0*PI*TIME/PERIOD)
    V(N,J,K)=-V(N-1,J,K)
    W(N,J,K)=-W(N-1,J,K)  

    U(N+1,J,K)=U(N,J,K)
    V(N+1,J,K)=V(N,J,K)
    W(N+1,J,K)=W(N,J,K)
  END DO
END DO    
RETURN
END SUBROUTINE BOUNDARY
!
!-------------------------------------------------------------
!
!SUBROUTINE PRESSURE_BOUNDARY
!USE VARIABLES
!MPLICIT NONE
!INTEGER :: I,J,K

!DO J=0,M,1
!  DO K=0,L,1
!     PRES(0,J,K)=PRES(1,J,K)
!     PRES(N,J,K)=PRES(N-1,J,K)
!  END DO
!END DO

!DO I=1,N,1
!  DO K=0,L,1
!     PRES(I,0,K)=PRES(I,1,K)
!     PRES(I,M,K)=PRES(I,M-1,K)                                    
!  END DO
!END DO

!DO I=1,N,1
!  DO J=0,M,1
!     PRES(I,J,0)=PRES(I,J,1)
!     PRES(I,J,L)=PRES(I,J,L-1)
!  END DO
!END DO
!RETURN
!END SUBROUTINE PRESSURE_BOUNDARY        
!
!--------------------------------------------------------------
!    GETTING  THE COORDINATE OF SOLID CENTER 
!    AND DETERMINDING ETA
!---------------------------------------------------------------
SUBROUTINE VOS                                                                 
USE VARIABLES
USE IBM
IMPLICIT NONE

    INTEGER :: I,J,K

    OPEN(UNIT=81,FILE='PARARA.F')
    READ(81,*)XC,YC,ZC,R
    CLOSE(81)
      
    OPEN(UNIT=88,FILE='ETA.dat')
    WRITE(88,*) 'VARIABLES=X,Y,Z,E'      
    WRITE(88,*) 'ZONE I=',N,',J=',M,',K=',L

DO K=1,L
  DO J=1,M
    DO I=1,N
       IF(Y(J)>YC .AND. X(I)>XC)THEN 
         DISTANCE=SQRT((((X(I)+X(I+1))/2.D0)-XC)**2+(((Y(J)+Y(J+1))/2.D0)-YC)**2)                      
           IF(DISTANCE > R) THEN
             ETA(I,J,K)=0.D0
           ELSE
             ETA(I,J,K)=1.D0
           END IF
        WRITE(88,*) X(I),Y(J),Z(K),ETA(I,J,K)                                
              
      ELSE IF(Y(J)>YC .AND. X(I)<=XC)THEN
        DISTANCE=SQRT((((X(I)+X(I-1))/2.D0)-XC)**2+(((Y(J)+Y(J+1))/2.D0)-YC)**2)                      
          IF(DISTANCE > R) THEN
            ETA(I,J,K)=0.D0
          ELSE
            ETA(I,J,K)=1.D0
          END IF
       WRITE(88,*) X(I),Y(J),Z(K),ETA(I,J,K) 

      ELSE IF(Y(J)<=YC .AND. X(I)>XC)THEN
        DISTANCE=SQRT((((X(I)+X(I+1))/2.D0)-XC)**2+(((Y(J)+Y(J-1))/2.D0)-YC)**2)                      
          IF(DISTANCE > R) THEN
            ETA(I,J,K)=0.D0
          ELSE
            ETA(I,J,K)=1.D0
          END IF
       WRITE(88,*) X(I),Y(J),Z(K),ETA(I,J,K)         
     
      ELSE IF(Y(J)<=YC .AND. X(I)<=XC)THEN
       	DISTANCE=SQRT((((X(I)+X(I-1))/2.D0)-XC)**2+(((Y(J)+Y(J-1))/2.D0)-YC)**2)                      
          IF(DISTANCE > R) THEN
            ETA(I,J,K)=0.D0
          ELSE
            ETA(I,J,K)=1.D0
          END IF
       WRITE(88,*) X(I),Y(J),Z(K),ETA(I,J,K) 
      END IF        
    END DO
  END DO
END DO           
    CLOSE(88)	    
    RETURN	  
END SUBROUTINE VOS
!
!------------------------------------------------------------
!
!------------------------------------------------------------
!   
SUBROUTINE QUICK
USE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION :: FUX,FUY,FUZ,FUC,FVX,FVY,FVZ,FVC,FWX,FWY,FWZ,FWC
DOUBLE PRECISION :: VISX,VISY,VISZ
DOUBLE PRECISION :: UE,UW
DOUBLE PRECISION :: VN,VS
DOUBLE PRECISION :: WT,WB
DOUBLE PRECISION :: X1,X2,X3,X4,X5,X6,X7,X8
DOUBLE PRECISION :: Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
DOUBLE PRECISION :: Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8
DOUBLE PRECISION :: A1,A2,A3,A4,A5,A6
DOUBLE PRECISION :: A7,A8,A9,A10,A11,A12
DOUBLE PRECISION :: VDX,VDY,VDZ
INTEGER :: I,J,K
DO K=1,L-1,1
  DO J=1,M-1,1
    DO I=1,N-2,1
      
       UE=0.5*(DX(I)*U(I+1,J,K)/DX1(I)+DX(I+1)*U(I,J,K)/DX1(I))
       UW=0.5*(DX(I)*U(I-1,J,K)/DX1(I-1)+DX(I-1)*U(I,J,K)/DX1(I-1))
       VN=0.5*(V(I,J,K)+V(I+1,J,K))
       VS=0.5*(V(I,J-1,K)+V(I+1,J-1,K))
       WT=0.5*(W(I,J,K)+W(I+1,J,K))
       WB=0.5*(W(I,J,K-1)+W(I+1,J,K-1))              

       X1=DX(I)/2.
       X2=DX(I+1)/2.
       X3=DX(I+1)+DX(I+2)/2.
       X4=DX(I-1)+DX(I)/2.
       X5=DX(I)/2.
       X6=DX(I-1)/2.
       X7=DX(I-1)+DX(I-2)/2.
       X8=DX(I)+DX(I+1)/2.   

       A1=X2*X4/(X4-X1)/(X1+X2)
       A2=X1*X4/(X1+X2)/(X2+X4)
       A3=-X1*X2/(X4-X1)/(X2+X4)
       A4=-X1*X3/(X1+X2)/(X2-X3)
       A5=X2*X3/(X1+X3)/(X1+X2)
       A6=X1*X2/(X1+X3)/(X2-X3)
       A7=X5*X7/(X7-X6)/(X5+X6)
       A8=X6*X7/(X5+X7)/(X5+X6)
       A9=-X5*X6/(X7-X6)/(X5+X7)
       A10=X6*X8/(X5+X6)/(X8-X5)
       A11=X5*X8/(X5+X6)/(X6+X8)
       A12=-X5*X6/(X6+X8)/(X8-X5)

       FUX=(DMAX1(UE,0.D0)*(A1*U(I,J,K)+A2*U(I+1,J,K)+A3*U(I-1,J,K))&
	-DMAX1(-UE,0.D0)*(A4*U(I+1,J,K)+A5*U(I,J,K)+A6*U(I+2,J,K))&
        -DMAX1(UW,0.D0)*(A7*U(I-1,J,K)+A8*U(I,J,K)+A9*U(I-2,J,K))&
        +DMAX1(-UW,0.D0)*(A10*U(I,J,K)+A11*U(I-1,J,K)+A12*U(I+1,J,K)))/DX(I)

       Y1=DY(J)/2.
       Y2=DY(J)/2.
       Y3=DY(J+1)+DY(J)/2.
       Y4=DY(J-1)+DY(J)/2.
       Y5=DY(J-1)/2.
       Y6=DY(J-1)/2.
       Y7=DY(J-2)+DY(J-1)/2.
       Y8=DY(J)+DY(J-1)/2.

       A1=Y2*Y4/(Y4-Y1)/(Y1+Y2)
       A2=Y1*Y4/(Y1+Y2)/(Y2+Y4)
       A3=-Y1*Y2/(Y4-Y1)/(Y2+Y4)
       A4=-Y1*Y3/(Y1+Y2)/(Y2-Y3)
       A5=Y2*Y3/(Y1+Y3)/(Y1+Y2)
       A6=Y1*Y2/(Y1+Y3)/(Y2-Y3)
       A7=Y5*Y7/(Y7-Y6)/(Y5+Y6)
       A8=Y6*Y7/(Y5+Y7)/(Y5+Y6)
       A9=-Y5*Y6/(Y7-Y6)/(Y5+Y7)
       A10=Y6*Y8/(Y5+Y6)/(Y8-Y5)
       A11=Y5*Y8/(Y5+Y6)/(Y6+Y8)
       A12=-Y5*Y6/(Y6+Y8)/(Y8-Y5)

       FUY=(DMAX1(VN,0.D0)*(A1*U(I,J,K)+A2*U(I,J+1,K)+A3*U(I,J-1,K))&
	-DMAX1(-VN,0.D0)*(A4*U(I,J+1,K)+A5*U(I,J,K)+A6*U(I,J+2,K))&
        -DMAX1(VS,0.D0)*(A7*U(I,J-1,K)+A8*U(I,J,K)+A9*U(I,J-2,K))&
	+DMAX1(-VS,0.D0)*(A10*U(I,J,K)+A11*U(I,J-1,K)+A12*U(I,J+1,K)))/DY1(J-1)

       Z1=DZ(K)/2.
       Z2=DZ(K)/2.
       Z3=DZ(K+1)+DZ(K)/2.
       Z4=DZ(K-1)+DZ(K)/2.
       Z5=DZ(K-1)/2.
       Z6=DZ(K-1)/2.
       Z7=DZ(K-2)+DZ(K-1)/2.
       Z8=DZ(K)+DZ(K-1)/2.

       A1=Z2*Z4/(Z4-Z1)/(Z1+Z2)
       A2=Z1*Z4/(Z1+Z2)/(Z2+Z4)
       A3=-Z1*Z2/(Z4-Z1)/(Z2+Z4)
       A4=-Z1*Z3/(Z1+Z2)/(Z2-Z3)
       A5=Z2*Z3/(Z1+Z3)/(Z1+Z2)
       A6=Z1*Z2/(Z1+Z3)/(Z2-Z3)
       A7=Z5*Z7/(Z7-Z6)/(Z5+Z6)
       A8=Z6*Z7/(Z5+Z7)/(Z5+Z6)
       A9=-Z5*Z6/(Z7-Z6)/(Z5+Z7)
       A10=Z6*Z8/(Z5+Z6)/(Z8-Z5)
       A11=Z5*Z8/(Z5+Z6)/(Z6+Z8)
       A12=-Z5*Z6/(Z6+Z8)/(Z8-Z5)
 
       FUZ=(DMAX1(WT,0.D0)*(A1*U(I,J,K)+A2*U(I,J,K+1)+A3*U(I,J,K-1))&
	-DMAX1(-WT,0.D0)*(A4*U(I,J,K+1)+A5*U(I,J,K)+A6*U(I,J,K+2))&
	-DMAX1(WB,0.D0)*(A7*U(I,J,K-1)+A8*U(I,J,K)+A9*U(I,J,K-2))&
	+DMAX1(-WB,0.D0)*(A10*U(I,J,K)+A11*U(I,J,K-1)+A12*U(I,J,K+1)))/DZ(K)

       FUC=0.D0 
       VDX=(U(I+1,J,K)-U(I,J,K))/DX(I)-(U(I,J,K)-U(I-1,J,K))/DX(I-1)
       VDX=VDX/DX1(I-1)
       VDY=(U(I,J+1,K)-U(I,J,K))/DY1(J)-(U(I,J,K)-U(I,J-1,K))/DY1(J-1)
       VDY=VDY/DY(J)
       VDZ=(U(I,J,K+1)-U(I,J,K))/DZ1(K)-(U(I,J,K)-U(I,J,K-1))/DZ1(K-1)
       VDZ=VDZ/DZ(K)
       VISX=(VDX+VDY+VDZ)/RE
       FUU(I,J,K)=DT*(-FUX-FUY-FUZ-FUC+VISX)

!       FUU(I,J,K)=(PRES(I,J,K)-PRES(I+1,J,K))/DX1(I)-FUX-FUY-FUZ-FUC+VISX
!             VISX=((U(I+1,J,K)-2.*U(I,J,K)+U(I-1,J,K))/DX1(I)/DX1(I)
!     1            +(U(I,J+1,K)-2.*U(I,J,K)+U(I,J-1,K))/DY(J)/DY(J)
!     2            +(U(I,J,K+1)-2.*U(I,J,K)+U(I,J,K-1))/DZ(K)/DZ(K))/RE

!                  print *,y1,y2,y3,y4,y5,y6,y7,y8
!                  print *,z1,z2,z3,z4,z5,z6,z7,z8
!                 print *,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
!                 print *,i,j,k,fux,fuy,fuz
!                 pause
    END DO
  END DO
END DO

DO K=1,L-1,1
  DO J=1,M-2,1
    DO I=1,N-1,1  
    
       VN=0.5*(DY(J+1)*V(I,J,K)/DY1(J)+DY(J)*V(I,J+1,K)/DY1(J))
       VS=0.5*(DY(J)*V(I,J-1,K)/DY1(J-1)+DY(J-1)*V(I,J,K)/DY1(J-1))
       UE=0.5*(U(I,J,K)+U(I,J+1,K))
       UW=0.5*(U(I-1,J,K)+U(I-1,J+1,K))
       WT=0.5*(W(I,J,K)+W(I,J+1,K))
       WB=0.5*(W(I,J,K-1)+W(I,J+1,K-1))

       X1=DX(I)/2.
       X2=DX(I)/2.
       X3=DX(I+1)+DX(I)/2.
       X4=DX(I-1)+DX(I)/2.
       X5=DX(I-1)/2.
       X6=DX(I-1)/2.
       X7=DX(I-2)+DX(I-1)/2.
       X8=DX(I)+DX(I-1)/2.

       A1=X2*X4/(X4-X1)/(X1+X2)
       A2=X1*X4/(X1+X2)/(X2+X4)
       A3=-X1*X2/(X4-X1)/(X2+X4)
       A4=-X1*X3/(X1+X2)/(X2-X3)
       A5=X2*X3/(X1+X3)/(X1+X2)
       A6=X1*X2/(X1+X3)/(X2-X3)
       A7=X5*X7/(X7-X6)/(X5+X6)
       A8=X6*X7/(X5+X7)/(X5+X6)
       A9=-X5*X6/(X7-X6)/(X5+X7)                 
       A10=X6*X8/(X5+X6)/(X8-X5)
       A11=X5*X8/(X5+X6)/(X6+X8)             
       A12=-X5*X6/(X6+X8)/(X8-X5)

       FVX=(DMAX1(UE,0.D0)*(A1*V(I,J,K)+A2*V(I+1,J,K)+A3*V(I-1,J,K))&
	-DMAX1(-UE,0.D0)*(A4*V(I+1,J,K)+A5*V(I,J,K)+A6*V(I+2,J,K))&
	-DMAX1(UW,0.D0)*(A7*V(I-1,J,K)+A8*V(I,J,K)+A9*V(I-2,J,K))&
	+DMAX1(-UW,0.D0)*(A10*V(I,J,K)+A11*V(I-1,J,K)+A12*V(I+1,J,K)))/DX1(I-1)

       Y1=DY(J)/2.
       Y2=DY(J+1)/2.
       Y3=DY(J+1)+DY(J+2)/2.
       Y4=DY(J-1)+DY(J)/2.
       Y5=DY(J)/2.
       Y6=DY(J-1)/2.
       Y7=DY(J-1)+DY(J-2)/2.
       Y8=DY(J)+DY(J+1)/2.

       A1=Y2*Y4/(Y4-Y1)/(Y1+Y2)
       A2=Y1*Y4/(Y1+Y2)/(Y2+Y4)
       A3=-Y1*Y2/(Y4-Y1)/(Y2+Y4)
       A4=-Y1*Y3/(Y1+Y2)/(Y2-Y3)
       A5=Y2*Y3/(Y1+Y3)/(Y1+Y2)
       A6=Y1*Y2/(Y1+Y3)/(Y2-Y3)
       A7=Y5*Y7/(Y7-Y6)/(Y5+Y6)
       A8=Y6*Y7/(Y5+Y7)/(Y5+Y6)
       A9=-Y5*Y6/(Y7-Y6)/(Y5+Y7)
       A10=Y6*Y8/(Y5+Y6)/(Y8-Y5)
       A11=Y5*Y8/(Y5+Y6)/(Y6+Y8)
       A12=-Y5*Y6/(Y6+Y8)/(Y8-Y5)

       FVY=(DMAX1(VN,0.D0)*(A1*V(I,J,K)+A2*V(I,J+1,K)+A3*V(I,J-1,K))&
	-DMAX1(-VN,0.D0)*(A4*V(I,J+1,K)+A5*V(I,J,K)+A6*V(I,J+2,K))&
	-DMAX1(VS,0.D0)*(A7*V(I,J-1,K)+A8*V(I,J,K)+A9*V(I,J-2,K))&
	+DMAX1(-VS,0.D0)*(A10*V(I,J,K)+A11*V(I,J-1,K)+A12*V(I,J+1,K)))/DY(J)

       Z1=DZ(K)/2.
       Z2=DZ(K)/2.
       Z3=DZ(K+1)+DZ(K)/2.
       Z4=DZ(K-1)+DZ(K)/2.
       Z5=DZ(K-1)/2.
       Z6=DZ(K-1)/2.
       Z7=DZ(K-2)+DZ(K-1)/2.
       Z8=DZ(K)+DZ(K-1)/2.

       A1=Z2*Z4/(Z4-Z1)/(Z1+Z2)
       A2=Z1*Z4/(Z1+Z2)/(Z2+Z4)
       A3=-Z1*Z2/(Z4-Z1)/(Z2+Z4)
       A4=-Z1*Z3/(Z1+Z2)/(Z2-Z3)
       A5=Z2*Z3/(Z1+Z3)/(Z1+Z2)
       A6=Z1*Z2/(Z1+Z3)/(Z2-Z3)
       A7=Z5*Z7/(Z7-Z6)/(Z5+Z6)
       A8=Z6*Z7/(Z5+Z7)/(Z5+Z6)
       A9=-Z5*Z6/(Z7-Z6)/(Z5+Z7)                 
       A10=Z6*Z8/(Z5+Z6)/(Z8-Z5)
       A11=Z5*Z8/(Z5+Z6)/(Z6+Z8)             
       A12=-Z5*Z6/(Z6+Z8)/(Z8-Z5)

       FVZ=(DMAX1(WT,0.D0)*(A1*V(I,J,K)+A2*V(I,J,K+1)+A3*V(I,J,K-1))&
	-DMAX1(-WT,0.D0)*(A4*V(I,J,K+1)+A5*V(I,J,K)+A6*V(I,J,K+2))&
	-DMAX1(WB,0.D0)*(A7*V(I,J,K-1)+A8*V(I,J,K)+A9*V(I,J,K-2))&
	+DMAX1(-WB,0.D0)*(A10*V(I,J,K)+A11*V(I,J,K-1)+A12*V(I,J,K+1)))/DZ(K)

       FVC=0.D0 
       VDX=(V(I+1,J,K)-V(I,J,K))/DX1(I)-(V(I,J,K)-V(I-1,J,K))/DX1(I-1)
       VDX=VDX/DX(I)
       VDY=(V(I,J+1,K)-V(I,J,K))/DY(J)-(V(I,J,K)-V(I,J-1,K))/DY(J-1)
       VDY=VDY/DY1(J-1)
       VDZ=(V(I,J,K+1)-V(I,J,K))/DZ1(K)-(V(I,J,K)-V(I,J,K-1))/DZ1(K-1)
       VDZ=VDZ/DZ(K)
       VISY=(VDX+VDY+VDZ)/RE
       FVV(I,J,K)=DT*(-FVX-FVY-FVZ-FVC+VISY)

!       FVV(I,J,K)=(PRES(I,J,K)-PRES(I,J+1,K))/DY1(J)-FVX-FVY-FVZ-FVC+VISY
!             VISY=((V(I+1,J,K)-2.*V(I,J,K)+V(I-1,J,K))/DX(I)/DX(I)
!     1            +(V(I,J+1,K)-2.*V(I,J,K)+V(I,J-1,K))/DY1(J)/DY1(J)
!     2            +(V(I,J,K+1)-2.*V(I,J,K)+V(I,J,K-1))/DZ(K)/DZ(K))/RE
    END DO
  END DO
END DO

DO K=1,L-2,1
  DO J=1,M-1,1
    DO I=1,N-1,1  
    
       WT=0.5*(DZ(K+1)*W(I,J,K)/DZ1(K)+DZ(K)*W(I,J,K+1)/DZ1(K))
       WB=0.5*(DZ(K)*W(I,J,K-1)/DZ1(K-1)+DZ(K-1)*W(I,J,K)/DZ1(K-1))
       UE=0.5*(U(I,J,K)+U(I,J,K+1))
       UW=0.5*(U(I-1,J,K)+U(I-1,J,K+1))
       VN=0.5*(V(I,J,K)+V(I,J,K+1))
       VS=0.5*(V(I,J-1,K)+V(I,J-1,K+1))

       X1=DX(I)/2.
       X2=DX(I)/2.
       X3=DX(I+1)+DX(I)/2.
       X4=DX(I-1)+DX(I)/2.
       X5=DX(I-1)/2.
       X6=DX(I-1)/2.
       X7=DX(I-2)+DX(I-1)/2.
       X8=DX(I)+DX(I-1)/2.   

       A1=X2*X4/(X4-X1)/(X1+X2)
       A2=X1*X4/(X1+X2)/(X2+X4)
       A3=-X1*X2/(X4-X1)/(X2+X4)
       A4=-X1*X3/(X1+X2)/(X2-X3)
       A5=X2*X3/(X1+X3)/(X1+X2)
       A6=X1*X2/(X1+X3)/(X2-X3)
       A7=X5*X7/(X7-X6)/(X5+X6)
       A8=X6*X7/(X5+X7)/(X5+X6)
       A9=-X5*X6/(X7-X6)/(X5+X7)
       A10=X6*X8/(X5+X6)/(X8-X5)
       A11=X5*X8/(X5+X6)/(X6+X8)
       A12=-X5*X6/(X6+X8)/(X8-X5)

       FWX=(DMAX1(UE,0.D0)*(A1*W(I,J,K)+A2*W(I+1,J,K)+A3*W(I-1,J,K))&
	-DMAX1(-UE,0.D0)*(A4*W(I+1,J,K)+A5*W(I,J,K)+A6*W(I+2,J,K))&
	-DMAX1(UW,0.D0)*(A7*W(I-1,J,K)+A8*W(I,J,K)+A9*W(I-2,J,K))&
	+DMAX1(-UW,0.D0)*(A10*W(I,J,K)+A11*W(I-1,J,K)+A12*W(I+1,J,K)))/DX(I)

       Y1=DY(J)/2.
       Y2=DY(J)/2.
       Y3=DY(J+1)+DY(J)/2.
       Y4=DY(J-1)+DY(J)/2.
       Y5=DY(J-1)/2.
       Y6=DY(J-1)/2.
       Y7=DY(J-2)+DY(J-1)/2.
       Y8=DY(J)+DY(J-1)/2.

       A1=Y2*Y4/(Y4-Y1)/(Y1+Y2)
       A2=Y1*Y4/(Y1+Y2)/(Y2+Y4)
       A3=-Y1*Y2/(Y4-Y1)/(Y2+Y4)
       A4=-Y1*Y3/(Y1+Y2)/(Y2-Y3)
       A5=Y2*Y3/(Y1+Y3)/(Y1+Y2)
       A6=Y1*Y2/(Y1+Y3)/(Y2-Y3)
       A7=Y5*Y7/(Y7-Y6)/(Y5+Y6)
       A8=Y6*Y7/(Y5+Y7)/(Y5+Y6)
       A9=-Y5*Y6/(Y7-Y6)/(Y5+Y7)
       A10=Y6*Y8/(Y5+Y6)/(Y8-Y5)
       A11=Y5*Y8/(Y5+Y6)/(Y6+Y8)
       A12=-Y5*Y6/(Y6+Y8)/(Y8-Y5)

       FWY=(DMAX1(VN,0.D0)*(A1*W(I,J,K)+A2*W(I,J+1,K)+A3*W(I,J-1,K))&
	-DMAX1(-VN,0.D0)*(A4*W(I,J+1,K)+A5*W(I,J,K)+A6*W(I,J+2,K))&
	-DMAX1(VS,0.D0)*(A7*W(I,J-1,K)+A8*W(I,J,K)+A9*W(I,J-2,K))&
	+DMAX1(-VS,0.D0)*(A10*W(I,J,K)+A11*W(I,J-1,K)+A12*W(I,J+1,K)))/DY(J)

       Z1=DZ(K)/2.
       Z2=DZ(K+1)/2.
       Z3=DZ(K+1)+DZ(K+2)/2.
       Z4=DZ(K)+DZ(K-1)/2.
       Z5=DZ(K)/2.
       Z6=DZ(K-1)/2.
       Z7=DZ(K-1)+DZ(K-2)/2.
       Z8=DZ(K)+DZ(K+1)/2.

       A1=Z2*Z4/(Z4-Z1)/(Z1+Z2)
       A2=Z1*Z4/(Z1+Z2)/(Z2+Z4)
       A3=-Z1*Z2/(Z4-Z1)/(Z2+Z4)
       A4=-Z1*Z3/(Z1+Z2)/(Z2-Z3)
       A5=Z2*Z3/(Z1+Z3)/(Z1+Z2)
       A6=Z1*Z2/(Z1+Z3)/(Z2-Z3)
       A7=Z5*Z7/(Z7-Z6)/(Z5+Z6)
       A8=Z6*Z7/(Z5+Z7)/(Z5+Z6)
       A9=-Z5*Z6/(Z7-Z6)/(Z5+Z7)
       A10=Z6*Z8/(Z5+Z6)/(Z8-Z5)
       A11=Z5*Z8/(Z5+Z6)/(Z6+Z8)
       A12=-Z5*Z6/(Z6+Z8)/(Z8-Z5)

       FWZ=(DMAX1(WT,0.D0)*(A1*W(I,J,K)+A2*W(I,J,K+1)+A3*W(I,J,K-1))&
	-DMAX1(-WT,0.D0)*(A4*W(I,J,K+1)+A5*W(I,J,K)+A6*W(I,J,K+2))&
	-DMAX1(WB,0.D0)*(A7*W(I,J,K-1)+A8*W(I,J,K)+A9*W(I,J,K-2))&
	+DMAX1(-WB,0.D0)*(A10*W(I,J,K)+A11*W(I,J,K-1)+A12*W(I,J,K+1)))/DZ1(K-1)

       FWC=0.D0
       VDX=(W(I+1,J,K)-W(I,J,K))/DX(I)-(W(I,J,K)-W(I-1,J,K))/DX(I-1)
       VDX=VDX/DX1(I-1)
       VDY=(W(I,J+1,K)-W(I,J,K))/DY(J)-(W(I,J,K)-W(I,J-1,K))/DY(J-1)
       VDY=VDY/DY1(J-1)
       VDZ=(W(I,J,K+1)-W(I,J,K))/DZ1(K)-(W(I,J,K)-W(I,J,K-1))/DZ1(K-1)
       VDZ=VDZ/DZ(K)
       VISZ=(VDX+VDY+VDZ)/RE
       FWW(I,J,K)=DT*(-FWX-FWY-FWZ-FWC+VISZ)

!       FWW(I,J,K)=(PRES(I,J,K)-PRES(I,J,K+1))/DZ1(K)-FWX-FWY-FWZ-FWC+VISZ
!             VISZ=((W(I+1,J,K)-2.*W(I,J,K)+W(I-1,J,K))/DX(I)/DX(I)
!     1            +(W(I,J+1,K)-2.*W(I,J,K)+W(I,J-1,K))/DY(J)/DY(J)
!     2            +(W(I,J,K+1)-2.*W(I,J,K)+W(I,J,K-1))/DZ1(K)/DZ1(K))/RE
    END DO 
  END DO
END DO
RETURN
END SUBROUTINE QUICK
!
!-------------------------------------------------------------
!
SUBROUTINE SOLVER
USE OMP_LIB
USE VARIABLES      
USE IBM
USE COEFFICIENT
IMPLICIT NONE
DOUBLE PRECISION :: DD,TS,TOL,UTEMP,VTEMP,WTEMP
DOUBLE PRECISION :: TEMP,DMAX,AREA_CYL
DOUBLE PRECISION :: AP,AE,AW,AN,AS,AT,AB
INTEGER TITE
INTEGER I,J,K,PPPP, II, JJ, KK	
CHARACTER(LEN=11) :: FILENAME,REFILE,LASTFILE	
FILENAME='CA00000.dat'
OPEN(UNIT=26,FILE='FLOW.dat')
OPEN(UNIT=39,FILE='COEFFICIENT.dat')
WRITE(39,*)'VARIABLES ="Time","Lw","CDx","CDy,"CDz"' 

OPEN(UNIT=71,FILE='TotalFx.dat')
OPEN(UNIT=72,FILE='TotalFy.dat')
OPEN(UNIT=73,FILE='TotalFz.dat')

OPEN(UNIT=74,FILE='CDx.dat')
OPEN(UNIT=75,FILE='CDy.dat')
OPEN(UNIT=76,FILE='CDz.dat')
!AREA_CYL= WIDTH
TIME=0.D0      
REIT=1
DO K=0,L,1
  DO J=0,M,1
    DO I=0,N,1  	        
    
       UN1(I,J,K)=U(I,J,K)		     
       VN1(I,J,K)=V(I,J,K)		      
       WN1(I,J,K)=W(I,J,K) 
       FUUN1(I,J,K)=0.D0
       FUUN2(I,J,K)=0.D0
       FVVN1(I,J,K)=0.D0
       FVVN2(I,J,K)=0.D0
       FWWN1(I,J,K)=0.D0
       FWWN2(I,J,K)=0.D0
!       Fx(I,J,K)=0.D0
!       Fy(I,J,K)=0.D0
!       Fz(I,J,K)=0.D0
    END DO
  END DO      
END DO
IF (RERUN == 1) THEN		
  LASTFILE='LAST.dat'
OPEN(UNIT=22,FILE=LASTFILE)		
  READ(22,*)TIME		
  READ(22,*)REIT
DO K=0,L,1
  DO J=0,M,1
    DO I=0,N,1  
    
      READ(22,*)U(I,J,K),V(I,J,K),W(I,J,K),PRES(I,J,K)
    END DO
  END DO
END DO  
DO K=0,L,1
  DO J=0,M,1
    DO I=0,N,1
      READ(22,*)FUUN1(I,J,K),FUUN2(I,J,K),FVVN1(I,J,K),FVVN2(I,J,K),FWWN1(I,J,K),FWWN2(I,J,K)
    END DO
  END DO
END DO            
 CLOSE(22)
REIT=REIT+1
END IF
TOL=1.0D-5     
!     CONSTRUCTION OF COEFFICIENT MATRIX
!
!
!     LOOP 10 :TIME STEPPING
!
DO TITE=REIT,TTT,1
  TIME=TIME+DT
  WRITE(*,*)'TIME =',TIME
!
!       PREDICTION
!       SOLVING VELOCITY FIELD
!       USING QUICK SCHEME
!
CALL QUICK
DO K=1,L-1,1
  DO J=1,M-1,1
    DO I=1,N-2,1  
    
      IF (REIT == 1) THEN
        U(I,J,K)=U(I,J,K)+DT*(PRES(I,J,K)-PRES(I+1,J,K))/DX1(I)+FUU(I,J,K)
      ELSE
        IF (REIT == 2) THEN
          U(I,J,K)=U(I,J,K)+DT*(PRES(I,J,K)-PRES(I+1,J,K))/DX1(I)+1.5*FUU(I,J,K)-0.5*FUUN1(I,J,K)
        ELSE
          U(I,J,K)=U(I,J,K)+DT*(PRES(I,J,K)-PRES(I+1,J,K))/DX1(I)&
                &+(23.*FUU(I,J,K)-16.*FUUN1(I,J,K)+5.*FUUN2(I,J,K))/12.
        ENDIF    
      ENDIF    
      FUUN2(I,J,K)=FUUN1(I,J,K)
      FUUN1(I,J,K)=FUU(I,J,K)                
    END DO 
  END DO
END DO

DO K=1,L-1,1
  DO J=1,M-2,1
    DO I=1,N-1,1  
    
      IF (REIT == 1) THEN
        V(I,J,K)=V(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J+1,K))/DY1(J)+FVV(I,J,K)
      ELSE
        IF (REIT == 2) THEN
          V(I,J,K)=V(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J+1,K))/DY1(J)+1.5*FVV(I,J,K)-0.5*FVVN1(I,J,K)
        ELSE
          V(I,J,K)=V(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J+1,K))/DY1(J)&
                &+(23.*FVV(I,J,K)-16.*FVVN1(I,J,K)+5.*FVVN2(I,J,K))/12.
        ENDIF    
      ENDIF    
      FVVN2(I,J,K)=FVVN1(I,J,K)
      FVVN1(I,J,K)=FVV(I,J,K)
    END DO       
  END DO
END DO

DO K=1,L-2,1
  DO J=1,M-1,1
    DO I=1,N-1,1  
    
      IF (REIT == 1) THEN
        W(I,J,K)=W(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J,K+1))/DZ1(K)+FWW(I,J,K)
      ELSE
        IF (REIT == 2) THEN
          W(I,J,K)=W(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J,K+1))/DZ1(K)+1.5*FWW(I,J,K)-0.5*FWWN1(I,J,K)
        ELSE
          W(I,J,K)=W(I,J,K)+DT*(PRES(I,J,K)-PRES(I,J,K+1))/DZ1(K)&
                &+(23.*FWW(I,J,K)-16.*FWWN1(I,J,K)+5.*FWWN2(I,J,K))/12.
        ENDIF    
      ENDIF    
      FWWN2(I,J,K)=FWWN1(I,J,K)
      FWWN1(I,J,K)=FWW(I,J,K)
    END DO 
  END DO
END DO
CALL BOUNDARY
!
!--------------------------------------------------------------
!     U,V,W AND P ---- PREDICTION-CORRECTION
!---------------------------------------------------------------
!
DO PPPP=1,ITER,1
! 	PRINT *,'TIME=',TIME,', ITERATIVE NUMBER =',PPPP
!   
!       PRESSURE CORRECTION 
!  
    DMAX=0.D0
  !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(U,V,W,PRES,DP,DX,DY,DZ,DX1,DY1,DZ1,DT,OMEGA)
           
  N_THREADS = OMP_GET_NUM_THREADS()
  THREAD_NUM = OMP_GET_THREAD_NUM()

  I_WORK_CHUNK = (N-1)/N_THREADS
  J_WORK_CHUNK = (M-1)/N_THREADS
  K_WORK_CHUNK = (L-1)/N_THREADS

!write(*,*) i_work_chunk,j_work_chunk,k_work_chunk
           
  I_LBOUND = I_WORK_CHUNK*THREAD_NUM + 1
  I_UBOUND = I_WORK_CHUNK*(THREAD_NUM + 1) 
  J_LBOUND = J_WORK_CHUNK*THREAD_NUM + 1
  J_UBOUND = J_WORK_CHUNK*(THREAD_NUM + 1)
  K_LBOUND = K_WORK_CHUNK*THREAD_NUM + 1
  K_UBOUND = K_WORK_CHUNK*(THREAD_NUM + 1)

  IF(THREAD_NUM==(N_THREADS-1)) THEN
    I_UBOUND = I_UBOUND + MOD((N-1),N_THREADS)
    J_UBOUND = J_UBOUND + MOD((M-1),N_THREADS)
    K_UBOUND = K_UBOUND + MOD((L-1),N_THREADS)
  END IF
!write(*,*) thread_num,"i: ",i_lbound,i_ubound
!write(*,*) thread_num,"j: ",j_lbound,j_ubound
!write(*,*) thread_num,"k: ",k_lbound,k_ubound

   !$OMP DO
!  DO I=1,N-1,1
!    DO J=1,M-1,1            
!      DO K=1,L-1,1
 ! !$OMP SIMD
    DO K=K_LBOUND,K_UBOUND,1
      DO J=J_LBOUND,J_UBOUND,1
        DO I=I_LBOUND,I_UBOUND,1    
      
      	 DD=(U(I,J,K)-U(I-1,J,K))*DY(J)*DZ(K)+(V(I,J,K)-V(I,J-1,K))*DX(I)*DY(J)+(W(I,J,K)-W(I,J,K-1))*DX(I)*DZ(K)
	 AP=-DY(J)*DZ(K)/DX1(I)-DY(J)*DZ(K)/DX1(I-1)-DX(I)*DZ(K)/DY1(J)-DX(I)*DZ(K)/DY1(J-1)-DX(I)*DY(J)/DZ1(K)&
         	&-DX(I)*DY(J)/DZ1(K-1)
	 AE=DY(J)*DZ(K)/DX1(I)
         AW=DY(J)*DZ(K)/DX1(I-1)
	 AN=DX(I)*DZ(K)/DY1(J)
         AS=DX(I)*DZ(K)/DY1(J-1)
         AT=DX(I)*DY(J)/DZ1(K)
         AB=DX(I)*DY(J)/DZ1(K-1)
         DP(I,J,K)=-AP*PRES(I,J,K)-AE*PRES(I+1,J,K)-AW*PRES(I-1,J,K)-AN*PRES(I,J+1,K)-AS*PRES(I,J-1,K)&
         	&-AT*PRES(I,J,K+1)-AB*PRES(I,J,K-1)+DD/DT
         DP(I,J,K)=DP(I,J,K)/AP
!         PRES(I,J,K)=PRES(I,J,K)+OMEGA*DP(I,J,K)
      END DO             
    END DO
  END DO
   !$OMP END DO

   !$OMP DO
!  DO I=1,N-1,1
!    DO J=1,M-1,1            
!      DO K=1,L-1,1
 ! !$OMP SIMD
    DO K=K_LBOUND,K_UBOUND,1
      DO J=J_LBOUND,J_UBOUND,1
        DO I=I_LBOUND,I_UBOUND,1    
     
         PRES(I,J,K)=PRES(I,J,K)+OMEGA*DP(I,J,K)
      END DO             
    END DO
  END DO
   !$OMP END DO
   !$OMP END PARALLEL
!89       CONTINUE
DO K=1,L-1,1
  DO J=1,M-1,1
    DO I=1,N-2,1  
    
       U(I,J,K)=U(I,J,K)-DT*(PRES(I+1,J,K)-PRES(I,J,K))/DX1(I)
    END DO
  END DO
END DO

DO K=1,L-1,1
  DO J=1,M-2,1
    DO I=1,N-1,1  
    
       V(I,J,K)=V(I,J,K)-DT*(PRES(I,J+1,K)-PRES(I,J,K))/DY1(J)
    END DO
  END DO
END DO

DO K=1,L-2,1
  DO J=1,M-1,1
    DO I=1,N-1,1  
    
       W(I,J,K)=W(I,J,K)-DT*(PRES(I,J,K+1)-PRES(I,J,K))/DZ1(K)
    END DO
  END DO
END DO
CALL BOUNDARY
!
  DO K=1,L-1,1
    DO J=1,M-1,1
      DO I=1,N-1,1                
      
         IF (DABS(DP(I,J,K)) > DMAX) THEN
            DMAX=DABS(DP(I,J,K))
              II=I
              JJ=J
              KK=K
         ENDIF
      END DO             
    END DO
  END DO
!------------------------------------------
	IF (DMAX < TOL) GOTO 89
END DO
!------------------------------------------
89      totalFx=0.D0      
        totalFy=0.D0
        totalFz=0.D0
 
 DO J=1,M
   A(J)=0.D0
   B(J)=0.D0
   C(J)=0.D0
 END DO

 DO K=1,L-1,1
   DO J=1,M-1,1
     DO I=1,N-2,1                
                             
        Fx(I,J,K)= -ETA(I,J,K)*( USOLID - U(I,J,K) )/DT
!        totalFx = totalFx + ( Fx(I,J,K)*( DX(I)*DY(J)*DZ(K) ) )             
     END DO  
   END DO
 END DO           
 DO K=1,L-1,1
   DO J=1,M-2,1
     DO I=1,N-1,1                       
        
        Fy(I,J,K) = -ETA(I,J,K)*( VSOLID - V(I,J,K) )/DT
!        totalFy = totalFy + ( Fy(I,J,K)*( DX(I)*DY(J)*DZ(K) ) )
     END DO   
   END DO        
 END DO
 DO K=1,L-2,1
   DO J=1,M-1,1
     DO I=1,N-1,1                     
     
        Fz(I,J,K) = -ETA(I,J,K)*( WSOLID - W(I,J,K) )/DT
!        totalFz = totalFz + ( Fz(I,J,K)*( DX(I)*DY(J)*DZ(K) ) ) 
     END DO
   END DO        
 END DO
!-----------------------------------------------------------------
! 	USING SIMPSON'S 1/3 RULE
!-----------------------------------------------------------------
!CALCULATE DRAG USING SIMPSON'S 1/3 RULE(X-FIRST)
 DO K=1,L-1
   DO J=1,M-1
     DO I=3,N-3,2	!CASE BY CASE N-3
	A(J)=A(J)+((X(I)-X(I-2))/6.D0)*(Fx(I-2,J,K)+Fx(I,J,K)+4.D0*Fx(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-1,2		
    totalFx=totalFx+((Y(J)-Y(J-2))/6.D0)*(A(J)+A(J-2)+4.D0*A(J-1))	
	WRITE(71,*) TIME,totalFx
 END DO 

!CALCULATE DRAG USING SIMPSON'S 1/3 RULE(Y-FIRST)
 DO K=1,L-2
   DO J=1,M-2
     DO I=3,N-1,2	
	B(J)=B(J)+((X(I)-X(I-2))/6.D0)*(Fy(I-2,J,K)+Fy(I,J,K)+4.D0*Fy(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-1,2		
    totalFy=totalFy+((Y(J)-Y(J-2))/6.D0)*(B(J)+B(J-2)+4.D0*B(J-1))	
	WRITE(72,*) TIME,totalFy
 END DO 

!CALCULATE DRAG USING SIMPSON'S 1/3 RULE(Z-FIRST)
 DO K=1,L-3
   DO J=1,M-3
     DO I=3,N+1,2	
	C(J)=C(J)+((X(I)-X(I-2))/6.D0)*(Fz(I-2,J,K)+Fz(I,J,K)+4.D0*Fz(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-5,2		
    totalFz=totalFz+((Y(J)-Y(J-2))/6.D0)*(C(J)+C(J-2)+4.D0*C(J-1))
	WRITE(73,*) TIME,totalFz	
 END DO 
!-----------------------------------------------------------------
! 	USING TRAPEZOIDAL RULE
!-----------------------------------------------------------------
!CALCULATE DRAG USING TRAPEZOIDAL RULE(X-FIRST)
 DO K=1,L-1
   DO J=1,M-1
     DO I=2,N-1,1	!CASE BY CASE N-3
	A(J)=A(J)+((X(I)-X(I-2))/2.D0)*(Fx(I-2,J,K)+Fx(I,J,K)+2.D0*Fx(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-1,2		
    totalFx=totalFx+((Y(J)-Y(J-2))/2.D0)*(A(J)+A(J-2)+2.D0*A(J-1))
	WRITE(71,*) TIME,totalFx	
 END DO 

!CALCULATE DRAG USING TRAPEZOIDAL RULE(Y-FIRST)
 DO K=1,L-2
   DO J=1,M-2
     DO I=3,N-1,2	
	B(J)=B(J)+((X(I)-X(I-2))/2.D0)*(Fy(I-2,J,K)+Fy(I,J,K)+2.D0*Fy(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-1,2		
    totalFy=totalFy+((Y(J)-Y(J-2))/2.D0)*(B(J)+B(J-2)+2.D0*B(J-1))	
	WRITE(72,*) TIME,totalFy
 END DO 

!CALCULATE DRAG USING TRAPEZOIDAL RULE(Z-FIRST)
 DO K=1,L-3
   DO J=1,M-3
     DO I=3,N+1,2	
	C(J)=C(J)+((X(I)-X(I-2))/2.D0)*(Fz(I-2,J,K)+Fz(I,J,K)+2.D0*Fz(I-1,J,K))	
     END DO                            
   END DO            
 END DO

 DO I=3,M-5,2		
    totalFz=totalFz+((Y(J)-Y(J-2))/2.D0)*(C(J)+C(J-2)+2.D0*C(J-1))	
	WRITE(73,*) TIME,totalFz
 END DO 
!----------------------------------------------------------------------  
!			UPDATE U AND V 
!----------------------------------------------------------------------  
    CDx=2.D0*totalFx		                             
    CDy=2.D0*totalFy		
    CDz=2.D0*totalFz
!   
 DO K=1,L-1,1
   DO J=1,M-1,1
     DO I=1,N-2,1      
     
        U(I,J,K) = U(I,J,K) + ETA(I,J,K)*( USOLID - U(I,J,K) ) 
     END DO                            
   END DO            
 END DO

 DO K=1,L-1,1
   DO J=1,M-2,1
     DO I=1,N-1,1             
     
        V(I,J,K) = V(I,J,K) + ETA(I,J,K)*( VSOLID - V(I,J,K) )
     END DO   
   END DO 
 END DO

 DO K=1,L-2,1
   DO J=1,M-1,1
     DO I=1,N-1,1             
     
        W(I,J,K) = W(I,J,K) + ETA(I,J,K)*( WSOLID - W(I,J,K) )
     END DO   
   END DO 
 END DO
!
!-----------------------------------------------------------------
!    COEFFICIENTS: CDi and Center_LW            
!-----------------------------------------------------------------
!    CDx=2.D0*totalFx/AREA_CYL                              
!    CDy=2.D0*totalFy/AREA_CYL
!    CDz=2.D0*totalFz/AREA_CYL
  
    CENTRAL_LW=0.D0                                   
    DO I=0,N,1                      ! DO I=101,N,1 ! it need to change by case : 100 = XC*(N-1)/LENGTH +1           
!      DO J=1,M,1 
        IF(U(I,(M+1)/2,(L+1)/2)*U(I+1,(M+1)/2,(L+1)/2)<0) THEN           
          XL = X(I)-(U(I,(M+1)/2,(L+1)/2)/(U(I+1,(M+1)/2,(L+1)/2)-U(I,(M+1)/2,(L+1)/2)))*(X(I+1)-X(I))           
         !YL=Y(J)-(V(I,J)/(V(I,J+1)-V(I,J)))*(Y(J+1)-Y(J))                                  
          CENTRAL_LW=(XL-XC-R)/(2.D0*R)                 
        END IF
!      END DO
    END DO 

	PHI=ASIN(U(0,0,0))*180.D0/PI

    IF (MOD(TITE,OTF) == 0) THEN 
       WRITE(39,*) TIME,CENTRAL_LW,CDx,CDy,CDz,PHI,CENTRAL_LW*2.D0
    END IF 

    IF ((TIME/PERIOD)>=22.D0 .AND. (TIME/PERIOD)>=28.D0) THEN 
       IF (MOD(TITE,OTF) == 0) THEN 
          WRITE(74,*) TIME/PERIOD,CDx
          WRITE(75,*) TIME/PERIOD,CDy
          WRITE(76,*) TIME/PERIOD,CDz
       END IF
    END IF
!
!-------------------------------------------------------------------                         
! 
  TEMP=0.D0   
  TS=0.D0
DO K=0,L,1
  DO J=0,M,1
    DO I=0,N,1   
    
       UTEMP=DABS(U(I,J,K)-UN1(I,J,K))
       VTEMP=DABS(V(I,J,K)-VN1(I,J,K))
       WTEMP=DABS(W(I,J,K)-WN1(I,J,K))
       TS=TS+UTEMP*UTEMP+VTEMP*VTEMP+WTEMP*WTEMP
       IF ((UTEMP>VTEMP).AND.(UTEMP>WTEMP)) THEN
	  IF(UTEMP > TEMP) TEMP=UTEMP
          ELSE IF(VTEMP > TEMP) THEN
		  TEMP=VTEMP
	  ELSE IF(WTEMP > TEMP) THEN
	          TEMP=WTEMP
       END IF
       UN1(I,J,K)=U(I,J,K)
       VN1(I,J,K)=V(I,J,K)
       WN1(I,J,K)=W(I,J,K)
    END DO 
  END DO
END DO
  TS=SQRT(TS/2./(N+1)/(M+1)/(L+1))
  WRITE(*,*)'ITERATION NUMBER =',PPPP ! SHOW ITERATION NUMBER 
!  WRITE(*,*)'DMAX=',DMAX,',I=',II,',J=',JJ,',K=',KK
  WRITE(*,*)'THE MAX. VARIATION=',TEMP,',SQRT OF VARIATION=',TS,',DMAX=',DMAX
  IF (MOD(TITE,FT) == 0) THEN
	WRITE(26,*)TITE,TIME,PPPP,TEMP,TS
  END IF
!-----------------------------------------------------------------------
!     OUTPUT NUMERICAL RESULTS
!-----------------------------------------------------------------------
  IF (MOD(TITE,OT) == 0) THEN
	WRITE(FILENAME(3:7),'(I5)')TITE/OT
	IF(TITE/OT < 10)   FILENAME(3:6)='0000'
	IF(TITE/OT < 100)  FILENAME(3:5)='000'
	IF(TITE/OT < 1000) FILENAME(3:4)='00'
	IF(TITE/OT < 10000)FILENAME(3:3)='0'
	REFILE='RERUN.dat'
  OPEN(UNIT=20,FILE=FILENAME)
!        WRITE(20,*)'TITLE="OSCILLATORY FLOW"'
!        WRITE(20,*)'VARIABLES=X,Y,Z,U,V,W,P'
!        WRITE(20,*)'ZONE T="OSCILLATORY FLOW", I=',N,', J=',M,', K=',L,',F=POINT'
	DO K=-1,L+1,1
          DO J=-1,M+1,1
            DO I=-1,N+1,1  
	WRITE(20,98)U(I,J,K),V(I,J,K),W(I,J,K),PRES(I,J,K) 
!	WRITE(20,98)(X(I)+X(I+1))/2.,(Y(J)+Y(J+1))/2.,(Z(K)+Z(K+1))/2.,U(I,J,K),V(I,J,K),W(I,J,K),PRES(I,J,K)
98 	FORMAT(4(E15.7,1X))
            END DO 
          END DO
        END DO
  OPEN(UNIT=21,FILE=REFILE)
	WRITE(21,*)TIME
	WRITE(21,*)REIT
        DO K=-1,L+1,1
          DO J=-1,M+1,1
	    DO I=-1,N+1,1
	  WRITE(21,*)U(I,J,K),UN1(I,J,K),V(I,J,K),VN1(I,J,K),W(I,J,K),WN1(I,J,K),PRES(I,J,K)
            END DO 
          END DO
        END DO
    CLOSE(21)
    CLOSE(20)
  ENDIF
END DO

RETURN
END SUBROUTINE SOLVER
