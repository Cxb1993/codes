MODULE VARIABLES
IMPLICIT NONE
DOUBLE PRECISION :: t1,t2
DOUBLE PRECISION :: CLOCK
INTEGER, PARAMETER :: N=41,M=41
DOUBLE PRECISION :: X(-1:N+1),Y(-1:M+1)
DOUBLE PRECISION :: U(-1:N+1,-1:M+1) 
DOUBLE PRECISION :: V(-1:N+1,-1:M+1)
DOUBLE PRECISION :: UN1(0:N,0:M)
DOUBLE PRECISION :: VN1(0:N,0:M)
DOUBLE PRECISION :: PRES(-1:N+1,-1:M+1)
DOUBLE PRECISION :: ST(0:N, 0:M), VORT(0:N,0:M)
DOUBLE PRECISION :: FUU(0:N,0:M),FVV(0:N,0:M)
DOUBLE PRECISION :: FUUN1(0:N,0:M),FVVN1(0:N,0:M)
DOUBLE PRECISION :: FUUN2(0:N,0:M),FVVN2(0:N,0:M)
DOUBLE PRECISION :: DP(0:N,0:M)
DOUBLE PRECISION :: DX(-1:N),DX1(-1:N)
DOUBLE PRECISION :: DY(-1:M),DY1(-1:M)
DOUBLE PRECISION :: RE,DT,OMEGA,TIME
DOUBLE PRECISION :: PI  
INTEGER :: TTT,OT,FT,RERUN,REIT,ITER
END MODULE VARIABLES


PROGRAM FLOW
!
!       THIS PROGRAM IS FOR SOLVING UNIFORM FLOW PAST A SQUARE CYLINDER 
!       USING SOLA METHOD
!       
!       AUTHOR: M. J. CHERN   DATE:4 JULY,2001
!       
!       
!       STAGGERED GRID (MAC)
!
!       NON-UNIFORM GRID VERSION 
!         DATE: 14 MAY, 2002
!	
!
!        QUICK SCHEME
!        DATE: 18 SEPT., 2002
!
!
!       PROGRAM MODIFIED FOR AN OSCIALLATING FLOW PAST A SQUARE CYLINDER
!       DATE: 26 OCt. 2002
!
!       MODIFY BOUNDARY CONDITIONS 
!       DATE: 26 Nov. 2002
!
!       CONVERT TO FORTRAN 90 FORMAT
!       DATE: 20 SEPTEMBER 2011
!
!       NATURAL CONVECTION
!       DATE: 4 JAN. 2015
!     
!       SOLA Method is replaced by projection method
!       This program is for cavity flow simualtions.
!       DATE: 14 May, 2015
USE VARIABLES
IMPLICIT NONE
!
!       INPUT DATA OF THE FUNCTION U(X,Y) ON COLLOCATION POINTS
 OPEN(UNIT=43,FILE='CLOCK.DAT')
  
  CALL cpu_time(t1)


	CALL INPUT
	PRINT *,'INPUT PASS'
!
    CALL SOLVER  

    CALL cpu_time(t2)
    WRITE(43,*) 'Time taken by code was',t2-t1,'seconds'
    PRINT *,'SOLVER PASS'

CONTAINS
!-------------------------------------------------------------
!
!
!
	SUBROUTINE INPUT
	USE VARIABLES
        IMPLICIT NONE
	INTEGER :: I,J

	OPEN(UNIT=15,FILE='XX.DAT')
        DO I=-1,N+1,1
           READ(15,*)X(I)
        END DO
	CLOSE(15)

	OPEN(UNIT=16,FILE='YY.DAT')
        DO J=-1,M+1,1
           READ(16,*)Y(J)
        END DO
	CLOSE(16)

        DO I=-1,N
         DX(I)=X(I+1)-X(I)
        END DO
        DX1(0)=0.5*(X(2)-X(1))
        DO I=-1,N-1
           DX1(I)=0.5*(X(I+2)-X(I))
        END DO

        DO J=-1,M
           DY(J)=Y(J+1)-Y(J)
        END DO

        DO J=-1,M-1
           DY1(J)=0.5*(Y(J+2)-Y(J))

        END DO

	OPEN(UNIT=19,FILE='PARA.F')
        READ(19,*)TTT,OT,DT,FT,ITER,RE,OMEGA,RERUN
	CLOSE(19)
        TIME=0.D0    
	DO I=-1,N+1,1
	   DO  J=-1,M+1,1
           U(I,J)=0.D0          
           V(I,J)=0.D0
           PRES(I,J)=0.D0	   
       END DO
    END DO
	RETURN
	END SUBROUTINE INPUT
!-------------------------------------------------------------
!
!
!
	SUBROUTINE BOUNDARY
	USE VARIABLES
    IMPLICIT NONE
	INTEGER :: I,J
!       DEFINE U(1,Y),V(1,Y)

DO J=0,M,1    
U(N-1,J)=0.D0                                 
U(N,J)=0.D0
V(N,J)=-V(N-1,J)
V(N+1,J)=V(N,J)
END DO
!
!       DEFINE U(0,Y),V(0,Y)
!
DO J=0,M,1
U(0,J)=0.D0
V(0,J)=-V(1,J)                                     
U(-1,J)=0.D0
V(-1,J)=V(0,J)
END DO
!
!       DEFINE U(X,0),V(X,0)
!
DO I=1,N
U(I,0)=-U(I,1)
V(I,0)=0.D0                                    
U(I,-1)=U(I,0)
V(I,-1)=0.D0
END DO
!
!       DEFINE U(X,1),V(X,1)
!
DO I=0,N,1
U(I,M)=2.0-U(I,M-1)
V(I,M-1)=0.D0
U(I,M+1)=U(I,M)
V(I,M)=V(I,M-1)
END DO

	RETURN
	END SUBROUTINE BOUNDARY
!-------------------------------------------------------------
!
!
!
SUBROUTINE PRESSURE_BOUNDARY
USE VARIABLES
IMPLICIT NONE
INTEGER :: I,J

DO J=0,M,1
PRES(0,J)=PRES(1,J)
PRES(N,J)=PRES(N-1,J)
END DO

DO I=1,N
PRES(I,0)=PRES(I,1)
PRES(I,M)=PRES(I,M-1)                                    
END DO

RETURN
END SUBROUTINE PRESSURE_BOUNDARY

!-------------------------------------------------------------
!
!
!
        SUBROUTINE QUICK
        USE VARIABLES
        IMPLICIT NONE
        DOUBLE PRECISION :: FUX,FUY,FUC,FVX,FVY,FVC
        DOUBLE PRECISION :: VISX,VISY,VDX,VDY
        DOUBLE PRECISION :: UE,UE1,UNN,UW,UW1,US1
        DOUBLE PRECISION :: VN,VE1,VNN,VS,VW1,VS1
        DOUBLE PRECISION :: X1,X2,X3,X4,X5,X6,X7,X8
        DOUBLE PRECISION :: Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
        DOUBLE PRECISION :: A1,A2,A3
        INTEGER :: I,J
        DO  I=1,N-2,1
           DO  J=1,M-1,1            
              X1=DX(I)/2.
              X2=DX(I+1)/2.
              X3=DX(I+1)+DX(I+2)/2.
              X4=DX(I-1)+DX(I)/2.
              X5=DX(I)/2.
              X6=DX(I-1)/2.
              X7=DX(I-1)+DX(I-2)/2.
              X8=DX(I)+DX(I+1)/2.
              UE=0.5*(DX(I)*U(I+1,J)/DX1(I)+DX(I+1)*U(I,J)/DX1(I))
              IF (UE > 0.0) THEN
                 A1=X2*X4/(X4-X1)/(X1+X2)
                 A2=X1*X4/(X1+X2)/(X2+X4)
                 A3=-X1*X2/(X4-X1)/(X2+X4)
                 UE1=A1*U(I,J)+A2*U(I+1,J)+A3*U(I-1,J)
              ELSE
                 A1=-X1*X3/(X1+X2)/(X2-X3)
                 A2=X2*X3/(X1+X3)/(X1+X2)
                 A3=X1*X2/(X1+X3)/(X2-X3)
                 UE1=A1*U(I+1,J)+A2*U(I,J)+A3*U(I+2,J)
              ENDIF
              UW=0.5*(DX(I)*U(I-1,J)/DX1(I-1)+DX(I-1)*U(I,J)/DX1(I-1))
              A1=X6*X8/(X5+X6)/(X8-X5)
              A2=X5*X8/(X5+X6)/(X6+X8)
              A3=-X5*X6/(X6+X8)/(X8-X5)
              IF (UW > 0.0) THEN
                     A1=X5*X7/(X7-X6)/(X5+X6)
                     A2=X6*X7/(X5+X7)/(X5+X6)
                     A3=-X5*X6/(X7-X6)/(X5+X7)
                     UW1=A1*U(I-1,J)+A2*U(I,J)+A3*U(I-2,J)
              ELSE
                     A1=X6*X8/(X5+X6)/(X8-X5)
                     A2=X5*X8/(X5+X6)/(X6+X8)
                     A3=-X5*X6/(X6+X8)/(X8-X5)
                     UW1=A1*U(I,J)+A2*U(I-1,J)+A3*U(I+1,J)
              ENDIF
              FUX=(UE*UE1-UW*UW1)/DX(I)
               Y1=DY(J)/2.
               Y2=DY(J)/2.
               Y3=DY(J+1)+DY(J)/2.
               Y4=DY(J-1)+DY(J)/2.
               Y5=DY(J-1)/2.
               Y6=DY(J-1)/2.
               Y7=DY(J-2)+DY(J-1)/2.
               Y8=DY(J)+DY(J-1)/2.
               VN=0.5*(V(I,J)+V(I+1,J))
              IF (VN > 0.D0) THEN
                 A1=Y2*Y4/(Y4-Y1)/(Y1+Y2)
                 A2=Y1*Y4/(Y1+Y2)/(Y2+Y4)
                 A3=-Y1*Y2/(Y4-Y1)/(Y2+Y4)
                 UNN=A1*U(I,J)+A2*U(I,J+1)+A3*U(I,J-1)
              ELSE
                A1=-Y1*Y3/(Y1+Y2)/(Y2-Y3)
                A2=Y2*Y3/(Y1+Y3)/(Y1+Y2)
                A3=Y1*Y2/(Y1+Y3)/(Y2-Y3)
                UNN=A1*U(I,J+1)+A2*U(I,J)+A3*U(I,J+2)
              ENDIF
              VS=0.5*(V(I,J-1)+V(I+1,J-1))
              IF (VS > 0.D0) THEN
                    A1=Y5*Y7/(Y7-Y6)/(Y5+Y6)
                    A2=Y6*Y7/(Y5+Y7)/(Y5+Y6)
                    A3=-Y5*Y6/(Y7-Y6)/(Y5+Y7)
                    US1=A1*U(I,J-1)+A2*U(I,J)+A3*U(I,J-2)
              ELSE
                    A1=Y6*Y8/(Y5+Y6)/(Y8-Y5)
                    A2=Y5*Y8/(Y5+Y6)/(Y6+Y8)
                    A3=-Y5*Y6/(Y6+Y8)/(Y8-Y5)
                    US1=A1*U(I,J)+A2*U(I,J-1)+A3*U(I,J+1)
              ENDIF
             FUY=(VN*UNN-VS*US1)/DY1(J-1)
             FUC=0.D0 
             VDX=(U(I+1,J)-U(I,J))/DX(I)-(U(I,J)-U(I-1,J))/DX(I-1)
             VDX=VDX/DX1(I-1)
             VDY=(U(I,J+1)-U(I,J))/DY1(J)-(U(I,J)-U(I,J-1))/DY1(J-1)
             VDY=VDY/DY(J)
             VISX=(VDX+VDY)/RE
             FUU(I,J)=DT*(-FUX-FUY-FUC+VISX)
	       END DO
       	    END DO   


        DO  I=1,N-1,1
           DO  J=1,M-2,1             
              X1=DX(I)/2.                               
              X2=DX(I)/2.                                
              X3=DX(I+1)+DX(I)/2.                         
              X4=DX(I-1)+DX(I)/2.                            
              X5=DX(I-1)/2.                                     
              X6=DX(I-1)/2.                              
              X7=DX(I-2)+DX(I-1)/2.                       
              X8=DX(I)+DX(I-1)/2.                        
              UE=0.5*(U(I,J)+U(I,J+1))
              IF (UE > 0.D0) THEN
                  A1=X2*X4/(X4-X1)/(X1+X2)
                  A2=X1*X4/(X1+X2)/(X2+X4)
                  A3=-X1*X2/(X4-X1)/(X2+X4)
                  VE1=A1*V(I,J)+A2*V(I+1,J)+A3*V(I-1,J)
              ELSE
                  A1=-X1*X3/(X1+X2)/(X2-X3)
                  A2=X2*X3/(X1+X3)/(X1+X2)
                  A3=X1*X2/(X1+X3)/(X2-X3)
                  VE1=A1*V(I+1,J)+A2*V(I,J)+A3*V(I+2,J)
              ENDIF
              UW=0.5*(U(I-1,J)+U(I-1,J+1))
              IF (UW > 0.D0) THEN
                  A1=X5*X7/(X7-X6)/(X5+X6)
                  A2=X6*X7/(X5+X7)/(X5+X6)
                  A3=-X5*X6/(X7-X6)/(X5+X7)                 
                  VW1=A1*V(I-1,J)+A2*V(I,J)+A3*V(I-2,J)                 
              ELSE
                  A1=X6*X8/(X5+X6)/(X8-X5)
                  A2=X5*X8/(X5+X6)/(X6+X8)             
                  A3=-X5*X6/(X6+X8)/(X8-X5)
                  VW1=A1*V(I,J)+A2*V(I-1,J)+A3*V(I+1,J)
              ENDIF
              FVX=(UE*VE1-UW*VW1)/DX1(I-1)                                           
              Y1=DY(J)/2.
              Y2=DY(J+1)/2.
              Y3=DY(J+1)+DY(J+2)/2.
              Y4=DY(J-1)+DY(J)/2.                                                  
              Y5=DY(J)/2.
              Y6=DY(J-1)/2.
              Y7=DY(J-1)+DY(J-2)/2.
              Y8=DY(J)+DY(J+1)/2.
              VN=0.5*(DY(J+1)*V(I,J)/DY1(J)+DY(J)*V(I,J+1)/DY1(J))
              IF (VN > 0.D0) THEN
                 A1=Y2*Y4/(Y4-Y1)/(Y1+Y2)
                 A2=Y1*Y4/(Y1+Y2)/(Y2+Y4)
                 A3=-Y1*Y2/(Y4-Y1)/(Y2+Y4)
                 VNN=A1*V(I,J)+A2*V(I,J+1)+A3*V(I,J-1)
              ELSE
                 A1=-Y1*Y3/(Y1+Y2)/(Y2-Y3)
                 A2=Y2*Y3/(Y1+Y3)/(Y1+Y2)
                 A3=Y1*Y2/(Y1+Y3)/(Y2-Y3)
                 VNN=A1*V(I,J+1)+A2*V(I,J)+A3*V(I,J+2)
              ENDIF
              VS=0.5*(DY(J)*V(I,J-1)/DY1(J-1)+DY(J-1)*V(I,J)/DY1(J-1))
              IF (VS > 0.D0) THEN
                    A1=Y5*Y7/(Y7-Y6)/(Y5+Y6)
                    A2=Y6*Y7/(Y5+Y7)/(Y5+Y6)
                    A3=-Y5*Y6/(Y7-Y6)/(Y5+Y7)
                    VS1=A1*V(I,J-1)+A2*V(I,J)+A3*V(I,J-2)
              ELSE
                    A1=Y6*Y8/(Y5+Y6)/(Y8-Y5)
                    A2=Y5*Y8/(Y5+Y6)/(Y6+Y8)
                    A3=-Y5*Y6/(Y6+Y8)/(Y8-Y5)
                    VS1=A1*V(I,J)+A2*V(I,J-1)+A3*V(I,J+1)
              ENDIF
             FVY=(VN*VNN-VS*VS1)/DY(J)                                                
             FVC=0.D0                                                                 
             VDX=(V(I+1,J)-V(I,J))/DX(I)-(V(I,J)-V(I-1,J))/DX(I-1)
             VDX=VDX/DX1(I-1)
             VDY=(V(I,J+1)-V(I,J))/DY1(J)-(V(I,J)-V(I,J-1))/DY1(J-1)
             VDY=VDY/DY(J)
             VISY=(VDX+VDY)/RE                                                         
             FVV(I,J)=DT*(-FVX-FVY-FVC+VISY)
             
 	      END DO
 	    END DO  

       
        RETURN
        END SUBROUTINE QUICK
!
!
!
!-------------------------------------------------------------
    SUBROUTINE SOLVER
    USE VARIABLES
    IMPLICIT NONE
	DOUBLE PRECISION :: DD
	DOUBLE PRECISION :: TEMP,DMAX
    DOUBLE PRECISION :: UTEMP, VTEMP
    DOUBLE PRECISION :: TS, TOL
    DOUBLE PRECISION :: AP, AN, AS, AW, AE
	INTEGER :: TITE
	INTEGER :: I,J,PPPP, II, JJ
	CHARACTER(LEN=11) :: FILENAME,REFILE,LASTFILE
  
	FILENAME='CA00000.DAT'                                          
   
    TIME=0.D0      
	REIT=1
    DO  I=0,N,1
       DO  J=0,M,1
		 UN1(I,J)=U(I,J)
		 VN1(I,J)=V(I,J)
         FUUN1(I,J)=0.D0
         FUUN2(I,J)=0.D0
         FVVN1(I,J)=0.D0
         FVVN2(I,J)=0.D0
	   END DO   
	END DO
	IF (RERUN == 1) THEN
		LASTFILE='LAST.DAT'
		OPEN(UNIT=22,FILE=LASTFILE)
		READ(22,*)TIME
		READ(22,*)REIT
		DO  I=0,N,1
		DO  J=0,M,1
		   READ(22,*)U(I,J),V(I,J),PRES(I,J)  
		END DO
	    END DO
        DO  J=0,M,1
		DO  I=0,N,1
		  READ(22,*)FUUN1(I,J),FUUN2(I,J),FVVN1(I,J),FVVN2(I,J)
        END DO        
        END DO
		CLOSE(22)
        REIT=REIT+1
	ENDIF
	TOL=1.0D-4
        OPEN(UNIT=26,FILE='flow.dat')              
!
!
!       LOOP 10 :TIME STEPPING
!
	DO  TITE=REIT,TTT,1
	TIME=TIME+DT
	PRINT *,'TIME =',TIME
!
!       PREDICTION
!       SOLVING VELOCITY FIELD
!       USING QUICK SCHEME
!
    CALL QUICK

        DO  I=1,N-2,1
           DO  J=1,M-1,1
             IF (REIT == 1) THEN
               U(I,J)=U(I,J)+FUU(I,J)
             ELSE
               IF(REIT == 2) THEN
                 U(I,J)=U(I,J)+1.5*FUU(I,J)-0.5*FUUN1(I,J)
               ELSE
                 U(I,J)=U(I,J)+(23.*FUU(I,J)-16.*FUUN1(I,J)+5.*FUUN2(I,J))/12.
               ENDIF    
             ENDIF    
             FUUN2(I,J)=FUUN1(I,J)
             FUUN1(I,J)=FUU(I,J)
	       END DO   
	    END DO  
        DO  I=1,N-1,1
           DO  J=1,M-2,1 
             IF (REIT == 1) THEN            
		        V(I,J)=V(I,J)+FVV(I,J)
             ELSE
               IF (REIT == 2) THEN
                 V(I,J)=V(I,J)+1.5*FVV(I,J)-0.5*FVVN1(I,J)
               ELSE
                 V(I,J)=V(I,J)+(23.*FVV(I,J)-16.*FVVN1(I,J)+5.*FVVN2(I,J))/12.
               ENDIF    
             ENDIF        
             FVVN2(I,J)=FVVN1(I,J)
             FVVN1(I,J)=FVV(I,J)         
           END DO   
	    END DO
        CALL BOUNDARY
!
!     U,V,W AND P ---- PROJECTION METHOD
!
	DO PPPP=1,ITER,1
!
!       PRESSURE CORRECTION 
!
	DMAX=0.D0

	DO  I=1,N-1,1
	   DO  J=1,M-1,1             
	      DD=(U(I,J)-U(I-1,J))*DY(J)+(V(I,J)-V(I,J-1))*DX(I)
          AP=-DY(J)/DX1(I)-DY(J)/DX1(I-1)-DX(I)/DY1(J)-DX(I)/DY1(J-1)
          AE=DY(J)/DX1(I)
          AW=DY(J)/DX1(I-1)
          AN=DX(I)/DY1(J)
          AS=DX(I)/DY(J-1)
	      DP(I,J)=-AP*PRES(I,J)-AE*PRES(I+1,J)-AW*PRES(I-1,J)-AN*PRES(I,J+1)-AS*PRES(I,J-1)+DD/DT
          DP(I,J)=DP(I,J)/AP
           PRES(I,J)=PRES(I,J)+OMEGA*DP(I,J)
                IF (DABS(DP(I,J)) > DMAX) THEN
                     DMAX=DABS(DP(I,J))
                     II=I
                     JJ=J
	            ENDIF             
	    END DO
      END DO
      CALL PRESSURE_BOUNDARY
	  IF (DMAX < TOL) GOTO 89
      END DO

89       CONTINUE
DO  I=1,N-2,1
    DO  J=1,M-1,1
        U(I,J)=U(I,J)-DT*(PRES(I+1,J)-PRES(I,J))/DX1(I)
    END DO
END DO

DO  I=1,N-1,1
    DO  J=1,M-2,1
        V(I,J)=V(I,J)-DT*(PRES(I,J+1)-PRES(I,J))/DY1(J)
    END DO
END DO

CALL BOUNDARY

              
         
!---------------------------------------------------------------------------------------------


      TEMP=0.D0   
      TS=0.D0
	DO  I=0,N,1
	   DO  J=0,M,1
		 UTEMP=DABS(U(I,J)-UN1(I,J))
		 VTEMP=DABS(V(I,J)-VN1(I,J))
		 TS=TS+UTEMP*UTEMP+VTEMP*VTEMP
		 IF (UTEMP > VTEMP) THEN
		      IF(UTEMP > TEMP) TEMP=UTEMP
		 ELSE IF(VTEMP > TEMP) THEN
		      TEMP=VTEMP
		 ENDIF
		 UN1(I,J)=U(I,J)
		 VN1(I,J)=V(I,J)         
       END DO
    END DO    
      TS=SQRT(TS/2./(N+1)/(M+1))
      WRITE(*,*)'UV-P ITERATION NUMBER =',PPPP
      WRITE(*,*)'DMAX= ',DMAX,' ,I= ',II,' ,J= ',JJ
      WRITE(*,*)'THE MAX. VARIATION=', TEMP,' ,SQRT OF VARIATION=',TS
      IF (MOD(TITE,FT) == 0) THEN
         WRITE(26,*)TITE,TIME,PPPP,TEMP,TS
      ENDIF
!
!     OUTPUT NUMERICAL RESULTS
!
		IF (MOD(TITE,OT) == 0) THEN
                WRITE(FILENAME(3:7),'(I5)')TITE/OT
                IF(TITE/OT .LT. 10)   FILENAME(3:6)='0000'
                IF(TITE/OT .LT. 100)  FILENAME(3:5)='000'
                IF(TITE/OT .LT. 1000) FILENAME(3:4)='00'
                IF(TITE/OT .LT. 10000)FILENAME(3:3)='0'
		REFILE='RERUN.DAT'
		OPEN(UNIT=20,FILE=FILENAME)
                DO I=0,N
                   DO J = 0,M
                     ST(I,J)=0.D0
                     VORT(I,J)=0.D0
                   END DO
                END DO
                DO I=1,N-1
                  ST(I,0)=0.D0
                  VORT(I,J)=0.D0
                  DO J=1,M-1
                     ST(I,J)=ST(I,J-1)+(Y(J+1)-Y(J))*U(I,J)
                     VORT(I,J)=-(U(I-1,J)-U(I-1,J-1))/(0.5*(Y(J+1)-Y(J-1)))+(V(I,J-1)-V(I-1,J-1))/(0.5*(X(I+1)-X(I-1)))
                  END DO
                END DO
                WRITE(20,*)'TITLE="CAVITY FLOW"'
                WRITE(20,*)'VARIABLES=X,Y,U,V,P,STREAM,VORTICITY'
                WRITE(20,*)'ZONE T="CAVITY FLOW", I=',N,' , J=',M,' F=POINT'
		DO  J=0,M-1,1
		DO  I=0,N-1,1  
		  WRITE(20,98)(X(I)+X(I+1))/2.,(Y(J)+Y(J+1))/2., U(I,J),V(I,J),PRES(I,J),ST(I,J),VORT(I,J)
98                FORMAT(7(E15.7,1X))
                END DO
                END DO
		OPEN(UNIT=21,FILE=REFILE)
		WRITE(21,*)TIME
		WRITE(21,*)TITE
		DO  J=0,M,1
		DO  I=0,N,1
		  WRITE(21,*)U(I,J),V(I,J),PRES(I,J)
                END DO        
                END DO
        DO  J=0,M,1
		DO  I=0,N,1
		  WRITE(21,*)FUUN1(I,J),FUUN2(I,J),FVVN1(I,J),FVVN2(I,J)
        END DO        
        END DO
		CLOSE(21)
	CLOSE(20)
	ENDIF
    END DO
    
    CLOSE(26)
    RETURN
    END SUBROUTINE SOLVER        
END PROGRAM FLOW
