MODULE VARIABLES
IMPLICIT NONE
INTEGER, PARAMETER :: N=41,M=121
INTEGER, PARAMETER :: Imax1=85,Imax2=40,      Imax3=40,Imax4=85        ! total= all this +1
INTEGER, PARAMETER :: Jmax1=70,Jmax2=40,      Jmax3=40,Jmax4=70
DOUBLE PRECISION :: X(-1:N+1),Y(-1:M+1)
DOUBLE PRECISION :: U(0:N,0:M) 
DOUBLE PRECISION :: V(0:N,0:M)
DOUBLE PRECISION :: UN1(0:N,0:M)
DOUBLE PRECISION :: VN1(0:N,0:M)
DOUBLE PRECISION :: PRES(0:N,0:M)
DOUBLE PRECISION :: FUU(0:N,0:M),FVV(0:N,0:M)
DOUBLE PRECISION :: FUUN1(0:N,0:M),FVVN1(0:N,0:M)
DOUBLE PRECISION :: FUUN2(0:N,0:M),FVVN2(0:N,0:M)
DOUBLE PRECISION :: DX(-1:N),DX1(-1:N)
DOUBLE PRECISION :: DY(-1:M),DY1(-1:M)
DOUBLE PRECISION :: RE,DT,OMEGA,FREQ,TIME
DOUBLE PRECISION :: CD,CDP,CDV,CL,CLP,CLV
DOUBLE PRECISION :: lamda=0.78,pi=3.141592654
DOUBLE PRECISION :: DX_left,DX_right
DOUBLE PRECISION :: DY_upper,DY_lower
DOUBLE PRECISION :: RATIOx(-1:N+1),RATIOy(-1:M+1)

INTEGER :: TTT,OT,FT,CONDIF,RERUN,REIT,ITER
END MODULE VARIABLES


PROGRAM MESH 
USE VARIABLES
IMPLICIT NONE

DOUBLE PRECISION :: LENGTH,HEIGHT
DOUBLE PRECISION :: DXX, DYY
INTEGER :: I,J
         
      LENGTH=1.D0
      HEIGHT=3.D0
      DXX=LENGTH/(N*1.0)
      DYY=HEIGHT/(M*1.0)
      
!---------------------------------X--------------------------------------------------------------------------------------
OPEN(UNIT=11,FILE='XX.DAT')

DO I=-1,N+1
X(I)=DXX*I
WRITE(11,*)X(I),I
   END DO           
      CLOSE(11)

     
    OPEN(UNIT=12,FILE='YY.DAT')

  


  DO J=-1,M+1
    Y(J)=DYY*J
         WRITE(12,*)Y(J),J
      END DO
      CLOSE(12)




!--------------------------------------------------------------------------------------------------
 OPEN(UNIT=13,FILE='mesh.dat')
      WRITE(13,*)'VARIABLES=X,Y'
      WRITE(13,*)'ZONE I=',N,' ,J=',M
      DO  J=1,M
         DO  I=1,N
            WRITE(13,*)X(I),Y(J)
         END DO
      END DO
      CLOSE(13)
END PROGRAM MESH
