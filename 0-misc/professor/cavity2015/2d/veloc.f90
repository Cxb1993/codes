!
!
!       AUTHOR:M J CHERN      17 MAY,1996
        PROGRAM STREAM
        PARAMETER (N=401,M=201)
        DOUBLE PRECISION :: U(0:N,0:M),V(0:N,0:M),P(0:N,0:M), TEM(0:N,0:M)
        DOUBLE PRECISION :: Fx(0:N,0:M),Fy(0:N,0:M)                           ! MY PLUS
        DOUBLE PRECISION :: QF(0:N,0:M)
        DOUBLE PRECISION :: Y(-1:M+1),X(-1:N+1)
        DOUBLE PRECISION :: PI
        INTEGER I,J
        PI=DCOS(-1.D0)
        OPEN(UNIT=15,FILE='CA00008.DAT')
        OPEN(UNIT=16,FILE='YY.DAT')
        OPEN(UNIT=17,FILE='VELRE100.DAT')
        OPEN(UNIT=18,FILE='XX.DAT')
        OPEN(UNIT=19,FILE='FX.DAT')
        OPEN(UNIT=20,FILE='FY.DAT')
        OPEN(UNIT=21,FILE='FXY.DAT')
        OPEN(UNIT=22,FILE='Q.DAT')
        DO J=0,M,1
                DO  I=0,N,1
                  READ(15,*) U(I,J),V(I,J),P(I,J),TEM(I,J),FX(I,J),FY(I,J),QF(I,J)
                END DO
       END DO
        DO 300 J=-1,M+1,1
                READ(16,*)Y(J)
300     CONTINUE
        DO 800 I=-1,N+1,1
                READ(18,*)X(I)
800     CONTINUE
        WRITE(17,*)'TITLE="VELOCITY VECTOR FIELD"'
        WRITE(17,*)'VARIABLES=X,Y,U,V,P, TEM'
        WRITE(17,*)'ZONE T="XYZ", I=',N,' ,J=',M,' ,F=POINT'
        DO 600 J=0,M-1,1
                DO 700 I=0,N-1,1
        WRITE(17,98)(X(I)+X(I+1))/2.,(Y(J)+Y(J+1))/2.,U(I,J),V(I,J),P(I,J), TEM(I,J)
700             CONTINUE
600     CONTINUE

        WRITE(19,*)'TITLE="FX FIELD"'
        WRITE(19,*)'VARIABLES=X,Y,FX'
        WRITE(19,*)'ZONE T="XYZ", I=',N-1,' ,J=',M-1,' ,F=POINT'
        DO 900 J=1,M-1,1
                DO 1000 I=1,N-1,1
        WRITE(19,99)X(I+1),(Y(J)+Y(J+1))/2.,FX(I,J)
1000             CONTINUE
900     CONTINUE
        WRITE(20,*)'TITLE="FY FIELD"'
        WRITE(20,*)'VARIABLES=X,Y,FY'
        WRITE(20,*)'ZONE T="XYZ", I=',N-1,' ,J=',M-1,' ,F=POINT'
        DO J=1,M-1,1
           DO I=1,N-1,1
            WRITE(20,99)(X(I)+X(I+1))/2., Y(J+1), FY(I,J)
           END DO
        END DO
        WRITE(21,*)'TITLE="FXY FIELD"'
        WRITE(21,*)'VARIABLES=X,Y,FX,FY'
        WRITE(21,*)'ZONE T="XYZ", I=',N-1,' ,J=',M-1,' ,F=POINT'
        DO J=1,M-1,1
           DO I=1,N-1,1
            WRITE(21,97)(X(I)+X(I+1))/2.,(Y(J)+Y(J+1))/2., FX(I,J), FY(I,J)
           END DO
        END DO

        WRITE(22,*)'TITLE="Q FIELD"'
        WRITE(22,*)'VARIABLES=X,Y,q'
        WRITE(22,*)'ZONE T="XYZ", I=',N-1,' ,J=',M-1,' ,F=POINT'
        DO J=1,M-1,1
            DO I=1,N-1,1
               WRITE(22,99)(X(I)+X(I+1))/2.,(Y(J)+Y(J+1))/2.+5.0, QF(I,J)
            END DO
        END DO
98      FORMAT(6(E15.7,1X))
99      FORMAT(3(E15.7,1X))
97      FORMAT(4(E15.7,1X))
        CLOSE(17)
        CLOSE(16)
        CLOSE(15)
        CLOSE(19)
        CLOSE(20)
        CLOSE(21)
        CLOSE(22)
        END PROGRAM STREAM
