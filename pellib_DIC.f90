!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE getEl_DIC_M_S(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,PTLD_Vlcty)
           !Description: This subroutine builds elementary 
           		!matrices for Mass, Stiffness. 
           		!For other purposes, this could 
           		!be updated to compute decompositions
           		!for example, DIC for u0,g0 etc.
           !Arguments out
           		!ELM: Elementary Mass 
           		!ELK: Elementary stiffness 

IMPLICIT NONE
double precision ELCOR
integer NODEL,NDIME,NORDER
DIMENSION ELCOR(NODEL,NDIME)

INTEGER NGAUSS, IGAUSS, JGAUSS,I,J, IELEM
DOUBLE PRECISION WG(NGAUSS), XG(NGAUSS),&
                 XI, ETA, WTX, WTY, WTXY,&
                 X,Y
                 
double precision SF(NODEL), SFDL(NDIME,NODEL)
double precision JCBT(NDIME,NDIME),&
                 JCBTI(NDIME,NDIME), DETJ


double precision n(NDIME)
integer EDG_TYPE

!definitions for DIC,M,S
double precision SFDG(NDIME,NODEL),&
                 r_mod,K_W,Omega
double complex	 ELM(NODEL,NODEL),&
		 ELK(NODEL,NODEL),&
		 FZ_Phi,FZ_Vlcty,&
		 PTLD_Phi(NODEL),&
		 PTLD_Vlcty(NODEL)


  
  !Initialize elementary matrices and vectors
  ELM = dcmplx(0.0D0,0.0D0)
  ELK = dcmplx(0.0D0,0.0D0)
  PTLD_Phi = dcmplx(0.0D0,0.0D0)
  PTLD_Vlcty = dcmplx(0.0D0,0.0D0)
  
  !Get integration points.
  call GAULEG(NGAUSS, XG, NGAUSS, WG, NGAUSS)
  
  !Integrate inside domain.
    DO IGAUSS = 1,NGAUSS
      XI = XG(IGAUSS)
      WTX = WG(IGAUSS)
      DO JGAUSS = 1,NGAUSS
        ETA = XG(JGAUSS)
        WTY = WG(JGAUSS)
        WTXY = WTX*WTY
        !Get shape functions
        call getSF(XI, ETA, SF, NODEL, SFDL, NDIME,NORDER)
        !Get Jacobians.
          !Get JCBT = transpose(Jacobian)
          CALL MATMUL(SFDL, NDIME, NODEL, ELCOR, NODEL, NDIME, JCBT,&
          NDIME, NDIME, NDIME, NDIME, NODEL)
          !Invert JCBT
          CALL invJCBT(JCBT, NDIME, JCBTI, DETJ)


        !Update weights
        WTXY = WTXY*DETJ
        
        !Get global partial derivatives, i.e. (d/dx,d/dy) 
        call MATMUL(JCBTI,NDIME,NDIME,SFDL,NDIME,NODEL,SFDG,&
        	    NDIME,NODEL,NDIME,NODEL,NDIME)
        
        !Get (X,Y) from (XI,ETA)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        ENDDO        

        !Build Mass and Stiffness
        DO I = 1,NODEL
          DO J = 1,NODEL
          ELM(J,I) = ELM(J,I) + SF(J)*SF(I)*WTXY
          ELK(J,I) = ELK(J,I) + (SFDG(1,J)*SFDG(1,I) + &
                                 SFDG(2,J)*SFDG(2,I))*WTXY
          endDO
        endDO

        !Build RHS for DIC
          !Calculate load          
          r_mod = sqrt(X*X + Y*Y)
          FZ_Phi = CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
					 - Omega*(0.0D0)))
          FZ_Vlcty = CDEXP(DCMPLX(0.0D0,1.0D0)*(K_W*r_mod&
					 - Omega*(0.0D0)))&
			*(DCMPLX(0.0D0,-1.0D0)*Omega)
					 
          !Calculate inner product of load with weight functions
          DO J = 1,NODEL
            PTLD_Phi(J) = PTLD_Phi(J) + FZ_Phi*SF(J)*WTXY
            PTLD_Vlcty(J) = PTLD_Vlcty(J) + FZ_Vlcty*SF(J)*WTXY
          endDO
        
      endDO
    endDO
 
  
end SUBROUTINE getEl_DIC_M_S
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
