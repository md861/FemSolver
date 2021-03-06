!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getLn_norm(IStep,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                   TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&
                   ELCOR,NGAUSS,&
                   K_W,Omega,dt,Cp,error1,error2)
  
IMPLICIT NONE
  integer IStep,TOTELS,ofst_ELDAT,NODEL,TOTCOORD,NDIME,NORDER
  double precision COORD(TOTCOORD,NDIME),ELCOR(NODEL,NDIME)
  integer ELNODDAT(NODEL),ELNDS(NODEL)
  
  integer ELDAT(TOTELS, ofst_ELDAT + NODEL)              

  integer TOTDOF,NANGL,NODEF,&
          GLDF(TOTELS,NODEL),ELDF(NODEL*NANGL)
  double complex Phi(TOTDOF),&
                 PhiNum,PhiExa


INTEGER NGAUSS
DOUBLE PRECISION WG(NGAUSS), XG(NGAUSS),&
                 XI, ETA, WTX, WTY, WTXY,&
                 X,Y
                 
double precision SF(NODEL), SFDL(NDIME,NODEL)
double precision JCBT(NDIME,NDIME),&
                 JCBTI(NDIME,NDIME), DETJ

integer IGAUSS,JGAUSS,I,J,IELEM,P,Q

double precision r_mod,K_W,Omega,dt,Cp,T

double precision DELTAPHIG1,MAGNIPHIG1,&
                 DELTAPHIG2,MAGNIPHIG2,&
                 DELTAPHIE1,MAGNIPHIE1,&
                 DELTAPHIE2,MAGNIPHIE2,&
                 error1,error2


  open(1171,file='error_1_data')
  open(1172,file='error_2_data')
  
  
  NODEF = NODEL*NANGL
  !Get integration points.
  call GAULEG(NGAUSS, XG, NGAUSS, WG, NGAUSS)
  T = IStep*dt
  
  !Initialize global numerator and denominator for error
  DELTAPHIG1=0.0D0
  MAGNIPHIG1=0.0D0
  DELTAPHIG2=0.0D0
  MAGNIPHIG2=0.0D0
  
  DO IELEM = 1,TOTELS
    call GetELNODE(ELNODDAT,NODEL,NORDER)
    ELNDS = ELDAT(IELEM, ofst_ELDAT + ELNODDAT)
    DO P = 1,NODEL
      DO Q = 1,NDIME
        ELCOR(P,Q) = COORD(ELNDS(P),Q)
      endDO
    endDO
    
    call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)

    !Initialize elementary numerator and denominator for error
    DELTAPHIE1=0.0D0
    MAGNIPHIE1=0.0D0
    DELTAPHIE2=0.0D0
    MAGNIPHIE2=0.0D0
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

        
        !Get (X,Y) from (XI,ETA)
        X = 0.0D0
        Y = 0.0D0
        DO I=1,NODEL
          X = X + SF(I)*ELCOR(I,1)
          Y = Y + SF(I)*ELCOR(I,2)
        ENDDO        
                
        !Interpolate numerical
        PhiNum = dcmplx(0.0D0,0.0D0)
        Do I=1,NODEF
          PhiNum = PhiNum + Phi(ELDF(I))*SF(I)
        endDO
        !Calculate exact
        PhiExa = dcmplx(0.0D0,0.0D0)
        r_mod = sqrt(X*X + Y*Y)
        PhiExa = CDEXP(DCMPLX(0.0D0,1)*(K_W*r_mod&
			-Omega*T))
        
        !Calculate local numerators and denominators
        DELTAPHIE1=DELTAPHIE1+CDABS(PHINUM-PHIEXA)**(1)*WTXY
        MAGNIPHIE1=MAGNIPHIE1+CDABS(PHIEXA)**(1)*WTXY
        DELTAPHIE2=DELTAPHIE2+CDABS(PHINUM-PHIEXA)**(2)*WTXY
        MAGNIPHIE2=MAGNIPHIE2+CDABS(PHIEXA)**(2)*WTXY
        
      endDO
    endDO
    !Calculate global numerators and denominators
    DELTAPHIG1=DELTAPHIG1+DELTAPHIE1
    MAGNIPHIG1=MAGNIPHIG1+MAGNIPHIE1
    DELTAPHIG2=DELTAPHIG2+DELTAPHIE2
    MAGNIPHIG2=MAGNIPHIG2+MAGNIPHIE2
  endDO
  error1 = ((DELTAPHIG1)**(1.0D0/1))/((MAGNIPHIG1)**(1.0D0/1))
  error2 = ((DELTAPHIG2)**(1.0D0/2))/((MAGNIPHIG2)**(1.0D0/2))
  write(1171,*) IStep,error1
  write(1172,*) IStep,error2
end subroutine getLn_norm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
