use module_data_read

implicit none
!Definitions
  integer NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,NGAUSS
  integer, allocatable :: ELDAT(:,:)
  integer ofst_ELDAT,ELTYPE
  double precision, allocatable :: COORD(:,:),ELCOR(:,:)
  integer, allocatable :: ELEDGDAT(:,:), ELNODDAT(:), ELNDS(:)
  integer, allocatable :: FACT(:,:,:),EFACT(:,:)  
  integer IStep, NPOINTPLOT, IELEM, I, J,counter
  
  integer TOTNOD
  integer, allocatable ::  GLDF(:,:),MapCoord2Coef(:),ELDF(:)
  integer NANGL
  !definitions for book-keeping 
  character (len=52) char_name
  integer n_record_IStep
  !definitions for DIC,M,S
  integer TOTDOF
  double precision K_W,Omega,PI
  double complex, allocatable :: ELM(:,:),ELK(:,:),&!Elementary Mass & Stiffness
  				 GM(:,:),GS(:,:),&  !Global Mass & Stiffness
  				 InvGM(:,:),&
  				 PTLD_Phi(:),PTLD_Vlcty(:),&
  				 RHS_Phi(:),RHS_Vlcty(:),&
  				 Phi(:),Vlcty(:)
  !definitions for getRHS
  double precision Cp	!Wave speed
  double precision dt
  double complex, allocatable :: RHS(:)
  !definitions for Implicit formulation
  double complex, allocatable :: PrevPhi(:),&
                                 temp_vec(:),InvLHS(:,:),&
                                 LHS(:,:)
  integer n_istep
  
  !definitions for other postprocessing
  double precision error1,error2
  
 
    
  !Get input data  
    call getData(NDIME,TOTELS,NODEL,TOTCOORD,NODEDGE,NORDER,NMARK,&
               ELDAT,ofst_ELDAT,COORD,ELEDGDAT,FACT,ELTYPE,GLDF,&
               TOTNOD,MapCoord2Coef)  
               
  !book-keeping
  char_name = "case_default"
  call system("mkdir "//char_name)	!Used for storing all results
  call system("mkdir PLOTS")		!Used for plots
  open(1002, file='logfile.txt', status="old", position="append", action="write")
               
  !Set-up problem related parameters
  NGAUSS = 2!CEILING((dble(NORDER) + 1)/2) + 5	!for p-FEM (the +1 is not 
  						!enough for p>2, I guess 
  						!for the wave integration, 
  						!because x,y interpolation
  						!works fine with lower ip)
  write(*,*)'NGAUSS = ',NGAUSS
  write(1002,*)'NGAUSS = ',NGAUSS
  NANGL = 1
  write(1002,*)'NANGL = ',NANGL,' is 1 for FEM'
  TOTDOF = TOTNOD*NANGL
  write(1002,*)'TOTDOF = ',TOTDOF
  PI = 4.0D0*DATAN(1.0D0)
  K_W = 1.0D0*PI
  Omega = 1.0D0
  Cp = Omega/K_W
  dt = 0.001
  n_istep = 1000

  
  !allocate variables
    allocate(ELNODDAT(NODEL),ELNDS(NODEL),ELCOR(NODEL,NDIME))
    allocate(EFACT(NMARK,ELTYPE))
    allocate(ELDF(NODEL*NANGL))
    !allocate for DIC,M,S
    allocate (ELM(NODEL,NODEL),ELK(NODEL,NODEL))
    allocate (PTLD_Phi(NODEL),PTLD_Vlcty(NODEL))
    allocate (GM(TOTDOF,TOTDOF),InvGM(TOTDOF,TOTDOF))
    allocate (GS(TOTDOF,TOTDOF))
    allocate (RHS_Phi(TOTDOF),RHS_Vlcty(TOTDOF))
    allocate (Phi(TOTDOF),Vlcty(TOTDOF))
    !allocate for getRHS
    allocate (RHS(TOTDOF))
    !allocate for Implicit formulation
    allocate (PrevPhi(TOTDOF))
    allocate (temp_vec(TOTDOF))
    allocate (LHS(TOTDOF,TOTDOF))
    allocate (InvLHS(TOTDOF,TOTDOF))
  
  !~~~~~~~~~~~~~~~~~~~~DIC,M,S~~~~~~~~~~~~~~~~~~~~!Starts
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'Processing DIC,M,S'
  write(*,*)'~~~~~~~~~~~~~~~~~~~'  
  !Initialize matrices and vectors
  GM = dcmplx(0.0D0,0.0D0)
  GS = dcmplx(0.0D0,0.0D0)
  RHS_Phi = dcmplx(0.0D0,0.0D0)
  RHS_Vlcty = dcmplx(0.0D0,0.0D0)
  Phi = dcmplx(0.0D0,0.0D0)
  Vlcty = dcmplx(0.0D0,0.0D0)
  !Assemble
  Do IELEM = 1,TOTELS
    !get element co-ordinates
    call GetELNODE(ELNODDAT,NODEL,NORDER)
    ELNDS = ELDAT(IELEM, ofst_ELDAT + ELNODDAT)
    DO I = 1,NODEL
      DO J = 1,NDIME
        ELCOR(I,J) = COORD(ELNDS(I),J)
      endDO
    endDO
        
    !Build local system    
    call getEl_DIC_M_S(IELEM,ELCOR,NODEL,NDIME,&
		   NGAUSS,NORDER,K_W,Omega,&
		   ELM,ELK,PTLD_Phi,PTLD_Vlcty)
		       
    !Build global system
      !get indices that map from ELememnt Degree of Freedom {ELDF} to GLobal Degree of Freedom {GLDF table}
      call getELDF(IELEM,GLDF,TOTELS,NODEL,ELNODDAT,NORDER,ELDF,NANGL)
      !transfer ELementary matrices/vectors to Global matrices/vectors
      DO I = 1,NODEL
        DO J = 1,NODEL
          GM(ELDF(I),ELDF(J)) = GM(ELDF(I),ELDF(J)) + ELM(I,J)
          GS(ELDF(I),ELDF(J)) = GS(ELDF(I),ELDF(J)) + ELK(I,J)
        endDO
        RHS_Phi(ELDF(I)) = RHS_Phi(ELDF(I)) + PTLD_Phi(I)
        RHS_Vlcty(ELDF(I)) = RHS_Vlcty(ELDF(I)) + PTLD_Vlcty(I)
      endDO    
  endDO  
  !Solve for Decomposition of Phi0,Vlcty0
    !Invert Mass
    call inv(GM,TOTDOF,InvGM,TOTDOF)
    !Multiply inverted Mass with RHS to obtain Phi0,Vlcty0
    CALL MATMULCPLX(InvGM,TOTDOF,TOTDOF,RHS_Phi,TOTDOF,1,Phi,TOTDOF,1,TOTDOF,1,TOTDOF)
    CALL MATMULCPLX(InvGM,TOTDOF,TOTDOF,RHS_Vlcty,TOTDOF,1,Vlcty,TOTDOF,1,TOTDOF,1,TOTDOF)
    !get Phi(-1) from Phi0,Vlcty0: Phi(-1) = Phi0 - dt*Vlcty0
    PrevPhi = Phi - (dt*Vlcty)	!This is useful for implicit formulation

  !~~~~~~~~~~~~~~~~~~~~DIC,M,S~~~~~~~~~~~~~~~~~~~~!Stops
  write(*,*)'~~~~~~~~~~~~~~~~~~~'
  write(*,*)'LHS for implicit-wave'
  write(*,*)'~~~~~~~~~~~~~~~~~~~'  
  
  !Update Stiffness for wave problem
  GS = (Cp*Cp)*GS
  
  !Compute LHS for implicit (Backward-Euler)
    !LHS = [M] + dt*dt*[S]
    LHS = GM + (dt*dt*GS)
    !Invert LHS
    call inv(LHS,TOTDOF,InvLHS,TOTDOF)
  
  
  !~~~~~~~~~~~~~~~~~~~~TimeLoop~~~~~~~~~~~~~~~~~~~~!Starts
  Do Istep = 1,n_istep
    write(*,*)'~~~~~~~~~~~~~~~~~~~'
    write(*,*)'Time-step ', IStep
    write(*,*)'~~~~~~~~~~~~~~~~~~~'  
    !getRHS
      call getRHS(TOTELS,NODEL,NORDER,NGAUSS,ofst_ELDAT,ELDAT,&
                 NDIME,COORD,TOTCOORD,NMARK,ELTYPE,FACT,GLDF,NANGL,&
                 IStep,dt,Cp,K_W,Omega,RHS,TOTDOF)
                 
    !Compute RHS for implicit (Backward-Euler)
      !RHS = [M](2Un-1 - Un-2) + dt*dt*<r,w>
      temp_vec = dcmplx(0.0D0,0.0D0)
      call MATMULCPLX(GM,TOTDOF,TOTDOF,&
      		    ((dcmplx(2.0D0,0.0D0)*Phi) - PrevPhi),TOTDOF,1,&
                    temp_vec,TOTDOF,1,TOTDOF,1,TOTDOF)
      RHS = temp_vec + ((dt*dt)*RHS)
    
    !Solver for Phi (and Vlcty)    
      !Save history of Phi(t-1)
      PrevPhi = Phi
      !compute Phi(t)
      call MATMULCPLX(InvLHS,TOTDOF,TOTDOF,RHS,TOTDOF,1,&
    		   Phi,TOTDOF,1,TOTDOF,1,TOTDOF)
  
  
    !Post Processing      
      !Errors
      call getLn_norm(IStep,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                   TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,GLDF,&
                   ELCOR,NGAUSS,&
                   K_W,Omega,dt,Cp,error1,error2)
      write(*,*)'L2 error % =',error2*100
      write(*,*)'L1 error % =',error1*100
      !Plots      
      NPOINTPLOT = 10	!# of integration points for plots
      n_record_IStep = CEILING(dble(n_istep)/20.0D0)
      if(((IStep .EQ. 1).OR.(mod(IStep,n_record_IStep).EQ. 0)).OR.(IStep .EQ. n_istep))Then
        call plotPara(IStep,NPOINTPLOT,TOTELS,ELDAT,COORD,ofst_ELDAT,NODEL,&
                   TOTCOORD,NDIME,NORDER,NANGL,TOTDOF,Phi,Vlcty,GLDF,&
                   K_W,Omega,dt)
      endif                   
      !Exports
      !open(1020, file = 'realMASS')
      !open(1021, file = 'imagMASS')
      !Do I = 1,TOTDOF
      !  DO J = 1,TOTDOF
      !    write(1020,1102,advance='no')dreal(GM(I,J))
      !	write(1021,1102,advance='no')dimag(GM(I,J))
      !	1102 format(f16.7)
      !  endDO
      !  write(1020,*)
      !  write(1021,*)
      !endDO
      !close(1020)
      !close(1021)  
  endDO !{Time}
  !~~~~~~~~~~~~~~~~~~~~TimeLoop~~~~~~~~~~~~~~~~~~~~!Stops
  
  !book-keeping	
  call system("cp dat "//char_name)
  call system("cp -r PLOTS "//char_name)
  call system("cp error_1_data "//char_name)
  call system("cp error_2_data "//char_name)
  call system("cp logfile.txt "//char_name)
  call system("rm dat")
  call system("rm -r PLOTS")
  call system("rm error_1_data")
  call system("rm error_2_data")
  call system("rm logfile.txt")
                   
END
