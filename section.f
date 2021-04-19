      SUBROUTINE SootSecRate(Pressure,T, RHO,SootSec,
     1                      XS,P1,SOURCETERM) ! #HTE
        !This subroutine takes the nucleation and surface HACA and
        !condensation growth rates and determines the source terms for
        !soot mass fraction and particles for the 2-eq model 
        !for a given CV, defined by the current I and J global
        !variables. These values are stored in SOURCETERM    
        IMPLICIT NONE
        INCLUDE 'secstrt.h'
C
        LOGICAL CUT_OFF, SINTER, LCOAG

C N        DOUBLE PRECISION  tao_i_i, tao_i_ip1
        DOUBLE PRECISION XS(MSECTION)
	DOUBLE PRECISION P1(MSECTION)
C N        COMMON / SOOT / XS
C N        COMMON /SOOT_FRAG/ tao_i_i, tao_i_ip1
        DOUBLE PRECISION Pressure, T, RHO, HAB
C N        INTEGER  NCF(*), IPAH(*)
        DOUBLE PRECISION SootSec(2*MSection),
     1                   SootGasRates(MSection,NumSootReaction),
     2                   NDim(NUMDIM),SootNucRates(NumNuc+NUMDIM)
        DOUBLE PRECISION SOURCETERM(2*MSection),SOURCEDIM(NUMDIM),
     1                   SOURCEH(MSECTION), NHSec(MSECTION) ! #HTE
          !local variable declarations
        DOUBLE PRECISION  COAGRATE(2*MSection), NucRate(2*MSection),
     +                    FragRate(2*MSection), SurfRate(2*MSection), 
     +                    SintRate(2*MSection), ObliRate(2*MSection), 
     +                    PNX(MSection), S_frag(MSection), 
     +                    DPRIME(MSECTION), CSINT(MSECTION),
     +                    CS1(MSection,MSection,MSection), 
     +                    CNS1(MSection,MSection,MSection), 
     +                    CS2(MSection,MSection),
     +                    CNS2(MSection,MSection),
     +                    CCOAG1(MSection,MSection,3*MSection),
     +                    CCOAG2(MSection,3*MSection) ! #HTE
C       --------    #HTE
        DOUBLE PRECISION PHX(MSection), PAHH(NumNuc),
     +                   HCOAGRATE(MSection), HNucRate,
     +                   HFragRate(MSection), HSurfRate(MSection),
     +                   HSintRate(MSection), HObliRate(MSection),
     1                   CHS1(MSection,MSection,MSection),
     +                   CHS2(MSection,MSection), HCarbRate(MSection),
     +                   Htemp, HOUTO2, HOUTOH, HINPAH, HINC2H2
        DOUBLE PRECISION ACarb, ECarb
        PARAMETER (ACarb = 1.85D9, ECarb = 58.D0)
C 
        DOUBLE PRECISION CoagEff
        PARAMETER (CoagEff = 1.0D0)!2.5D-1 
        DOUBLE PRECISION TotalNucRate, BetaCoag, AINPAH, AINC2H2, AINO2,
     1                   AINOH, AOUTPAH, AOUTC2H2, AOUTO2, AOUTOH, temp,
     2                   total_rate_for_frag, specific_ox_rate, Area,
     3                   AINPAHR, AOUTPAHR, NucDRate(NUMDIM)
        INTEGER L, R, Sec, INIT
        DOUBLE PRECISION W, WL, WU, UTOT, DUMMY
        DOUBLE PRECISION A_frag
        DOUBLE PRECISION UPR,VPR,BETAC,TcoagEff 
        DOUBLE PRECISION PAH_C(NPAH)
        DOUBLE PRECISION, External:: ppMassEffec
        INTEGER NP
	!DOUBLE PRECISION, DIMENSION (:), POINTER :: ptr
C
C
        !done local variable declarations


        !initialize all matrices
        DO INIT=1, 2*MSection
            COAGRATE(INIT) = 0.d0
            SOURCETERM(INIT) = 0.0D0
        ENDDO
	!allocate(ptr(2*MSECTION))
        !determine the number of primary particles per aggregate, and
        !store in PNX
        DO L=1, MSection
            IF(SootSec(L).LE.smallnum .OR.
     1         SootSec(L+MSection).LE.smallnum) THEN
                PNX(L) = 1.D0
                !WRITE(*,*) 'IF',L ,PNX(L)
            ELSE
                PNX(L) = SootSec(L+MSection)/(SootSec(L))
                !WRITE(*,*) 'ELSE', L, PNX(L)
            ENDIF
        ENDDO
		P1=PNX
        !done determining the number of primary particles per aggregate.              
        !*****************************************************************          
            
        !***********************************************************************
        !second, determine coagulation rates  
            !coag. of sections r, s <= L, mass->L
        
        IF(CoagEff .EQ. 0)THEN
                GOTO 9000
        ELSEIF(LCOAG)THEN
        ! Ali: Skip the ELSEIF
        !GO TO 2000


            DO L=1, 2*MSection
            DO Sec=1, MSection
            DO R=1, MSection
                CCOAG1(R,Sec,L) = 0.D0
            ENDDO
                CCOAG2(Sec,L) = 0.D0
            ENDDO
            ENDDO
C
            DO L=1,MSection
                W=EXP(XS(L))
                IF (L.EQ.1) THEN
                    WL=W
                ELSE
                    WL=EXP(XS(L-1))
                ENDIF
                IF (L.EQ.MSection) THEN
                    WU=W*FS
                ELSE
                    WU=EXP(XS(L+1))
                ENDIF
                IF(PNX(L).NE.0.0)THEN
                       DPRIME(L) = (6.0*W/DensityP/PNX(L)/PI)**(1./3.)*1.E+7
                ELSE
                       
                       DPRIME(L) = (6.0*W/DensityP/PI)**(1./3.)*1.E+7
                ENDIF
            DO R=1,L
            DO Sec=R,L
                UTOT = EXP(XS(R))+ EXP(XS(Sec))

                UPR=ppMassEffec(SootSec(R),SootSec(R+Msection),
     1                              EXP(XS(R)))
                VPR=ppMassEffec(SootSec(Sec),SootSec(Sec+Msection),
     1                              EXP(XS(Sec)))

                IF(UTOT .GT. WL .AND. UTOT .LT. W) THEN
!                    CALL BETA(Pressure, T, DUMMY, PNX(R), R, 
!     1                        DUMMY,PNX(Sec),Sec,BETACOAG)
                 BETACOAG = BETAC(XS(R),XS(SEC),UPR,VPR,T,Pressure) !BETA FUNCTION

                    CS1(R,Sec,L)=(WL-UTOT)/(WL-W)*BETACOAG !eta

                    IF(R.EQ.Sec) THEN
                        CS1(R,Sec,L)=CS1(R,Sec,L)/2.d0 !1-delta/2
                    ENDIF
                ELSEIF(UTOT .GE. W .AND. UTOT .LT. WU) THEN   
!                    CALL BETA(Pressure, T, DUMMY, PNX(R), R,
!     1                        DUMMY, PNX(Sec), Sec, BETACOAG)
                 BETACOAG = BETAC(XS(R),XS(SEC),UPR,VPR,T,Pressure) !beta

                    CS1(R,Sec,L)=(WU-UTOT)/(WU-W)*BETACOAG !eta

                    IF(R.EQ.Sec) THEN
                        CS1(R,Sec,L)=CS1(R,Sec,L)/2.d0 !1-delta/2
                    ENDIF
                ELSE
                        CS1(R,Sec,L)=0. !before 0.
                ENDIF
                !IF(CUT_OFF .AND. DPRIME(L) .LE. DCutOff)THEN
                !     CNS1(R,Sec,L)=CS1(R,Sec,L) !pp from agg.
                !ELSE
                     CNS1(R,Sec,L)=CS1(R,Sec,L)*(PNX(Sec)+PNX(R))*W/UTOT !pp
                !ENDIF
C
C
                CCOAG1(R,Sec,L) = CS1(R,Sec,L) 
                CCOAG1(R,Sec,L+MSection) = CNS1(R,Sec,L)
            ENDDO
            ENDDO
            ENDDO 
            ! Ali: running a block test: test passed
          
            DO L=1,MSection !second term in the right hand side
                DO Sec=L,MSection          
                UPR=ppMassEffec(SootSec(L),SootSec(L+Msection),
     1                              EXP(XS(L)))
                VPR=ppMassEffec(SootSec(Sec),SootSec(Sec+Msection),
     1                              EXP(XS(Sec)))
!                    CALL BETA(Pressure, T, DUMMY, PNX(L), L, 
!     1                        DUMMY, PNX(Sec), Sec, BETACOAG)
                    BETACOAG = BETAC(XS(L),XS(SEC),UPR,VPR,T,Pressure)
                    CS2(Sec,L)=0.d0-BETACOAG !negative beta
                    CS2(L,Sec)=CS2(Sec,L)
                    CNS2(Sec,L)=CS2(Sec,L)*PNX(L)
                    CNS2(L,Sec)=CS2(L,Sec)*PNX(Sec)
C
                ENDDO
            ENDDO
C           
            ! Ali: running a block test: test failed
            !***************************************
            DO L=1,MSection
                DO Sec=1,MSection
                CCOAG2(Sec,L) = CS2(Sec,L)
                CCOAG2(Sec,L+MSection) = CNS2(Sec,L)

                ENDDO
            ENDDO
            !**************************************
2000    CONTINUE
        !WRITE(*,*) 'ELSEIF SKIPPED'
        ELSE
C
        ! Ali: skip ELSE  
            DO L=1,MSection
            DO R=1,L
            DO Sec=R,L
                CS1(R,Sec,L) = CCOAG1(R,Sec,L)
                CNS1(R,Sec,L) = CCOAG1(R,Sec,L+MSection)
            ENDDO
            ENDDO
            ENDDO
C
            DO L=1,MSection
                DO Sec=1,MSection
                  CS2(Sec,L)  = CCOAG2(Sec,L)
                  CNS2(Sec,L) = CCOAG2(Sec,L+MSection)
                ENDDO
            ENDDO
3000    CONTINUE
        !WRITE(*,*) 'ELSE SKIPPED'
        END IF
C
        !Ali: skip the main loop
        !WRITE(*,*) 'Mail loop skipped!'
        !GO TO 8000
            DO L=1, MSection
            !Ali: skip the inner nested loop
            !WRITE(*,*) 'Skip the inner nested loop'
            !GO TO 4000
!*----------------------------------------------------------------
                DO R=1, L
                DO Sec=R, L
                    CoagRate(L)=CoagRate(L)+CS1(R,Sec,L)*
     1                          SootSec(R)*SootSec(Sec)*RHO**2.d0
                    CoagRate(MSection+L)=CoagRate(MSection+L)+
     1                                   CNS1(R,Sec,L)*SootSec(R)*
     2                                   SootSec(Sec)*RHO**2.d0
                        
               ! WRITE(*,*) 'R,Sec,CS1(R,Sec,L), SS(R),SS(Sec),RHO'
       ! WRITE(*,*) R,Sec, CS1(R,Sec,L), SootSec(R), SootSec(Sec),RHO
                ENDDO
                ENDDO        
!*---------------------------------------------------------------
4000    CONTINUE


            !Ali: skip the inter loop
            !WRITE(*,*) 'Skip the inner loop'
            !GO TO 5000
                DO Sec=1,MSection
                    CoagRate(L)=CoagRate(L)+CS2(Sec,L)*SootSec(Sec)*
     1                          SootSec(L)*RHO**2.d0
                    CoagRate(MSection+L)=CoagRate(MSection+L)+
     1                                   CNS2(Sec,L)*SootSec(Sec)*
     2                                   SootSec(L)*RHO**2.d0
                
                ENDDO
5000     CONTINUE
                 
            ENDDO        
8000    CONTINUE
        

9000    CONTINUE
	! done determining coagulation rates
	!***********************************************
        DO L=1,2*MSection
            SOURCETERM(L) =  CoagEff*CoagRate(L) 
            if(CoagRate(L).GT.0.0.OR.CoagRate(L).LT.0.0)then
                 !WRITE(*,*) 'error',L, CoagRate(L)
!*----------------------------------------------------------------
                DO R=1, L
                DO Sec=R, L
                !WRITE(*,*) 'R,Sec,CS1(R,Sec,L), SS(R),SS(Sec),RHO'
        !WRITE(*,*) R,Sec, CS1(R,Sec,L), SootSec(R), SootSec(Sec),RHO
                ENDDO
                ENDDO
!*---------------------------------------------------------------

            endif     
        ENDDO    
C
      end subroutine SootSecRate
C**************************************************************************************************************************************************************
c ----------------------------------------------------------------------CF
C
	DOUBLE PRECISION FUNCTION ppMassEffec(N,NP,U)
	 IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
	 DOUBLE PRECISION N, NP
	IF (N.LE.1E-20.OR.NP.LE.1.E-20) THEN
	 ppMassEffec=U
	ELSE
	 ppMassEffec=U*N/NP
	ENDIF
	RETURN
	END


c **************************************************************************************************************************************************************
      double precision function viscosity(T)
      !This function calculates the viscosity of air as a function of
      !temperature (K) (from AK, myf 1997)

          !local variable declarations
            implicit none
            double precision T  
            !done local variable declarations
        
            viscosity = 1.458D-05 * (T**1.5D0) / (T + 110.4D0)
            return
      end
C**************************************************************************************************************************************************************
C**************************************************************************************************************************************************************
      double precision function freePath(P, T)
      !This function calculates the mean free path for air in cm (from
      !AK, myf 1997)
    
          !local variable declaration
            implicit none
            double precision P,T
            !done local variable declarations
        
            freePath = 2.3701D-02 * T / P 
        
            return
      end
C**************************************************************************************************************************************************************
      DOUBLE PRECISION FUNCTION BETAC(X,Y,UPRIME,VPRIME,T,P)
      !Fuchs beta used for aggdriv.f, T = 1830, P=1, 12ms, 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER(I-N)
      DOUBLE PRECISION ONE_THIRD, UPRIME, VPRIME,TGAS,FREEMP,VISCOS,
     1    Viscosity, freePath, DI, DJ, DFI, DFJ, GIJ, GI, GJ,
     2    CIJ, CBARI, CBARJ, LMDAI, LMDAJ, BETACO,BETACM,KN,DF
     3    CCFI,CCFJ,FPAIR,LP,KNUD,NX,NY,DFM,KFM,RMI,RMJ,RGI,RGJ
     4    DMI,DMJ,MWair,RCONST,AIRDENSITY,VISCOSAIR,NP,DPI,DPJ
     5    DGI,DGJ
     

      INCLUDE 'secstrt.h'
      ONE_THIRD = 1.0D0/3.0d0
      VISCOS = viscosity(T)*1.0D-1 !kg/ms!1.846D-5
      FPAIR=6.7D-8!1.4D-7freePath(T,P)*1.0D-2 !mean free path for air in m
      TGAS = T
      PDENSITY = densityP*1.0D3 !KG/M^3

      !!Ali Equations
      !MWair = 0.0289647 !kg/mol
      !RCONST = 8.205736D-5 !M^3*ATM/K/MOL
      !AIRDENSITY = MWair*(1/(RCONST*T))
      !VISCOSAIR=1.425D-6*(T**(0.5039)/(1+(108.3/T))) !kg/ms
      !FPAIR=(VISCOSAIR/AIRDENSITY)*
      !1    SQRT((PI*MWair)/(2*BOLTZMANN*1.0D-7*T))
      !print*,FPAIR

      !Compute the geometric properties of the cluster (AGGREGATE)
      U = EXP(X)*1.0D-3 !kg
      V = EXP(Y)*1.0D-3
      UPRIME = UPRIME*1.0D-3
      VPRIME = VPRIME*1.0D-3
      !print*, 'masses of clusters/particles',U,V,UPRIME,VPRIME
      
      !RADIUS OF THE PRIMARY PARTICLES
      RXP=((3.*UPRIME/4./PI/PDENSITY)**ONE_THIRD) !m
      RYP=((3.*VPRIME/4./PI/PDENSITY)**ONE_THIRD)
      DPI=RXP*2 !m
      DPJ=RYP*2
      !!RXP=6.7D-5
      !!DO WHILE (RXP.GE.4.47d-10)
            
      !!CONSTANTS
       !IF (KNUD.LE.1.0D1) THEN 
        !KN=1.3D0
        !DF=1.78D0
        !DFM=2.27D0
        !KFM=1.11D0
       !ELSE IF (KNUD.GE.1.0D1) THEN 
        KFM=1.11D0
        KN=1.4D0
        DF=1.91D0
        DFM=2.15D0
       !END IF 
       !print*, FPAIR 
       !print*, KNUD 
       !!RXP=FPAIR/KNUD ! PP radius in m
       KNUD=FPAIR/RXP 
       !RYP=RXP
       !print*, RXP

       !NUMBER OF PARTICLES PER AGGREGATE
       IF (U.GE.UPRIME) THEN
        NX = U/UPRIME !g/g
       ELSE 
        NX=1
       ENDIF

       IF (V.GE.VPRIME) THEN
        NY=V/VPRIME
       ELSE
        NY=1
       ENDIF
       !PRINT*, NY,NX  
       !UPRIME=(RXP**3.0)*(4.0*PI*PDENSITY/3.0)
       !VPRIME=(RYP**3.0)*(4.0*PI*PDENSITY/3.0)
       !PRINT*, UPRIME
       !NP=1
       
       !!!This whole section is from the paper
       !Computing the radii of gyration
       !RGI=RXP*(NX/KN)**(1/DF)
       !RGJ=RYP*(NY/KN)**(1/DF)
       !DI = 2.0*RGI  !Collision diameter in m
       !DJ = 2.0*RGJ

       !Computing the mobility radii
       !RMI=RXP*(NX/KFM)**(1/DFM)
       !RMJ=RYP*(NY/KFM)**(1/DFM)    
       !DMI = 2.0*RMI ! MOBILITY diameter in m
       !DMJ = 2.0*RMJ
       !full coalescense
       !DI = DMI !COLLISION = MOBILITY diameter in m
       !DJ = DMJ
       !!!!!!!!!!!!!

       !MOBILITY IN M
       DMI=DPI*(NX**0.45)!m
       DMJ=DPJ*(NY**0.45)
       !GYRATION IN M
       DGI=DMI/(NX**(-0.2)+0.4)
       DGJ=DMJ/(NY**(-0.2)+0.4)
       !COLLISION DIAMETER IN M
       IF (DMI.GE.DGI) THEN
        DI=DMI
       ELSE
        DI=DGI
       ENDIF
       IF (DMJ.GE.DGJ) THEN
        DJ=DMJ
       ELSE
        DJ=DGJ  
       ENDIF    

       !1 COMPUTING THE MEAN PARTICLE VELOCITY
       CBARI = SQRT((8.0*BOLTZMANN*1.0D-7*TGAS)/(PI*U)) !M/S !UPRIME
       CBARJ = SQRT((8.0*BOLTZMANN*1.0D-7*TGAS)/(PI*V))
    
       !2 COMPUTING THE AVERAGE MEAN PARTICLE VELOCITY
       CIJ = SQRT((CBARI**2.0+CBARJ**2.0))
       !PRINT*, CIJ
    
       !3 COMPUTING THE CUNNINGHAM CORRECTION FACTOR
       !from wikipedia
       !CCFI=1+(2*FPAIR/DI)*(1.257+0.4*EXP((-0.55*DI)/FPAIR))
       !CCFJ=1+(2*FPAIR/DJ)*(1.257+0.4*EXP((-0.55*DJ)/FPAIR))
       !Ali 
       CCFI=1+(2*FPAIR/DI)*(1.21+0.4*EXP((-0.78*DI)/FPAIR)) !COLLISION DIAMETER?
       CCFJ=1+(2*FPAIR/DJ)*(1.21+0.4*EXP((-0.78*DJ)/FPAIR))
     
       !4 COMPUTING PARTICLE DIFFUSIVITY
       DFI = (BOLTZMANN*1.0D-7*TGAS*CCFI)/(3*PI*VISCOS*DMI) !m^2/s
       DFJ = (BOLTZMANN*1.0D-7*TGAS*CCFJ)/(3*PI*VISCOS*DMJ)

       !5 CALCULATING MEAN FREE PATH OF PARTICLES
       LMDAI = (8.0*DFI)/(PI*CBARI) !M
       LMDAJ = (8.0*DFJ)/(PI*CBARJ)

       !6 COMPUTING MEAN DISTANCE IN M
       GI = (1.0/(3.0*DI*LMDAI))*((DI+LMDAI)**3.0- 
     1 	(DI**2.0+LMDAI**2.0)**1.5)-DI
       GJ = (1.0/(3.0*DJ*LMDAJ))*((DJ+LMDAJ)**3.0- 
     1 		(DJ**2.0+LMDAJ**2.0)**1.5)-DJ

       !7 COMPUTING THE AVERAGE MEAN DISTANCE M
       GIJ = SQRT(GI**2.0+GJ**2.0)
    
       !8 COMPUTING BETA COLLISION FREQUENCY FOR CONTINUUM REGIME
       BETACO = 2.0*PI*(DI+DJ)*(DFI+DFJ) !M^3/S
       !print*, 'continuum beta', BETACO

       !COMPUTING THE KNUDSEN NUMBER 
       !KNDSN=FPAIR/RXP

       !9 COMPUTING FUCHS COLLISION FREQUENCY
       BETACM = BETACO*(((DI+DJ)/(DI+DJ+2*GIJ))+
     1       ((8*(DFI+DFJ))/(CIJ*(DI+DJ))))**(-1)
       BETAC=BETACM*1.0D6 !making it cm^3/s
      ! write(*,"('frequency function',',',ES10.2E2,',',ES10.2E2,
      !1      ',',ES10.2E2)") 
      !2      KNUD,BETAC,RXP

        !RXP=RXP-1.0D-10
        !KNUD = KNUD+1.0D-1

      !END DO 
      !print*, 'END'
      RETURN
      END
