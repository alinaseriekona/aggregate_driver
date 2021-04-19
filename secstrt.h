C       CVS Revision:$Revision: 3.0 $  created $Date: 2017/12/14 $
C       this is section.f file secstrt.h V.1.0 January 2017;
C       it contains pointers for the soot kinetics
C       subroutines' data storage arrays
C
C
C     include file for PRESS section.f, dated: Jan. 2, 2017
C
      INTEGER
     1   NPAH,  NUMNUC, NUMDIM, NumSootReaction,
     2   IC2H2, IO2,    IH2O,   IH,     IH2,    IC2H6,  IN2,
     3   ICH4,  IO,     IOH,    ICO,    ICO2,   IC2H4,  IC6H6,
     4   IHE,   IAR,    MSECTION

C        PARAMETER(NumSootReaction = 10)
      COMMON /SECSTRT/
C
C     Integer constants
C
     1   NPAH, NumSootReaction,
     2   IC2H2, IO2,    IH2O,   IH,     IH2,    IC2H6,  IN2,
     3   ICH4,  IO,     IOH,    ICO,    ICO2,   IC2H4,  IC6H6,
     4   IHE,   IAR
C
C    SOOTSECTION MODEL PARAMETERS
C
      DOUBLE PRECISION
     1   AV,    AMU,    C_MW,   C_MASS, MaxTemp,        BOLTZMANN,
     2   PI,    smallnum,       STFN_BLTZ,      ONETHIRD,
     3   FS,    DFRCT,  AKF,    C_aggl, c_min_mono,     DensityP,
     4   FVOL,  CHI,    GasCon, NumCIncep, DCutOff, SFAC

	COMMON/SECSTRT/ SFAC

C
C     Double precision gas constants
C
      PARAMETER(
     1   AV = 6.022D+23, AMU = 1.D0/AV, C_MW = 12.011,
     2   MaxTemp = 3000., BOLTZMANN = 1.3806488D-16, PI = 3.141592D0)
      PARAMETER(STFN_BLTZ = 5.670373D-5)
      PARAMETER(
     1    smallnum = 1.E-20, ONETHIRD = 0.333333,
     2   GasCon = 1.987D-3, C_MASS = C_MW*AMU)
      PARAMETER(DCutOff = 10.0D0)
C
C     Double precision soot constants
C
      PARAMETER(
     1   FS = 1.9D0, DFRCT = 1.8D0, AKF = 1.37D0, C_aggl = 3.0D0,
     2   c_min_mono = 90000D0, FVOL = 1.43D0, DensityP = 1.9D0,!1.9
     3   CHI = 2.3D+15, NumCIncep = 410.D0,
     4   MSECTION =5, NUMNUC = 3, NUMDIM=NumNuc*(NumNuc+1)/2)
C
C     END include file for section.f
C
