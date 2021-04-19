	program aggdrive
	
	implicit none
	include 'secstrt.h'
	integer :: i,j,N,K,L
	double precision :: PNX(MSECTION), DIAM(MSECTION)
	DOUBLE PRECISION :: MM(MSECTION), MMA(MSECTION), W(MSECTION) !mass matrix will store the different mass bins
	double precision :: SSEC(2*MSECTION) !aggregate, primary particle
	DOUBLE PRECISION :: SUMAGG = 0.0D0, SUMPP = 0.0D0, SUMDIAM
	DOUBLE PRECISION :: DIAM_MEAN
	double precision :: SST(2*MSECTION) !source terms
	DOUBLE PRECISION :: DM = 1.0D-7, t
	!double precision, parameter :: tres=12.0D-3,tstep = 1.0D-6 !residence time and time step red -5, -6
	double precision, parameter :: tres=2.0D-6,tstep = 1.0D-6 !residence time and time step red -5, -6
	DOUBLE PRECISION, PARAMETER :: P=1, TEMP=1830 !PRESSURE IN ATM AND TEMPERATURE IN K
C        density=1, residence time=10,timestep=1, sourceterm=1
c	
         print*, 'MSECTION', MSECTION
c	validation using Temp =1830, tres 12.0D-3, tstep=1.0D-6, NumCincep=410, DFRCT


c 1   FS = 1.9D0, DFRCT = 1.8D0, AKF = 1.37D0, C_aggl = 3.0D0,
c     2   c_min_mono = 90000D0, FVOL = 1.43D0, DensityP = 1.9D0,!1.9
c     3   CHI = 2.3D+15, NumCIncep = 410.D0,
c	 4   MSECTION =1, NUMNUC = 3, NUMDIM=NumNuc*(NumNuc+1)/2)
	 


	!PRINT*, 'AGGREGATE DRIVER PROGRAM'
	!PRINT*, 'RESIDENCE TIME', TRES
	!PRINT*, 'TEMPERATURE', TEMP
	!PRINT*, 'TIME STEP', TSTEP
	!Initializing the arrays, - each section (I) has 10 aggregates and 10 PP's
	do I = 1,MSECTION
	 SSEC(I)= 0.0D0 !2.6234E12 !2.86E12 9E9
	 SSEC(I+MSECTION)= 0.0D0 !2.6234E12
	 if (I>MSECTION-24) then
	  SSEC(I)=0.0D0 !1.0D-20 
	  SSEC(I+MSECTION)= 0.0D0 !1.0D-20
	 endif
!	 write(*,"('IAGG',I3,',',ES13.2E3,',',2X,'IPP',I3,',',
!	1  ES13.2E3,',')") I,SSEC(I),I,SSEC(I+MSECTION)

	 !MM(I)=0.0
	 MMA(I)=0.0
	 SST(I)=0.0
	end do
	t=0.0
	DO K=1, MSECTION
	 SUMAGG=SSEC(K)+SUMAGG
	 SUMPP=SSEC(K+MSECTION)+SUMPP
 	ENDDO
!	write(*,"('Time =',',',1X,ES12.2E3,','2X,'SUM AGG',',',2X,
!	1  ES13.2E3,','2X,'SUM PP',',',2X,ES13.2E3, ',')") t, SUMAGG, SUMPP
	PRINT*, ' '

	!print*,"there are",MSECTION,"sections in this program"
        !print*,"the initial aggs:", SSEC(1:MSECTION)
	!print*, 'the initial sum of primary', SUMAGG
	!print*,"the initial primary:", SSEC(MSECTION+1:2*MSECTION)
	!print*, 'the initial sum of primary', SUMPP
C	!Creating the mass bins
        !Mass bin using the array from INITSOOT subroutine
	MMA(1)=LOG(C_MASS*NumCIncep)
        W(1)=EXP(MMA(1))
		DIAM(1)=(((6*W(1))/(DensityP*PI))**(1./3.))*1.0D7
	!MM(1)=DensityP*(PI*(4.0/3.0)*(DIAM/2)**3)
	DO I=2,MSECTION
	 !MM(I)=MM(I-1)*FS
	 MMA(I)=MMA(I-1)+LOG(FS)
	 W(I)=EXP(MMA(I))
	 DIAM(I)=(((6*W(I))/(DensityP*PI))**(1./3.))*1.0D7
	 PRINT*, 'DIAMETERS', I, DIAM(I)

        END DO
		!PRINT*, 'the mass matrix is: ', MM
	!PRINT*,'the log mass matrix is', MMA

	!calling sootsecrate subroutine to get source terms SST
	CALL SOOTSECRATE(P,TEMP,1.0D0,SSEC,MMA,PNX,SST)
	!CALL SOOTSECRATE(P,TEMP,DensityP,SSEC,MM,SST)
c
	!PRINT*, 'PRIMARY PARTICLES/AGGREGATE:',PNX
	!PRINT*,'the source terms are: ', SST
C
	!Solving equation for each section through residence time
	!do J=1,MSECTION
	!print*,"section",J
	  write(*,"('Time =',',',1X,ES12.2E3,',',3X,'SUM AGG',',', 
	1 ES12.2E3,',',3X,'SUM PP',',',ES12.2E3,',',
	2 2X,'MEAN DIAMETER',',',ES12.2E3)") t,SUMAGG,SUMPP,DIAM_MEAN
	 
         t=tstep
	 do
	 SUMAGG=0.0D0
	 SUMPP=0.0D0
	 SUMDIAM=0.0D0
	 DIAM_MEAN=0.0D0
	  if (t>tres) exit !start at t=0 and end at residence time
	  !write(*,"('Time =',',', 1X,ES12.2E3)") t  
	  do J=1,MSECTION !obtain the next iteration in each bin for t
           !print*, 'Section', J, SST(J)
           SSEC(J)=SSEC(J)+tstep*SST(J)
	   SSEC(J+MSECTION)=SSEC(J+MSECTION)+tstep*SST(J+MSECTION)
	   if (SSEC(J)<=0.0D0) then !if the particles are 0 or lower, the bins are empty
	    SSEC(J)=1.0D-20 
	   endif
	   if (SSEC(J+MSECTION)<=0.0D0) then
	    SSEC(J+MSECTION)=1.0D-20 
	   endif
C	   write(*,"('AGG',I3,2X,ES13.2E3,',',3X,'PP',I3,2X,
C	1  ES13.2E3,',',3X,'MASS',2X,
C	2   ES13.2E3)") J,SSEC(J),J,SSEC(J+MSECTION),W(J)
C	
c	   write(*,"('AGG',I3,',',ES13.2E3,',',3X,'sstAGG',',',ES13.2E3,
c	1   ',',3X,'PP',I3,',',ES13.2E3,',',3X,'sstPP',',',ES13.2E3,',',
c	2   3X,'MASS',I3,',',ES13.2E3,2X,',','DIAM',','ES13.2E3)") 
c	3   J,SSEC(J),SST(J),J,SSEC(J+MSECTION),
c	4   SST(J+MSECTION),J,W(J),DIAM(J)
C
	   SUMAGG=SUMAGG+SSEC(J)
	   SUMPP=SUMPP+SSEC(J+MSECTION)
	   SUMDIAM=SUMDIAM+(SSEC(J)*DIAM(J))
 	 end do
	 DIAM_MEAN = SUMDIAM/SUMAGG
	 CALL SOOTSECRATE(P,TEMP,1.0D0,SSEC,MMA,PNX,SST)
	  write(*,"('Time =',',',1X,ES12.2E3,',',3X,'SUM AGG',',', 
	1 ES12.2E3,',',3X,'SUM PP',',',ES12.2E3,',',
	2 2X,'MEAN DIAMETER',',',ES12.2E3)") t,SUMAGG,SUMPP,DIAM_MEAN
	 t=t+tstep
	
	print*, ' ' 
	! write(*,"('AGG',3X,I1,3X,ES8.2E3,',',3X,'PP',3X,I1,3X,ES8.2E3,
!	1  ',',3X,'Time',3X,ES8.2E3)") J,SSEC(J),J,SSEC(J+MSECTION),t 

	end do
	print*, 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        print*,'the new aggregates:        ', SSEC(1:MSECTION)
	print*,'the new primary particles: ', SSEC(MSECTION+1:2*MSECTION)
C
	end program aggdrive


