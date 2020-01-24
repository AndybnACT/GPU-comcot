!******************* MODIFICATION RECORDS *****************************
!
! 1.  COMCOT (CORNELL MULTIGRID COUPLED TSUNAMI MODEL) VERSION 1.4
!     ALL RIGHTS RESERVED.
! 2.  REVISED BY S. N. SEO BASED ON SHUTO'S MODEL (AUG. 10, 1993)
! 3.  VERSION 1.0 ARRANGED BY YONGSIK CHO (SEP. 15, 1993)
! 4.  NONLINEAR MODEL, AUTOMATIC INITIAL SURFACE INTERPOLATION,
! 5.  GENERAL GRID MATCHING, USER INTERFACE IS INTRODUCED
!     BY SEUNG-BUHM WOO (1999, 3.)
! 6.  GRID AND FAULT PARAMETER MODULES INTRODUCED BY TOM LOGAN (ARSC,
!	  NOV., 2002)
!
!.....MODIFICATIONS/IMPROVEMENTS (XIAOMING WANG)
!      1. FIXED BUGS IN CONNECTING BOUNDARY INTERPOLATION
!		  (ARRAY OVERFLOW), 03/21/2003
!      2. APPLY BOTH CARTESIAN AND SPHERICAL COORDINDATE SYSTEMS (FOR
!		  OUTEST LAYER), 05/20/2003
!      3. ADD BOTTOM MOTION IN LINEAR MODEL FOR TENTATIVE TRY, JUN20 2003
!      4. FIX ERRORS IN INITIAL SURFACE CALCULATION (SMYLIE'S MODEL)
!      5. MODIFICATION ON INITIAL SURFACE PROFILE CALCULATION,
!		  07/13/2003
!      6. ADD OKADA'S MODEL TO GENERATE INITIAL SURFACE PROFILE,
!			08/19/2003
!         APPLY MULTIPLE LAYERS OF THE SAME LEVEL (09/04/2003)
!      7. MIGRATION FROM FORTRAN 77 TO FORTRAN 90 (10/07/2003)
!         MAKE EVERYTHING PARAMETER-DRIVEN
!      8. COMBINE ALL SIMILAR SUBROUTINES (11/19/2003)
!      9. REORGANIZE NONLINEAR PARTS  (BASED ON IMPROVEMENTS BY
!			TOM LOGAN, ARSC, 12/29/2003)
!     10. USE ARBITRARY INITIAL SURFACE PROFILE VIA AN EXTERNAL
!		  FILE (01/23/2004)
!     11. REFINE MOVING BOUNDARY ALGRITHM (02/12/2004)
!     12. TESTING MOVING BOUNDARY AND REMOVE UNNECESSARY VARIABLES
!		  (02/28/2004)
!     13. ADD WAVE MAKER SUBROUTINE TO GENERATE SINE AND SOLITARY WAVES
!         ADD BOLIQUE INCIDENT WAVE INTO WAVE MAKER
!     14. ADD ARBITRAY TIME SERIES OF ELEVATION INTO WAVE MAKER
!		  (FOR BENCHMARK IN CATALINA LONG WAVE WORKSHOP 2004)
!     15. REMOVE SUBROUTINES: PQSTILL,PQTOTAL,INI_MV;
!		  REFINE CONMOME SUBROUTINE (MAY 29,2004)
!     16. IN SUBROUTINE JNQ, INTERPOLATE FLUX INSTEAD OF VELOCITY
!		  (VELOCITY DISCONTINUOUS ACROSS INTERFACE!!!)  (NOV 2004)
!     17. INCREASE RUNUP CUTOFF (ELMAX) FROM 20M TO 50M DUE TO 2004
!		  TSUNAMI (JAN 2005)
!     18. ADD SUBMARINE LANDSLIDE SUBROUTINE INTO COMCOT TO INCLUDE
!		  TRANSIENT MOTION OF SEAFLOOR, WHICH WAS USED TO SIMULATE 2004
!		  INDIAN OCEAN TSUNAMI BASED ON CHEN JI'S TRANSIENT FAULT MODEL
!		 (MARCH 2005)
!     19. REFINE WAVE MAKER SUBROUTINE, REMOVE BUGS FOUND
!         AND ADD THE CAPABILITY TO GENERATE FOCUSING WAVES
!		  (SEPTEMBER 2006)
!     20. IMPROVE THE MAPPING ALGORITHM OF INITIAL DISPLACEMENT FROM
!		  SPHERICAL COORDINATES TO CARTESIAN COORDINATES. THE ORIGINAL
!		  ONE IS VERY ROUGH AND WILL RESULT IN A LARGE LOCATION
!		  DISCREPANCY AS LATITUDE INCREASES.
!         (TOGETHER WITH TSO-REN WU, NCU, OCTOBER 2006)
!.........
!     COMCOT V1.6 WAS CREATED AND MAINTAINED BY XIAOMING WANG.
!	  FOR ADDITIONAL INFO AND BUGS/ERRORS FOUND, PLEASE CONTACT HIM AT
!	  XW46@CORNELL.EDU
!.........
!     21. SEVERAL BUGS WERE FOUND IN VERSION 1.6, E.G., IN BOTTOM
!		  FRICTION CALCULATION AND NONLINEAR SUBROUTINES. (JAN 2007)
!     22. TO ENHANCE THE STABILITY ALONG THE INTERFACE OF TWO GRIDS,
!		  LINEAR SWE IS ADOPTED TO SOLVE THE FIRST 3 GRIDS INSIDE THE
!		  INTERFACE BETWEEN THE GRID WITH LINEAR SWE AND THE GRID USING
!		  NONLINEAR SWE WHEN NESTED GRID IS IMPLEMENTED (MARCH 2007)
!     23. ADD ONE SUBROUTINE FOR COURANT NUMBER CHECK (MARCH, 2007)
!     23A. INCREASE OUTPUT ACCURACY FROM 15F8.2 TO 15F9.3 (MARCH, 2007)
!     24. ADD A NEW FUNCTION - SPATIAL VARIATION OF MANNING'S ROUGHNESS
!		  COEF.; NOW MANNING'S ROUGNESS CAN BE EITHER A CONSTANT OR
!		  VARIABLE IN SPACE (AS JANAKA WIJETUNGE REQUIRES, APR. 2007)
!     25. ENABLE AN OPTION FOR LATE START: THE SIMULATION WILL SAVE
!		  RESUMING SNAPSHOTS EVERY 1000 STEPS; SIMULATION COULD BE
!		  RESUMED FROM ANY OF THESE SAVED SNAPSHOTS. (APR. 2007)
!     26. ADD MANNING'S FORMULA (BOTTOM FRICTION) TO LINEAR MODEL
!		  (APR. 2007)
!     27. IT WAS ASSUMED THAT THE GRID CELL WAS A SQUARE AND CLOSE TO A
!		  SQUQRE WHEN THE MIXING OF SPHERICAL AND CARTESIAN COORD. WERE
!		  ADOPTED (FOR NESTED GRIDS). HOWEVER, LARGE DISCREPANCY OCCURS
!		  WHEN LATITUDE IS HIGH. FOR EXAMPLE, ONE ARC SECOND WILL BE
!		  LONGER IN EAST-WEST THAN IN NORTH-SOUTH AT HIGH LATITUDE.
!         NOW, DX AND DY ARE NO LONGER ASSUMED TO BE EQUAL IN CARTESIAN
!		  AND CALCULATED BASED ON LATITUDE. (MAY, 2007)

!.....DEC 20 2007: IMPROVE NUMERICAL DISPERSION OVER SLOWLY VARYING
!	               TOPOGRAPHY (CHO 1995, YOON,2002)
!.....JAN. 2008, ADD SPONGE LAYER SUPPORT.
!.....
!####### CONFIGURATION FILE COMCOT.CTL CHANGED FOR VERSION 1.7 ########
!####### OCT., 2008 (GNS)                                      ########
!//////////////////////////////////////////////////////////////////////
!
!     COMCOT V1.7 WAS CREATED AND MAINTAINED BY XIAOMING WANG.
!	  FOR ADDITIONAL INFO AND BUGS/ERRORS FOUND, PLEASE CONTACT HIM
!	  AT XW46@CORNELL.EDU / X.WANG@GNS.CRI.NZ
!
!//////////////////////////////////////////////////////////////////////
!%%%%%# THIS IS STILL A BETA VERSION UNDER DEVELOPMENT !!!!!!!!!!!
!     #. SOLVERS WERE REPACKAGED, NEW COMCOT.CTL INTRODUCED;
!     #. UP TO 12 LEVELS OF NESTED GRIDS CAN BE IMPLEMENTED;
!     #. AUTOMATIC NESTED-GRID MATCHING, SUB-LEVEL GRID REGION IS
!		 SPECIFIED BY COORDINATES;
!     #. MORE BATHYMETRY FORMAT INPUT; XYZ FORMAT PREFERRED;
!     #. MORE CONTROL ON BOUNDARY CONDITIONS;
!     #. OUTPUT OF TIMEHISTORY RECORDS AND MAX. FREE SURFACE
!		 DISPLACEMENT;
!     #. ADJUSTABLE GRID SIZE FOR SPHERICAL COORDINATES
!		 --> GENERATING SQUARE GRIDS;
!     #. TIME STEP SIZE RATIO OF PARENT TO CHILD GRID IS NO LONGER
!		 FIXED AS 2,
!        BUT DETERMINED BASED ON WATER DEPTH;
!     #. UP TO 99 FAULT PLANES CAN BE IMPLEMENTED AT DIFFERENT TIMES;
!     #. TSUNAMIS CAN BE GENERATED BY A COMBINATION OF FAULT PLANES AND
!		 SUBMARINE LANDSLIDES;
!     #. IMPROVED THE COUPLING SCHEME BETWEEN SPHERICAL COORDINATES AND
!		 CARTESIAN COORDINATES, CARTESIAN GRIDS CAN BE IMPLEMENTED IN
!		 HIGH LATITUDE WITHOUT DISTORTION (LESS EFFICIENT, BUT MORE
!		 ACCURATE IN COMPARISON WITH THE OLD SCHEME);
!
!//////////////////////////////////////////////////////////////////////
!*************************VARIABLE DECLARATION*************************
!*
!* NX   X-DIMENSION OF COMPUTATION DOMAIN
!* NY   Y-DIMENSION OF COMPUTATION DOMAIN
!* DX   GRID SIZE
!* DY   GRID SIZE
!* DT   TIME STEP
!* LAYCORD   - COORDINATES USED FOR CURRENT LAYER : SPHERICAL/CARTESIAN
!* LAYGOV    - GOVERNING EQ USED FOR CURRENT LAYER: LINEAR/NONLINEAR
!* LAYSWITCH - IF GRID LAYER IS INCLUDED IN A SIMULATION; 0:YES,1:NO
!* REL_SIZE  - GRID SIZE RATIO OF PARENT GRID TO ITS CHILD GRID
!*
!* Z(:,:,2) FREE SURFACE ELEVATION
!* M(:,:,2) VOLUME FLUX (OR DISCHARGE) IN X DIRECTION
!* N(:,:,2) VOLUME FLUX (OR DISCHARGE) IN Y DIRECTION
!* HM(:,:,2) VOLUME FLUX AT DISCHARGE POINT
!* HN(:,:,2) VOLUME FLUX AT DISCHARGE POINT
!*
!* DZ(:,:,2) - TOTAL WATER DEPTH (USED FOR DETERMINING MOVING BOUNDARY)
!* FILE IO-UNIT 23, 25, 666, 999, and 1001-9999 RESERVED.
!**********************************************************************

!----------------------------------------------------------------------
      PROGRAM COMCOT
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE WAVE_PARAMS
	  USE FAULT_PARAMS
	  USE LANDSLIDE_PARAMS
	  USE BCI_PARAMS
      TYPE (LAYER)  :: LO
!	  TYPE (LAYER), DIMENSION(12)  :: LA
	  TYPE (LAYER), ALLOCATABLE  :: LA(:)
	  TYPE (WAVE)       :: WAVE_INFO
!	  TYPE (FAULT), DIMENSION(99)  :: FAULT_INFO
	  TYPE (FAULT), ALLOCATABLE  :: FAULT_INFO(:)
	  TYPE (LANDSLIDE)  :: LANDSLIDE_INFO
	  TYPE (BCI)        :: BCI_INFO
      INTEGER           :: INI_SURF
	  INTEGER           :: BC_TYPE,OUT_OPT
	  INTEGER           :: KEND,IPRT,STAT
	  INTEGER           :: START_TYPE,START_STEP,ISTART
      REAL              :: TEND, T_INV, H_LIMIT
	  REAL, ALLOCATABLE :: TS_LOC(:,:)			! FOR TIMEHISTORY RECORDS
	  INTEGER, ALLOCATABLE :: TS_ID(:)			! FOR TIMEHISTORY RECORDS
	  REAL TS_DAT								! FOR TIMEHISTORY RECORDS
	  CHARACTER(LEN=20),ALLOCATABLE :: FNAME(:)	! FOR TIMEHISTORY RECORDS
	  CHARACTER(LEN=200) JOB							! FOR SIMULATION JOB DESCRIPTION
	  !REAL              :: VERSION=1.7
      INTEGER :: RSTAT=0
      REAL FI, ST
      integer(kind=8) :: t1, t2, trate, tmax
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
!
	  ALLOCATE(LA(NUM_GRID))
	  ALLOCATE(FAULT_INFO(NUM_FLT))
!----------------------------------------------------------------------
      write (*,*) '  '
      write (*,*) '************** GPU-COMCOT ******************'
      write (*,*) '*                                          *'
      write (*,*) '*            VERSION= G1.7                 *'
      write (*,*) '*                                          *'
	  write (*,*) '********************************************'
!----------------------------------------------------------------------
!///// READ USER-DEFINED OPTIONS (NO. OF LAYERS, COORDINATE, GOVERN EQN.)
      CALL READ_CONFIG (LO,LA,TEND,T_INV,INI_SURF,WAVE_INFO,		&
	                   FAULT_INFO,LANDSLIDE_INFO,BCI_INFO,			&
					   START_TYPE,START_TIME,BC_TYPE,OUT_OPT,JOB)

!///// READING BATH/TOPO DATA /////////////////////////////////////////
      CALL READ_BATHYMETRY(LO,LA)

!///// READ INITIAL CONDITION /////////////////////////////////////////
	  IF (BC_TYPE.NE.3) THEN
	     CALL GET_INI_SURF(LO,LA,INI_SURF,WAVE_INFO,FAULT_INFO,		&
													LANDSLIDE_INFO)
	  ENDIF

!///// ADJUST BATHYMETRY, TIME STEP AND CALC DEPENDENT PARAMETERS /////
	  CALL ADJUST_BATHYMETRY (LO,LA)
	  CALL CR_CHECK (LO,LA)
	  CALL ALPHA_CALC (LO,LA)
	  !CALCULATE PARAMETERS USED FOR SPHERICAL COORD.
      CALL SPHERICAL_PARAMETERS(LO,LA)

!///// SET UP BOUNDARY CONDITIONS /////////////////////////////////////
	  IF (BC_TYPE.EQ.1) CALL SPONGE_COEF (LO)
	  IF (BC_TYPE.EQ.2) CALL BC_WALL (LO,WAVE_INFO)
	  IF (BC_TYPE.EQ.3) CALL GET_BC_DATA (BCI_INFO,LO)

!/////DETERMINE STARTING TIME STEP # //////////////////////////////////
!      WRITE(*,*) START_STEP
	  START_STEP = NINT(START_TIME/LO%DT)
	  ISTART = START_STEP
      IF (START_TYPE .EQ. 0) ISTART = 1

!//////// DISPLAY INPUT PARAMETERS ////////////////////////////////////
      CALL SCREEN_DISP (LO,LA,FAULT_INFO,LANDSLIDE_INFO,WAVE_INFO,	&
						INI_SURF,TEND,T_INV,KEND,IPRT,START_TYPE,	&
						START_TIME,BC_TYPE,BCI_INFO,OUT_OPT,JOB)
!//////// CREATE DATA FILES FOR TIME SERIES RECORDS OUTPUT ////////////
!	  OPEN DATA FILE TIME.DAT TO WRITE TIME SEQUENCE
	  OPEN (UNIT=99,FILE='time.dat',STATUS='UNKNOWN')
!     TIME HISTORY RECORDS CAN BE CREATED FOR UP TO 9999 LOCATIONS
      IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) THEN
	     NUM_REC = -1
	     TS_DAT = 0.0
         OPEN (UNIT=25,FILE='ts_location.dat',STATUS='OLD',IOSTAT=STAT)
         IF (STAT /=0) THEN
            PRINT *,"ERROR:: CAN'T OPEN TS_LOCATION.DAT; EXITING."
            STOP
         END IF
		 WRITE (*,*) '    OUTPUT TIME HISTORY RECORDS AT:'
         DO WHILE (RSTAT == 0)
		    NUM_REC = NUM_REC + 1
			READ (25,*, IOSTAT = RSTAT) TEMP1,TEMP2
	     ENDDO
		 REWIND(25)
		 ALLOCATE(TS_LOC(NUM_REC,2))
		 ALLOCATE(TS_ID(NUM_REC))
		 ALLOCATE(FNAME(NUM_REC))
		 TS_LOC = 0.0
		 TS_ID = 0
	     FNAME = ' '
		 DO K = 1,NUM_REC
		    READ (25,*) TS_LOC(K,1),TS_LOC(K,2)
			IF (LO%LAYCORD.EQ.0 .AND. TS_LOC(K,1).LT.0.0) THEN
			   TS_LOC(K,1) = 360.0 + TS_LOC(K,1)
			ENDIF
			WRITE (*,*) '       ',TS_LOC(K,1),TS_LOC(K,2)
			CALL GAUGE_LAYERID (ID,TS_LOC(K,1),TS_LOC(K,2),LO,LA)
			TS_ID(K) = ID
		 ENDDO
         CLOSE(25)
		 !OPEN DATA FILES TO WRITE TIME HISTORY DATA
         DO NIC = 1,NUM_REC
            WRITE (FNAME(NIC),1) NIC
 1          FORMAT('ts_record',I4.4,'.dat')
            IO_UNIT = 1000+NIC
            OPEN (IO_UNIT,FILE=FNAME(NIC),STATUS='UNKNOWN')
         ENDDO
	  END IF
!......................................................................
	  !PAUSE

      CALL GCOMCOT_INIT_GPU()
      CALL GCOMCOT_INIT_LAYER(LO%ID, LO%PARENT, LO%LEVEL,           &
                              LO%R1, LO%R2, LO%R3, LO%R4, LO%R5,    &
                              LO%R6, LO%R11, LO%H, LO%Z, LO%NX, LO%NY)

!////////////////////// SIMULATION BEGINS /////////////////////////////
      WRITE (*,*) '    '
	  WRITE (*,*) '***************** OUTPUT RESULTS ******************'
      WRITE (*,*) 'TIMESTEP          MINUTE'
	  TIME = DBLE(ISTART-1)*LO%DT
!	  T_MINUTES = TIME/60.0

!.....OUTPUT TIME SEQUENCE AND TIME SERIES DATA AT T = 0.0
	  IF (ISTART.EQ.1) THEN
		 WRITE (99,*) TIME
		 CALL ALL_PRT (ISTART-1,LO,LA)
		 IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) THEN
            DO NIC = 1,NUM_REC
               WRITE(*,*) "[WARNING] THE FUNCTION [GET_TS] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT----->[will be updated soon]"
		       CALL GET_TS (TS_DAT,TS_LOC(NIC,1),TS_LOC(NIC,2),		&
												LO,LA,TS_ID(NIC))
!*			  WRITE(*,*) TS_LOC(NIC,1),TS_LOC(NIC,2),TS_DAT
			   IO_UNIT = 1000+NIC
		       WRITE(IO_UNIT,'(F10.5)') TS_DAT
            ENDDO
		 ENDIF
	  ENDIF

      DO K = ISTART,KEND
         TIME = TIME+LO%DT
		 T_MINUTES = TIME/60.0
         IF (MOD(K,10) .EQ. 0) THEN
            WRITE (*,'(I9,8X,F8.3)') K,T_MINUTES
         ENDIF

!.... .. CALL HOT START FUNCTION
         IF (K .EQ. START_STEP) THEN
            WRITE(*,*) "[WARNING] THE FUNCTION [HOT_START] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
		    CALL HOT_START (LO,LA,START_TYPE,START_STEP)
		 ENDIF
!*		 ! CREATE 5 BACKUPS DURING THE ENTIRE SIMULATION DURATION
!*		 IF (K.GT.1 .AND. MOD(K,FLOOR(KEND/5.)) .EQ. 0) THEN
!*		    CALL HOT_START (LO,LA,0,K)
!*		 ENDIF

!....... CALL MULTIPLE FAULT PLANE MODEL
         IF (INI_SURF .EQ. 0 .OR. INI_SURF .EQ. 4) THEN
		    CALL GET_FLOOR_DEFORM (LO,LA,FAULT_INFO,TIME)
		 ENDIF
!... ...CALL LANDSLIDE MODEL
		 IF ( INI_SURF .EQ. 3 .OR. INI_SURF .EQ. 4 ) THEN
		    IF ( TIME.GE.LANDSLIDE_INFO%T(1) .AND.					&
		         TIME.LE.LANDSLIDE_INFO%T(LANDSLIDE_INFO%NT) ) THEN
               WRITE(*,*) "[WARNING] THE FUNCTION [LAND_SLIDE] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
		       CALL LAND_SLIDE(LO,LANDSLIDE_INFO,TIME)
		    ENDIF
	     ENDIF

!////////SOLVE MASS CONSERVATION EQN FOR LAYER 1 (THE OUTEST LAYER)
         !!....CALL SEDIMENT TRANSPORT MODEL
         IF (LO%LAYCORD.NE.0 .AND. LO%SEDI_SWITCH .EQ. 0) THEN
            WRITE(*,*) "[WARNING] THE FUNCTION [SED_TRANSPORT] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
		    CALL SED_TRANSPORT (LO)
		 ENDIF

         ! WRITE(*,*) "COMPUTING MASS EQ ON CPU"
         ! call system_clock ( t1, trate, tmax )
		 ! CALL MASS (LO)
         ! call system_clock ( t2, trate, tmax )
         ! write ( *, * ) 'Serial Code Elapsed real time = ', real ( t2 - t1 ) / real ( trate )
         ! WRITE (*,*) "COMPUTING MASS EQ ON GPU"
         CALL MASS_LAUNCH(LO%Z(:,:,1),LO%Z(:,:,2),LO%H(:,:), LO%ID)

!.......SOLVE RADIATION OPEN BOUNDARY
         ! CALL OPEN (LO)
         CALL OPENBD_LAUNCH(LO%Z(:,:,2))

!.......CALL WAVE MAKER TO GENERATE WAVES
         IF (INI_SURF .EQ. 2) THEN
            WRITE(*,*) "[WARNING] THE FUNCTION [WAVE_MAKER] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
            CALL WAVE_MAKER (TIME,LO,WAVE_INFO)
		 ENDIF

!.......APPLY FACTS INPUT CONDITION
         IF (BC_TYPE.EQ.3) THEN
		    IF (TIME.GE.BCI_INFO%T(1) .AND. 						&
							TIME.LT.BCI_INFO%T(BCI_INFO%NT)) THEN
               WRITE(*,*) "[WARNING] THE FUNCTION [BC_INPUT] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
		       CALL BC_INPUT (BCI_INFO,LO,TIME)
			ENDIF
		 ENDIF
!//////////////////////////////////////////////////////////////////////
!        CALL SUBLAYER COUPLED MODEL
!//////////////////////////////////////////////////////////////////////
         IF (LO%NUM_CHILD .GE. 1) THEN
             WRITE(*,*) "[WARNING] THE MULTI-LAYERS FUNCTIONALITY IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
             CALL ALL_GRID (LO,LA)
         ENDIF
!//////////////////////////////////////////////////////////////////////
!        SOLVE MOMENTUM CONSERVATION EQN FOR LAYER 1
!//////////////////////////////////////////////////////////////////////
         ! WRITE(*,*) "COMPUTING MOMT EQ ON CPU"
         ! call system_clock ( t1, trate, tmax )
         ! CALL MOMENT (LO)
         ! call system_clock ( t2, trate, tmax )
         ! write ( *, * ) 'Serial Code Elapsed real time = ', real ( t2 - t1 ) / real ( trate )
         !
         ! WRITE(*,*) "COMPUTING MOMT EQ ON GPU"
         CALL MOMT_LAUNCH(LO%M(:,:,2), LO%N(:,:,2), LO%Z(:,:,2), LO%ID)

!.......USE SPONGE LAYER .....
         IF (BC_TYPE.EQ.1) THEN
             WRITE(*,*) "[WARNING] THE FUNCTION [SPONGE_LAYER] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
             CALL SPONGE_LAYER(LO)
         ENDIF
!*!.......WALL BOUNDARY
         IF (BC_TYPE.EQ.2) THEN
             WRITE(*,*) "[WARNING] THE FUNCTION [BC_WALL] IS CURRENTLY NOT SUPPORTED ON GPU COMCOT"
             CALL BC_WALL (LO,WAVE_INFO)
         ENDIF
!......................................................................
!       OUTPUT RESULTS AT GIVEN TIME INTERVAL (MINUTE)
!......................................................................
         IF (MOD(K,IPRT) .EQ. 0) THEN
            CALL ALL_PRT (K,LO,LA)
         ENDIF


		 IF (OUT_OPT.EQ.0 .OR. OUT_OPT.EQ.2) THEN
            CALL MAXAMP_LAUNCH(LO%ID)
		    CALL MAX_AMP (LO,LA,TIME,TEND)
		 ENDIF


!.......UPDATE VARIABLES OF LAYER 01 (LO) FOR NEXT TIME STEP
         CALL CHANGE()
!.......OUTPUT TIME SEQUENCE AND TIME HISTORY RECORDS AT T = K*LO%DT
!!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
!          WRITE (99,*) TIME                                                !
! 		 IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) THEN                           !
!             DO NIC = 1,NUM_REC                                            !
! 		       CALL GET_TS (TS_DAT,TS_LOC(NIC,1),TS_LOC(NIC,2),		&       !
! 												LO,LA,TS_ID(NIC))           !
! !*			  WRITE(*,*) TS_LOC(NIC,1),TS_LOC(NIC,2),TS_DAT             !
! 			   IO_UNIT = 1000+NIC                                           !
! 		       WRITE(IO_UNIT,'(F10.5)') TS_DAT                              !
!             ENDDO                                                         !
! 		 END IF                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       END DO

!..... CLOSE ALL OPENNED DATA FILES
	  CLOSE(99)
	  IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) THEN
	     DO NIC = 1,NUM_REC
	        CLOSE(1000+NIC)
	     ENDDO
      END IF

      END
!//////////////// END OF MAIN PROGRAM /////////////////////////////////

!//////////////////////////////////////////////////////////////////////
      BLOCK DATA
         COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					   NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
         DATA ELMAX    / -50.0 /			! UPPER LIMIT OF MAXIMUM RUNUP !
         DATA GRAV     / 9.807 /			! GRAVITY
         DATA PI       / 3.14159265358979 /
	     DATA R_EARTH  / 6378000.0 /		! RADIUS OF THE EARTH
	     DATA GX       / 1.0E-5 /			! FOR INUNDATION, DZ=0.0 IF DZ < GX
	     DATA EPS      / 1.0E-10 /			! Z,U,V = 0.0 IF < EPS
	     DATA ZERO     / 0.0 /
	     DATA ONE      / 1.0 /
	     DATA NUM_GRID / 12 /				! MAX. NUMBER OF SUB-LEVEL GRIDS
		 DATA NUM_FLT  / 999 /				! TOTAL NO. OF FAULT SEGMENTS
	     DATA V_LIMIT  / 20.0 /				! UPPER LIMIT OF MAXIMUM VELOCITY !
	     DATA RAD_DEG  / 0.01745329252 /	! ARC RADIAN OF 1 DEGREE
	     DATA RAD_MIN  / 0.000290888208665721 /  ! ARC RADIAN OF 1 MINUTE
      END
!//////////////////////////////////////////////////////////////////////


!----------------------------------------------------------------------
      SUBROUTINE CHANGE()
!.....TRANSFER INFORMATION (FREE SURFACE ELEVATION, VOLUME FLUX IN
!      X, Y DIRECTIONS)FROM LAST STEP TO NEXT STEP (FOR OUTEST LAYER)
!----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
!     USE LAYER_PARAMS
!     TYPE (LAYER)	:: LO
!       IF (LO%LAYGOV .GT. 1) THEN
!           LO%M0(:,:)   = LO%M(:,:,1)
!           LO%N0(:,:) = LO%N(:,:,1)
!       ENDIF
!       LO%Z(:,:,1) = LO%Z(:,:,2)
!       LO%M(:,:,1) = LO%M(:,:,2)
!       LO%N(:,:,1) = LO%N(:,:,2)
!       ! DO J=1,LO%NY
!       !    DO I=1,LO%NX
!       !       LO%Z(I,J,1) = LO%Z(I,J,2)
!       !       LO%M(I,J,1) = LO%M(I,J,2)
!       !       LO%N(I,J,1) = LO%N(I,J,2)
! 		! 	IF (LO%LAYGOV.GT.1) THEN
! 		!        LO%M0(I,J)  = LO%M(I,J,1)
! 		!        LO%N0(I,J)  = LO%N(I,J,1)
! 		! 	ENDIF
!       !    END DO
!       ! END DO
! !.....UPDATE BATHYMETRY IF TRANSIENT SEAFLOOR MOTION IS ENABLED
! 	  IF (LO%INI_SWITCH.EQ.3 .OR. LO%INI_SWITCH.EQ.4) THEN
! 	     LO%HT(:,:,1) = LO%HT(:,:,2)
! 	     LO%H(:,:) = LO%H(:,:) + LO%HT(:,:,2) - LO%HT(:,:,1)
! 	  ENDIF
!!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL CUDA_UPDATE_LAYER(1)
!
      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE ALLOC (L,SWITCH)
!     CREATED BY TOM LOGAN, UAF (2005)
!     UPDATED ON SEP 17 2006 (XIAOMING WANG,CORNELL UNIVERSITY)
! SWITCH:
!       1 -- TOP-LEVEL GRID (OUTEST LAYER)
!       OTHER -- SUB-LEVEL GRID LAYERS
!       ...
! LINCHK:
!		0 -- NONLINEAR TERMS WON'T BE CALCULATED IN NSWE
!		1 -- NONLINEAR TERMS WILL BE CALCULATED IN NSWE
! MODSCM:
!		0 -- MODIFIED SCHEME IS IMPLETMENTED (IMAMURA,1987;CHO,1995)
!		1 -- MODIFIED SCHEME IS NOT IMPLETMENTED
!.....UPDATED ON NOV 03 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE(LAYER)               :: L
      INTEGER		:: SWITCH	! 1 - TOP LAYER; ELSE - OTHER LAYER
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      L%RX = L%DT/L%DX
	  L%RY = L%DT/L%DY
      L%GRX = GRAV*L%RX
	  L%GRY = GRAV*L%RY

      IF (SWITCH .EQ. 1) THEN
	      L%LINCHK = 1
	      L%MODSCM = 0
		  L%CORNERS(1) = 1
		  L%CORNERS(2) = L%NX
		  L%CORNERS(3) = 1
		  L%CORNERS(4) = L%NY
      ELSE
          L%LINCHK = 1
      	  L%MODSCM = 1
      END IF

      IF ( L%LAYSWITCH.EQ.0 ) THEN
		ALLOCATE(L%Z(L%NX,L%NY,2))
		ALLOCATE(L%M(L%NX,L%NY,2))
		ALLOCATE(L%N(L%NX,L%NY,2))
		ALLOCATE(L%H(L%NX,L%NY))
		ALLOCATE(L%XFLUX(L%NX,2))
		ALLOCATE(L%YFLUX(L%NY,2))

		ALLOCATE(L%HP(L%NX,L%NY))
		ALLOCATE(L%HQ(L%NX,L%NY))
		ALLOCATE(L%DZ(L%NX,L%NY,2))
		ALLOCATE(L%DEFORM(L%NX,L%NY))
		IF (L%LAYCORD.EQ.0) THEN
           ALLOCATE(L%R0(L%NX,L%NY))
		   ALLOCATE(L%R1(L%NX,L%NY))
		   ALLOCATE(L%R11(L%NX,L%NY))
		   ALLOCATE(L%R2(L%NX,L%NY))
		   ALLOCATE(L%R21(L%NX,L%NY))
		   ALLOCATE(L%R22(L%NX,L%NY))
		   ALLOCATE(L%R3(L%NX,L%NY))
		   ALLOCATE(L%R4(L%NX,L%NY))
		   ALLOCATE(L%R5(L%NX,L%NY))
		   ALLOCATE(L%R6(L%NY))
		ENDIF

		IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4) THEN
		   ALLOCATE(L%HT(L%NX,L%NY,2))
		   L%HT = 0.0
		ENDIF

		IF (L%LAYGOV.GE.1) THEN
		   ALLOCATE(L%MASK(L%NX,L%NY))
		   L%MASK = 0
		ENDIF
		IF (L%LAYGOV.GT.1) THEN
		   ALLOCATE(L%M0(L%NX,L%NY))
		   ALLOCATE(L%N0(L%NX,L%NY))
!		   ALLOCATE(L%MASK(L%NX,L%NY))
		   L%M0 = 0.0
		   L%N0 = 0.0
!		   L%MASK = 0
		ENDIF

		IF (L%FRIC_SWITCH.EQ.2) THEN
		   ALLOCATE(L%FRIC_VCOEF(L%NX,L%NY))
		   L%FRIC_VCOEF = L%FRIC_COEF
		ENDIF
!*		IF (L%LAYGOV.EQ.1 .OR. L%LAYGOV.EQ.9) L%MODSCM = 1
		IF (MOD(L%LAYGOV,2).NE.0) L%MODSCM = 1

		ALLOCATE(L%Z_MAX(L%NX,L%NY))
		ALLOCATE(L%Z_MIN(L%NX,L%NY))

!*		IF (L%UPZ .EQ. .FALSE.) THEN
!*		   ALLOCATE(L%POS(L%NX,L%NY,2))
!*		   ALLOCATE(L%CXY(L%NX,L%NY,4))
!*		   L%POS = 0
!*		   L%CXY = 0.0
!*		ENDIF

		L%Z = 0.0
		L%M = 0.0
		L%N = 0.0
		L%H = 0.0
		L%HP = 0.0
		L%HQ = 0.0
		L%DZ = 0.0

		L%DEFORM = 0.0
		L%XFLUX = 0.0
		L%YFLUX = 0.0
!		L%NAME = ''
!.......UPZ - FLAG TO IDENTIFY COORDINATE SYSTEM
!			  .TRUE. - USE THE SAME COORD. AS ITS PARENT GRID
!			  .FALSE. - USE DIFFERENT COORD.
!		L%UPZ = .TRUE.

		L%Z_MAX = 0.0
		L%Z_MIN = 0.0
!.......L%X,L%Y,L%DX,L%DY,L%DT ARE ALLOCATED AND
!		INITIALIZED IN DX_CALC AND SUBGRID_MATCHING
	  END IF

	  IF (L%LAYSWITCH.EQ.0 .AND. L%LAYGOV.GE.2) THEN
		ALLOCATE(L%ALPHA(L%NX,L%NY))
		ALLOCATE(L%A1X(L%NX,L%NY))
		ALLOCATE(L%A2X(L%NX,L%NY))
		ALLOCATE(L%A1Y(L%NX,L%NY))
		ALLOCATE(L%A2Y(L%NX,L%NY))
		ALLOCATE(L%CF(L%NX,L%NY,4))
		ALLOCATE(L%CB(L%NX,L%NY,4))
		L%ALPHA  = 1.0
		L%CF = 0.0
		L%CB = 0.0
	  END IF

	  IF (SWITCH .EQ. 1 .AND. L%BC_TYPE.EQ.1) THEN
		 ALLOCATE(L%SPONGE_COEFX(L%NX,L%NY))
	     ALLOCATE(L%SPONGE_COEFY(L%NX,L%NY))
		 L%SPONGE_COEFX = 0.0
		 L%SPONGE_COEFY = 0.0
	  ENDIF

	  IF (L%FRIC_SWITCH.EQ.3) THEN
	     L%SEDI_SWITCH = 0
		 ALLOCATE(L%DH(L%NX,L%NY,2))
		 L%DH = 0.0
	  ELSE
         L%SEDI_SWITCH = 1
	  ENDIF

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE SPHERICAL_PARAMETERS(LO,LA)
!     CREATED BY TOM LOGAN, UAF (2004)
!  ******** PARAMETERS FOR SPHERICAL COORD. *******
!  Updated on Dec21 2008 (Xiaoming Wang, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      IF (LO%LAYSWITCH .EQ. 0 .AND. LO%LAYCORD .EQ. 0) THEN
		 WRITE (*,*) '    CALCULATING COEFFICIENTS FOR SPHERICAL COORD.'
		 CALL PARA(LO)
	  ENDIF
	  DO I=1,NUM_GRID
        IF (LA(I)%LAYSWITCH .EQ. 0 .AND. LA(I)%LAYCORD .EQ. 0) THEN
		   CALL PARA(LA(I))
		ENDIF
      END DO
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE PARA(LO)
!......................................................................
!DESCRIPTION:
!	  #. CALCULATE PARAMETERS FOR SPHERICAL COORD.
!NOTE:
!	  #. REVISED SIGNIFICANTLY BY TOM LOGAN, UAF (2005)
!	  #. UPDATED ON NOV.07, 2008 (XIAOMING WANG, GNS)
!               LINES STARTING WITH "!!" ARE REMOVED
!               LINES ENDED WITH "!!"   ARE ADDED
!	  #. UPDATED ON FEB06 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: LO
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA CORI_W/7.2722E-5/, OSIXTY/0.016666666667/
	  DATA TWLVTH/0.08333333333333/

!.....PARAMETERS FOR SPHERICAL COORD IN GRID LAYER LO
	  IF (LO%LAYCORD .EQ. 0) THEN
	     IF (LO%LAYGOV.LE.1 .OR. LO%LAYGOV.EQ.9) THEN
!*          DX = R_EARTH*LO%DX*RAD_MIN
!*          RR = LO%DT/DX
            RR = LO%DT/R_EARTH
!*          RS = 0.5*G*RR
            RS = GRAV*RR             !!
            RT = 0.5*LO%DT*CORI_W
            DO J=1,LO%NY
		       ANG_M = LO%Y(J)*RAD_DEG
		       ANG_N = (LO%Y(J)+0.5*LO%DEL_Y(J)*OSIXTY)*RAD_DEG
               COS_M = COS(ANG_M)
               COS_N = COS(ANG_N)
               SIN_M = SIN(ANG_M)
               SIN_N = SIN(ANG_N)
               DO I=1,LO%NX
 		          LO%R0(I,J) = RR/(LO%DEL_Y(J)*RAD_MIN)       !!
			      LO%R1(I,J) = RR/COS_M/(LO%DX*RAD_MIN)
			      LO%R11(I,J) = RR/COS_M/(LO%DEL_Y(J)*RAD_MIN)
			      IF (LO%LAYGOV.EQ.0) THEN
					 LO%R2(I,J) = RS/COS_M/(LO%DX*RAD_MIN)			&
											*AMAX1(LO%HP(I,J),ZERO)
				  ENDIF
				  IF (LO%LAYGOV.EQ.1) THEN
					 LO%R2(I,J) = RS/COS_M/(LO%DX*RAD_MIN)
				  ENDIF
			      LO%R21(I,J) = RR/COS_M/(LO%DX*RAD_MIN)
				  LO%R22(I,J) = RR/COS_N/(LO%DX*RAD_MIN)
                  LO%R3(I,J) = RT*SIN_M
                  IF (LO%LAYGOV.EQ.0) THEN
					 LO%R4(I,J) = RS/(LO%DEL_Y(J)*RAD_MIN)			&
											*AMAX1(LO%HQ(I,J),ZERO)
				  ENDIF
				  IF (LO%LAYGOV.EQ.1) THEN
					 LO%R4(I,J) = RS/(LO%DEL_Y(J)*RAD_MIN)
				  ENDIF
                  LO%R5(I,J) = RT*SIN_N
               END DO
               LO%R6(J) = COS_N
            END DO
		 ELSEIF (LO%LAYGOV.EQ.2 .OR. LO%LAYGOV.EQ.3) THEN
!*          DX = R_EARTH*LO%DX*RAD_MIN
!*          RR = LO%DT/DX
            RR = LO%DT/R_EARTH
!*          RS = 0.5*G*RR
            RS = GRAV*RR             !!
            RT = 0.5*LO%DT*CORI_W
            DO J=1,LO%NY
		       ANG_M = LO%Y(J)*RAD_DEG
		       ANG_N = (LO%Y(J)+0.5*LO%DEL_Y(J)*OSIXTY)*RAD_DEG
               COS_M = COS(ANG_M)
               COS_N = COS(ANG_N)
               SIN_M = SIN(ANG_M)
               SIN_N = SIN(ANG_N)
               DO I=1,LO%NX
				  IF (LO%MASK(I,J).EQ.0) THEN
				     DX = LO%DX
				     DY = LO%DEL_Y(J)
				  ENDIF
				  IF (LO%MASK(I,J).EQ.1) THEN
				     DX = LO%DX*LO%ALPHA(I,J)
				     DY = LO%DEL_Y(J)*LO%ALPHA(I,J)
				  ENDIF
 		          LO%R0(I,J) = RR/(DY*RAD_MIN)       !!
			      LO%R1(I,J) = RR/COS_M/(DX*RAD_MIN)
			      LO%R11(I,J) = RR/COS_M/(DY*RAD_MIN)
			      IF (LO%LAYGOV.EQ.2) THEN
					 LO%R2(I,J) = RS/COS_M/(DX*RAD_MIN)				&
										*AMAX1(LO%HP(I,J),ZERO)
				  ENDIF
				  IF (LO%LAYGOV.EQ.3) THEN
					 LO%R2(I,J) = RS/COS_M/(DX*RAD_MIN)
				  ENDIF
			      LO%R21(I,J) = RR/COS_M/(DX*RAD_MIN)
				  LO%R22(I,J) = RR/COS_N/(DX*RAD_MIN)
                  LO%R3(I,J) = RT*SIN_M
                  IF (LO%LAYGOV.EQ.2) THEN
					 LO%R4(I,J) = RS/(DY*RAD_MIN)					&
										*AMAX1(LO%HQ(I,J),ZERO)
				  ENDIF
				  IF (LO%LAYGOV.EQ.3) THEN
					 LO%R4(I,J) = RS/(DY*RAD_MIN)
				  ENDIF
                  LO%R5(I,J) = RT*SIN_N
               END DO
               LO%R6(J) = COS_N
            END DO
		 ENDIF
      ENDIF

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE SED_TRANSPORT (LO)
!.....TENTATIVE TRY TO INCLUDING SEDIMENT TRANSPORT MODEL (APR. 2007)
!     SEDIMENT FLUX IS COMPUTED ACCORDING TO RIBBERINK (1998)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  REAL D, DS, THETA_C, THETA_P,THETA_Q,G,RHO,RHO_S
	  REAL H(LO%NX,LO%NY),DH(LO%NX,LO%NY)
	  REAL QX(LO%NX,LO%NY),QY(LO%NX,LO%NY)
	  CHARACTER(LEN=30) FNAME

	  G = 9.807
	  RHO = 1000.0
	  RHO_S = 2650.0
	  DS = 0.0002
	  POROSITY = 0.3
	  THETA_C = 0.04
	  S = RHO_S/RHO

!	  Q(:,:) = 0.0
	  QX(:,:) = 0.0
	  QY(:,:) = 0.0

!      H(:,:) = -LO%H(:,:)
	  DH(:,:) = LO%H(:,:) + LO%Z(:,:,1)

	  DO I = 1,LO%NX
	     IM1 = I-1
		 IP1 = I+1
		 IF (IM1.LE.0) IM1=1
		 IF (IP1.GE.LO%NX) IP1=LO%NX

	     DO J = 1,LO%NY
		 	JM1 = J-1
		    JP1 = J+1
		    IF (JM1.LE.1) JM1=1
		    IF (JP1.GE.LO%NY) JP1=LO%NY
			!TOTAL WATER DEPTH AT DISCHARGE LOCATION P
            HP = 0.5*(DH(I,J)+DH(IP1,J))
			!TOTAL WATER DEPTH AT DISCHARGE LOCATION Q
			HQ = 0.5*(DH(I,J)+DH(I,JP1))

			!SEDIMENT FLUX AT DISCHARGE LOCATION P
		    IF (HP .GT. 0.05) THEN
		       COEF_P = RHO*G*(LO%FRIC_COEF**2)/(HP**2.33333)
			   !Y COMPONENT OF VOLUME FLUX AT DISCHARGE LOCATION P
			   XQQ = 0.25*(LO%N(I,J,1)+LO%N(IP1,J,1)				&
						+ LO%N(I,JM1,1)+LO%N(IP1,JM1,1))
			   !TOTAL VOLUME FLUX AT DISCHARGE LOCATION P
			   FLUX_P = SQRT(LO%M(I,J,1)**2+XQQ**2)
			   !X COMPONENT OF BOTTOM SHEAR STRESS AT DISCHARGE LOCATION P
               TPX = COEF_P*LO%M(I,J,1)*FLUX_P
			   !Y COMPONENT OF BOTTOM SHEAR STRESS AT DISCHARGE LOCATION P
               TPY = COEF_P*XQQ*FLUX_P
			   ! TOTAL SHEAR STRESS AT DISCHARGE LOCATION P
			   TPB = SQRT(TPX**2+TPY**2)
			   !SHIELDS PARAMETER AT DISCHARGE LOCATION P
			   THETA_P = TPB/((RHO_S-RHO)*G*DS)
			   IF (THETA_P .GT. THETA_C) THEN
				  !RIBBERINK (1998)
			      QS_P = 11.0*SQRT((S-1.0)*G*DS**3)					&
							* ((THETA_P-THETA_C)**1.65)
				  IF (QS_P.GT.0.02) QS_P = 0.02
				  !X COMPONENT OF SAND FLUX AT DISCHARGE LOCATION P
				  IF (FLUX_P .GT. 1.0E-10)							&
							QX(I,J) = QS_P * LO%M(I,J,1)/FLUX_P
			   ENDIF
			ENDIF

			!SEDIMENT FLUX AT DISCHARGE LOCATION Q
            IF (HQ .GT. 0.05) THEN
		       COEF_Q = RHO*G*(LO%FRIC_COEF**2)/(HQ**2.33333)
			   XPP = 0.25*(LO%M(I,J,1)+LO%M(I,JP1,1)+LO%M(IM1,J,1)	&
							+ LO%M(IM1,JP1,1))
			   FLUX_Q = SQRT(XPP**2+LO%N(I,J,1)**2)
               TQX = COEF_Q*XPP*FLUX_Q
               TQY = COEF_Q*LO%N(I,J,1)*FLUX_Q
			   TQB = SQRT(TQX**2+TQY**2)
			   THETA_Q = TQB/((RHO_S-RHO)*G*DS)
			   IF (THETA_Q .GT. THETA_C) THEN
			      !RIBBERINK (1998)
			      QS_Q = 11.0*SQRT((S-1.0)*G*DS**3)					&
							*((THETA_Q-THETA_C)**1.65)
				  IF (QS_Q.GT.0.02) QS_Q = 0.02
				  IF (FLUX_Q .GT. 1.0E-10)							&
							QY(I,J) = QS_Q*LO%N(I,J,1)/FLUX_Q
			   ENDIF

            ENDIF
         ENDDO
      ENDDO

      DO I = 2,LO%NX-1
	     DO J = 2,LO%NY-1
		    LO%DH(I,J,2) = LO%DH(I,J,1)-LO%RX*(1/(1-POROSITY))		&
							*(QX(I,J)-QX(I-1,J)+QY(I,J)-QY(I,J-1))
         ENDDO
	  ENDDO
      LO%H(:,:) = LO%H(:,:)-LO%DH(:,:,2)


	  RETURN
	  END
