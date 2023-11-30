!-----------------------------------------------------------------------
      SUBROUTINE SCREEN_DISP (LO,LA,FAULT_INFO,LANDSLIDE_INFO,		&
	                          WAVE_INFO,INI_SURF,TEND,T_INV,KEND,	&
							  IPRT,START_TYPE,START_TIME,BC_TYPE,	&
							  BCI_INFO,OUT_OPT,JOB)
!DESCRIPTION:
!	  #. DISPLAY INFO ON SCREEN BEFORE COMPUTATION;
!	  #. SAVE SIMULATION INFO TO SIMULATION_INFO.TXT FOR CHECKING;
!NOTES:
!	  #. CREATED BY XIAOMING WANG (SEP 17 2006)
!	  #. UPDATED ON ??? ?? 2008
!	  #. UPDATED ON MAR11 2009 (XIAOMING WANG, GNS)
!		 1. ADD MORE OUTPUT INFO FOR FACTS INPUT
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE WAVE_PARAMS
	  USE LANDSLIDE_PARAMS
	  USE FAULT_PARAMS
	  USE BCI_PARAMS
	  TYPE (LAYER)  :: LO
	  TYPE (LAYER),DIMENSION(NUM_GRID) :: LA
	  TYPE (WAVE)   :: WAVE_INFO
	  TYPE (LANDSLIDE)  :: LANDSLIDE_INFO
!*	  TYPE (FAULT)      :: FAULT_INFO
	  TYPE (FAULT),DIMENSION(NUM_FLT)	:: FAULT_INFO
	  TYPE (BCI) BCI_INFO
	  REAL TEND,T_INV,H_LIMIT, START_TIME
	  INTEGER KEND,IPRT,INI_SURF,START_TYPE, BC_TYPE,OUT_OPT
	  CHARACTER(LEN=200)	   :: JOB
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      KEND = NINT(TEND/LO%DT)
	  IPRT = NINT(T_INV/LO%DT)

!.....DISPLAY INPUT INFORMATION ON SCREEN
	  WRITE (*,*) '***********************************************************'
	  WRITE (*,*) '*           INPUT INFORMATION - COMCOT V1.7               *'
	  WRITE (*,*) '***********************************************************'
!*	  WRITE (*,*) JOB
      WRITE (*,*) '------------------- GENERAL INFORMATION -------------------'
 	  WRITE (*,*) '    TOTAL RUN TIME            (SECOND) :',TEND
	  WRITE (*,*) '    TIME INTERVAL FOR OUTPUT  (SECOND) :',T_INV
	  IF (T_INV .LT. LO%DT) THEN
	     T_INV = LO%DT
		 IPRT = NINT(T_INV/LO%DT)
		 WRITE (*,*) '    OUTPUT INTERVAL < DT, RESET TO', T_INV
		 WRITE (*,*) '   '
      ENDIF
	  WRITE (*,*) '    TIME STEP SIZE            (SECOND) :',LO%DT
	  WRITE (*,*) '    TOTAL STEPS TO RUN         (STEPS) :',KEND
	  WRITE (*,*) '    STEP INTERVAL FOR OUTPUT   (STEPS) :',IPRT
	  IF (OUT_OPT.EQ.0 .OR. OUT_OPT.EQ.2) WRITE (*,*) '    MAX. SURFACE DISPLACEMENT OUTPUT   : ENABLED'
	  IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) WRITE (*,*) '    TIME HISTORY RECORD OUTPUT         : ENABLED'
	  IF (OUT_OPT .EQ. 1) WRITE (*,*) '    MAX. SURFACE DISPLACEMENT OUTPUT   : DISABLED'
	  IF (OUT_OPT .EQ. 0) WRITE (*,*) '    TIME HISTORY RECORD OUTPUT         : DISABLED'
	  IF (START_TYPE.EQ.1) THEN
	     WRITE (*,*) '    HOT START FUNCTION                 : ENABLED'
		 WRITE (*,*) '       SIMULATION RESUMES FROM T =     :',START_TIME
		 ISTEP_LEFT = NINT((TEND-START_TIME)/LO%DT)
		 WRITE (*,*) '       TOTAL NUMBER OF TIME STEPS LEFT :',ISTEP_LEFT
      ENDIF
	  WRITE (*,*) '    SHORELINE LOCATED AT DEPTH CONTOUR :',LO%H_LIMIT
	  IF (BC_TYPE.EQ.0) WRITE (*,*) '    BOUNDARY CONDITION                 : RADIATION (OPEN)'
	  IF (BC_TYPE.EQ.1) WRITE (*,*) '    BOUNDARY CONDITION                 : ABSORBING'
	  IF (BC_TYPE.EQ.2) WRITE (*,*) '    BOUNDARY CONDITION                 : SOLID WALL'
	  IF (BC_TYPE.EQ.3) WRITE (*,*) '    BOUNDARY CONDITION                 : FACTS INPUT'
	  IF (ABS(LO%TIDE_LEVEL).GT.GX) WRITE (*,*) '    TIDAL LEVEL CORRECTION (METERS):',LO%TIDE_LEVEL
	  WRITE (*,*) ' '

      WRITE (*,*) '---------------- INITIAL CONDITION INFORMATION ---------------'
	  IF (INI_SURF .EQ. 0 .OR. INI_SURF .EQ. 9) THEN
	      WRITE (*,*) '    #USE BUILT-IN FAULT MODEL'
		  WRITE (*,*) '       TOTAL NUMBER OF FAULT SEGMENTS  :',FAULT_INFO(1)%NUM_FLT
		  DO K = 1,FAULT_INFO(1)%NUM_FLT
		     IF (FAULT_INFO(K)%SWITCH.EQ.0) THEN
			    WRITE (*,*) '     PARAMETERS FOR FAULT SEGMENT      :',K
		        WRITE (*,*) '       FAULT RUPTURE TIME     (SECOND) :',FAULT_INFO(K)%T0
		        WRITE (*,*) '       EPICENTER (LON, LAT)   (DEGREE) :',FAULT_INFO(K)%X0,FAULT_INFO(K)%Y0
		        WRITE (*,*) '       FOCAL DEPTH         (KILOMETER) :',FAULT_INFO(K)%HH/1000.0
		        WRITE (*,*) '       FAULT LENGTH        (KILOMETER) :',FAULT_INFO(K)%L/1000.0
		        WRITE (*,*) '       FAULT WIDTH         (KILOMETER) :',FAULT_INFO(K)%W/1000.0
		        WRITE (*,*) '       STRIKE ANGLE    (THETA, DEGREE) :',FAULT_INFO(K)%TH
		        WRITE (*,*) '       DIP ANGLE       (DELTA, DEGREE) :',FAULT_INFO(K)%DL
		        WRITE (*,*) '       SLIP ANGLE     (LAMBDA, DEGREE) :',FAULT_INFO(K)%RD
		        WRITE (*,*) '       DISLOCATION             (METER) :',FAULT_INFO(K)%D
			 ELSE
			    WRITE (*,*) '    USE DEFORMATION DATA FOR SEGMENT   :',K
		        WRITE (*,*) '       DEFORMATION DATA FILE NAME      :',FAULT_INFO(K)%DEFORM_NAME
			 ENDIF
		  ENDDO
	  ELSE IF (INI_SURF .EQ. 1) THEN
	      WRITE (*,*) '    #USE INITIAL SURFACE DEFORMATION FILE:',FAULT_INFO(1)%DEFORM_NAME
	  ELSE IF (INI_SURF.EQ.2) THEN
	      WRITE (*,*) '    #USE INCIDENT WAVE MODEL'
          IF (WAVE_INFO%MK_TYPE==1) WRITE (*,*) '         INCIDENT SOLITARY WAVE'
		  IF (WAVE_INFO%MK_TYPE==1) WRITE (*,*) '         CHARACTERISTIC WATER DEPTH (M):',WAVE_INFO%DEPTH
	      IF (WAVE_INFO%MK_TYPE==1) WRITE (*,*) '         CHARACTERISTIC WAVE HEIGHT (M):',WAVE_INFO%AMP
          IF (WAVE_INFO%MK_TYPE==2) WRITE (*,*) '         CUSTOMIZED TIMEHISTORY PROFILE:',WAVE_INFO%FORM_NAME
		  IF (WAVE_INFO%MK_TYPE==2) WRITE (*,*) '         TOTAL ENTRIES IN INPUT FILE   :',WAVE_INFO%FORM_LEN
		  IF (WAVE_INFO%MK_TYPE==2) WRITE (*,*) '         ENDING TIME OF WAVE INPUT     :',WAVE_INFO%T(WAVE_INFO%FORM_LEN)
	  ELSE IF (INI_SURF .EQ. 3) THEN
	      WRITE (*,*) '    #USE SUBMARINE LANDSLIDE MODEL'
		  IF (LANDSLIDE_INFO%OPTION .LE. 1) THEN
		  WRITE (*,*) '       OBTAIN SNAPSHOTS OF WATERDEPTH VARIATIONS FROM INPUT DATA FILE'
		  WRITE (*,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (*,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
          WRITE (*,*) '       X_START                         :',LANDSLIDE_INFO%X_START
          WRITE (*,*) '       X_END                           :',LANDSLIDE_INFO%X_END
          WRITE (*,*) '       Y_START                         :',LANDSLIDE_INFO%Y_START
          WRITE (*,*) '       Y_END                           :',LANDSLIDE_INFO%Y_END
		  WRITE (*,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (*,*) '       DATA FILE NAME OF DEPTH VARIATION SNAPSHOTS:',LANDSLIDE_INFO%FILENAME
		  ENDIF
		  IF (LANDSLIDE_INFO%OPTION .EQ. 2) THEN
		  WRITE (*,*) '       CALC. SNAPSHOTS OF WATERDEPTH VARIATIONS FROM BUILT-IN THEORY'
          WRITE (*,*) '       X START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%X_START
          WRITE (*,*) '       X END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%X_END
          WRITE (*,*) '       Y START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%Y_START
          WRITE (*,*) '       Y END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%Y_END
		  WRITE (*,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (*,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(1)
		  WRITE (*,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(LANDSLIDE_INFO%NT)
		  WRITE (*,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (*,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
		  WRITE (*,*) '       LANDSLIDE STARTS AT (X,Y)       :',LANDSLIDE_INFO%XS,LANDSLIDE_INFO%YS
		  WRITE (*,*) '       LANDSLIDE STOPS AT (X,Y)        :',LANDSLIDE_INFO%XE,LANDSLIDE_INFO%YE
          WRITE (*,*) '       TYPICAL SLOPE OF PATH   (DEGREE):',LANDSLIDE_INFO%SLOPE
          WRITE (*,*) '       LENGTH OF SLIDING MASS   (METER):',LANDSLIDE_INFO%A
          WRITE (*,*) '       WIDTH OF SLIDING MASS    (METER):',LANDSLIDE_INFO%B
          WRITE (*,*) '       THICKNESS OF SLIDING MASS(METER):',LANDSLIDE_INFO%THICKNESS
		  ENDIF
	  ELSE IF (INI_SURF .EQ. 4) THEN
	      WRITE (*,*) '    #USE FAULT MODEL + LANDSLIDE MODEL'
		  WRITE (*,*) '       TOTAL NUMBER OF FAULT SEGMENTS  :',FAULT_INFO(1)%NUM_FLT
		  DO K = 1,FAULT_INFO(1)%NUM_FLT
		     IF (FAULT_INFO(K)%SWITCH.EQ.0) THEN
			    WRITE (*,*) '     PARAMETERS FOR FAULT SEGMENT      :',K
		        WRITE (*,*) '       FAULT RUPTURE TIME     (SECOND) :',FAULT_INFO(K)%T0
		        WRITE (*,*) '       EPICENTER (LON, LAT)   (DEGREE) :',FAULT_INFO(K)%X0,FAULT_INFO(K)%Y0
		        WRITE (*,*) '       FOCAL DEPTH         (KILOMETER) :',FAULT_INFO(K)%HH/1000.0
		        WRITE (*,*) '       FAULT LENGTH        (KILOMETER) :',FAULT_INFO(K)%L/1000.0
		        WRITE (*,*) '       FAULT WIDTH         (KILOMETER) :',FAULT_INFO(K)%W/1000.0
		        WRITE (*,*) '       STRIKE ANGLE    (THETA, DEGREE) :',FAULT_INFO(K)%TH
		        WRITE (*,*) '       DIP ANGLE       (DELTA, DEGREE) :',FAULT_INFO(K)%DL
		        WRITE (*,*) '       SLIP ANGLE     (LAMBDA, DEGREE) :',FAULT_INFO(K)%RD
		        WRITE (*,*) '       DISLOCATION             (METER) :',FAULT_INFO(K)%D
			 ELSE
			    WRITE (*,*) '     USE DEFORMATION DATA FOR SEGMENT   :',K
		        WRITE (*,*) '       DEFORMATION DATA FILE NAME      :',FAULT_INFO(K)%DEFORM_NAME
			 ENDIF
		  ENDDO
	      WRITE (*,*) '       PARAMETERS FOR LANDSLIDE MODEL  :'
		  IF (LANDSLIDE_INFO%OPTION .LE. 1) THEN
		  WRITE (*,*) '       OBTAIN SNAPSHOTS OF WATERDEPTH VARIATIONS FROM INPUT DATA FILE'
		  WRITE (*,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (*,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
          WRITE (*,*) '       X_START                         :',LANDSLIDE_INFO%X_START
          WRITE (*,*) '       X_END                           :',LANDSLIDE_INFO%X_END
          WRITE (*,*) '       Y_START                         :',LANDSLIDE_INFO%Y_START
          WRITE (*,*) '       Y_END                           :',LANDSLIDE_INFO%Y_END
		  WRITE (*,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (*,*) '       DATA FILE NAME OF DEPTH VARIATION SNAPSHOTS:',LANDSLIDE_INFO%FILENAME
		  ENDIF
		  IF (LANDSLIDE_INFO%OPTION .EQ. 2) THEN
		  WRITE (*,*) '       CALC. SNAPSHOTS OF WATERDEPTH VARIATIONS FROM BUILT-IN THEORY'
          WRITE (*,*) '       X START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%X_START
          WRITE (*,*) '       X END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%X_END
          WRITE (*,*) '       Y START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%Y_START
          WRITE (*,*) '       Y END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%Y_END
		  WRITE (*,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (*,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(1)
		  WRITE (*,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(LANDSLIDE_INFO%NT)
		  WRITE (*,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (*,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
		  WRITE (*,*) '       LANDSLIDE STARTS AT (X,Y)       :',LANDSLIDE_INFO%XS,LANDSLIDE_INFO%YS
		  WRITE (*,*) '       LANDSLIDE STOPS AT (X,Y)        :',LANDSLIDE_INFO%XE,LANDSLIDE_INFO%YE
          WRITE (*,*) '       TYPICAL SLOPE OF PATH   (DEGREE):',LANDSLIDE_INFO%SLOPE
          WRITE (*,*) '       LENGTH OF SLIDING MASS   (METER):',LANDSLIDE_INFO%A
          WRITE (*,*) '       WIDTH OF SLIDING MASS    (METER):',LANDSLIDE_INFO%B
          WRITE (*,*) '       THICKNESS OF SLIDING MASS(METER):',LANDSLIDE_INFO%THICKNESS
		  ENDIF

	  ENDIF
	  IF (BC_TYPE.EQ.3) THEN
		  WRITE (*,*) '    USE FACTS INPUT BOUNDARY CONDITION!'
		  WRITE (*,*) '       GRID DIMENSION OF FACTS INPUT   :',BCI_INFO%NX,'*',BCI_INFO%NY
		  WRITE (*,*) '       TOTAL SNAPSHOTS OF FACTS INPUT  :',BCI_INFO%NT
		  WRITE (*,*) '       STARTING TIME OF FACTS INPUT    :',BCI_INFO%T(1)
		  WRITE (*,*) '       ENDING TIME OF FACTS INPUT      :',BCI_INFO%T(BCI_INFO%NT)
		  WRITE (*,*) '       DURATION OF FACTS INPUT         :',BCI_INFO%DURATION
		  WRITE (*,*) '       X_START OF FACTS INPUT          :',BCI_INFO%X(1)
		  WRITE (*,*) '       X_END OF FACTS INPUT            :',BCI_INFO%X(BCI_INFO%NX)
		  WRITE (*,*) '       Y_START OF FACTS INPUT          :',BCI_INFO%Y(1)
		  WRITE (*,*) '       Y_END OF FACTS INPUT            :',BCI_INFO%Y(BCI_INFO%NY)
		  WRITE (*,*) '       FILENAMES OF FACTS INPUT        :'
		  WRITE (*,*) BCI_INFO%FNAMEH
		  WRITE (*,*) BCI_INFO%FNAMEU
		  WRITE (*,*) BCI_INFO%FNAMEV
	  ENDIF
	  WRITE (*,*) ' '
	  WRITE (*,*) '--------------- 1ST-LEVEL GRID INFORMATION -----------------'
	  WRITE (*,*) '    #  GRID IDENTIFICATION NUMBER      :',LO%ID
	  IF (LO%LAYCORD .EQ. 0) WRITE (*,*) '       USE SPHERICAL COORDINATE SYSTEM'
      IF (LO%LAYCORD .EQ. 1) WRITE (*,*) '       USE CARTESIAN COORDINATE SYSTEM'
	  IF (LO%LAYGOV .EQ. 0)  WRITE (*,*) '       USE LINEAR SHALLOW WATER EQUATIONS'
      IF (LO%LAYGOV .EQ. 1)  WRITE (*,*) '       USE NONLINEAR SHALLOW WATER EQUATIONS'
	  IF (LO%LAYGOV .EQ. 2)  WRITE (*,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
	  IF (LO%LAYGOV .EQ. 3)  WRITE (*,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
	  IF (LO%LAYGOV .EQ. 4)  WRITE (*,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
	  IF (LO%LAYGOV .EQ. 5)  WRITE (*,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
	  IF (LO%FRIC_SWITCH .EQ. 1) WRITE (*,*) '       BOTTOM FRICTION                 : DISABLED'
	  IF (LO%FRIC_SWITCH .NE. 1) WRITE (*,*) '       BOTTOM FRICTION                 : ENABLED'
	  IF (LO%FRIC_SWITCH .EQ. 0) WRITE (*,*) '          USE CONSTANT ROUGHNESS COEF. :',LO%FRIC_COEF
	  IF (LO%FRIC_SWITCH .EQ. 2) WRITE (*,*) '          USE VARIABLE ROUGHNESS COEFFICIENTS.'
	  IF (LO%FLUXSWITCH .EQ. 0) WRITE (*,*)  '       VOLUME FLUX OUTPUT              : ENABLED'
	  IF (LO%FLUXSWITCH .EQ. 1) WRITE (*,*)  '       VOLUME FLUX OUTPUT              : DISABLED'
	  IF (LO%FLUXSWITCH .NE. 2) WRITE (*,*)  '       SURFACE DISPLACEMENT OUTPUT     : ENABLED'
	  IF (LO%FLUXSWITCH .EQ. 2) WRITE (*,*)  '       SURFACE DISPLACEMENT OUTPUT     : DISABLED'
	  WRITE (*,*) '       GRID LAYER POSITIONS            :'
	  IF (LO%LAYCORD .EQ. 0) THEN
	     WRITE (*,*) '         X_START              (DEGREE) :',LO%X(1)
	     WRITE (*,*) '         X_END                (DEGREE) :',LO%X(LO%NX)
	     WRITE (*,*) '         Y_START              (DEGREE) :',LO%Y(1)
	     WRITE (*,*) '         Y_END                (DEGREE) :',LO%Y(LO%NY)
		 WRITE (*,*) '       GRID DIMENSION          (NX*NY) :',LO%NX,'*',LO%NY
	     WRITE (*,*) '       GRID SIZE           (DX,MINUTE) :',LO%DX
	     WRITE (*,*) '       GRID SIZE           (DY,MINUTE) :',LO%DY
		 WRITE (*,*) '       TIME STEP SIZE     (DT, SECOND) :',LO%DT
	  ENDIF
	  IF (LO%LAYCORD .EQ. 1) THEN
	     WRITE (*,*) '         X_START               (METER) :',LO%X(1)
	     WRITE (*,*) '         X_END                 (METER) :',LO%X(LO%NX)
	     WRITE (*,*) '         Y_START               (METER) :',LO%Y(1)
	     WRITE (*,*) '         Y_END                 (METER) :',LO%Y(LO%NY)
		 WRITE (*,*) '       GRID DIMENSION          (NX*NY) :',LO%NX,'*',LO%NY
	     WRITE (*,*) '       GRID SIZE           (DX, METER) :',LO%DX
	     WRITE (*,*) '       GRID SIZE           (DY, METER) :',LO%DY
		 WRITE (*,*) '       TIME STEP SIZE     (DT, SECOND) :',LO%DT
	  ENDIF
	  WRITE (*,*) '       NUMBER OF CHILD GRID LAYERS     :',LO%NUM_CHILD
	  DO I = 1,NUM_GRID
	  IF (LA(I)%LAYSWITCH.EQ.0 .AND. LA(I)%PARENT.EQ.LO%ID) THEN
	  WRITE (*,*) '               CHILD GRID LAYER ID     :',LA(I)%ID
	  ENDIF
	  ENDDO
	  WRITE (*,*) '       BATHYMETRY DATA FILE NAME       :',LO%DEPTH_NAME

!*	  CALL CR_CHECK (LO,LA)
  	  DO I=1,NUM_GRID
	     IF (LA(I)%LAYSWITCH .EQ. 0) THEN
            WRITE (*,*) '--------------- SUB-LEVEL GRID INFORMATION -----------------'
		    WRITE (*,*) '    #  GRID IDENTIFICATION NUMBER      :',LA(I)%ID
	        IF (LA(I)%LAYCORD .EQ. 0) WRITE (*,*) '       USE SPHERICAL COORDINATE SYSTEM'
            IF (LA(I)%LAYCORD .EQ. 1) WRITE (*,*) '       USE CARTESIAN COORDINATE SYSTEM'
	        IF (LA(I)%LAYGOV .EQ. 0)  WRITE (*,*) '       USE LINEAR SHALLOW WATER EQUATIONS'
            IF (LA(I)%LAYGOV .EQ. 1)  WRITE (*,*) '       USE NONLINEAR SHALLOW WATER EQUATIONS'
			IF (LA(I)%LAYGOV .EQ. 2)  WRITE (*,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
            IF (LA(I)%LAYGOV .EQ. 3)  WRITE (*,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
			IF (LA(I)%LAYGOV .EQ. 4)  WRITE (*,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
            IF (LA(I)%LAYGOV .EQ. 5)  WRITE (*,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
	        IF (LA(I)%FRIC_SWITCH .EQ. 1) WRITE (*,*) '       BOTTOM FRICTION                 : DISABLED'
			IF (LA(I)%FRIC_SWITCH .NE. 1) WRITE (*,*) '       BOTTOM FRICTION                 : ENABLED'
			IF (LA(I)%FRIC_SWITCH .EQ. 0) WRITE (*,*) '          USE CONSTANT ROUGHNESS COEF. :',LO%FRIC_COEF
			IF (LA(I)%FRIC_SWITCH .EQ. 2) WRITE (*,*) '          USE VARIABLE MANNING ROUGHNESS COEFFICIENTS.'
			IF (LA(I)%FLUXSWITCH .EQ. 0) WRITE (*,*)  '       VOLUME FLUX OUTPUT              : ENABLED'
			IF (LA(I)%FLUXSWITCH .EQ. 1) WRITE (*,*)  '       VOLUME FLUX OUTPUT              : DISABLED'
			IF (LA(I)%FLUXSWITCH .NE. 2) WRITE (*,*)  '       SURFACE DISPLACEMENT OUTPUT     : ENABLED'
			IF (LA(I)%FLUXSWITCH .EQ. 2) WRITE (*,*)  '       SURFACE DISPLACEMENT OUTPUT     : DISABLED'
			WRITE (*,*) '       PARENT GRID ID                  :',LA(I)%PARENT
			WRITE (*,*) '       GRID SIZE RATIO                 :',LA(I)%REL_SIZE
			WRITE (*,*) '       TIME STEP SIZE RATIO            :',LA(I)%REL_TIME
			WRITE (*,*) '       POSITIONS IN ITS PARENT LAYER   :'
			IF (LA(I)%LAYCORD .EQ. 0) THEN
	           WRITE (*,*) '         X_START              (DEGREE) :',LA(I)%X(2)
	           WRITE (*,*) '         X_END                (DEGREE) :',LA(I)%X(LA(I)%NX)
	           WRITE (*,*) '         Y_START              (DEGREE) :',LA(I)%Y(2)
	           WRITE (*,*) '         Y_END                (DEGREE) :',LA(I)%Y(LA(I)%NY)
	           WRITE (*,*) '         I_START               (INDEX) :',LA(I)%CORNERS(1)
	           WRITE (*,*) '         I_END                 (INDEX) :',LA(I)%CORNERS(2)
	           WRITE (*,*) '         J_START               (INDEX) :',LA(I)%CORNERS(3)
	           WRITE (*,*) '         J_END                 (INDEX) :',LA(I)%CORNERS(4)
			   WRITE (*,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
			   WRITE (*,*) '       GRID SIZE           (DX,MINUTE) :',LA(I)%DX
			   WRITE (*,*) '       GRID SIZE           (DY,MINUTE) :',LA(I)%DY
			   WRITE (*,*) '       TIME STEP SIZE     (DT, SECOND) :',LA(I)%DT
			ENDIF
			IF (LA(I)%LAYCORD .EQ. 1) THEN
			   IF (LA(I)%UPZ) THEN
	           WRITE (*,*) '         X_START               (METER) :',LA(I)%X(2)
			   WRITE (*,*) '         X_END                 (METER) :',LA(I)%X(LA(I)%NX)
			   WRITE (*,*) '         Y_START               (METER) :',LA(I)%Y(2)
			   WRITE (*,*) '         Y_END                 (METER) :',LA(I)%Y(LA(I)%NY)
	           WRITE (*,*) '         I_START               (INDEX) :',LA(I)%CORNERS(1)
	           WRITE (*,*) '         I_END                 (INDEX) :',LA(I)%CORNERS(2)
	           WRITE (*,*) '         J_START               (INDEX) :',LA(I)%CORNERS(3)
	           WRITE (*,*) '         J_END                 (INDEX) :',LA(I)%CORNERS(4)
			   WRITE (*,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
			   ELSE
			   WRITE (*,*) '         X_START_S            (DEGREE) :',LA(I)%XT(2)
	           WRITE (*,*) '         X_END_S              (DEGREE) :',LA(I)%XT(LA(I)%NX)
			   WRITE (*,*) '         Y_START_S            (DEGREE) :',LA(I)%YT(2)
	           WRITE (*,*) '         Y_END_S              (DEGREE) :',LA(I)%YT(LA(I)%NY)
	           WRITE (*,*) '         X_START_C             (METER) :',LA(I)%X(2)
	           WRITE (*,*) '         X_END_C               (METER) :',LA(I)%X(LA(I)%NX)
	           WRITE (*,*) '         Y_START_C             (METER) :',LA(I)%Y(2)
	           WRITE (*,*) '         Y_END_C               (METER) :',LA(I)%Y(LA(I)%NY)
			   WRITE (*,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
			   ENDIF
			   WRITE (*,*) '       GRID SIZE           (DX, METER) :',LA(I)%DX
			   WRITE (*,*) '       GRID SIZE           (DY, METER) :',LA(I)%DY
			   WRITE (*,*) '       TIME STEP SIZE     (DT, SECOND) :',LA(I)%DT
			ENDIF
			WRITE (*,*) '       NUMBER OF CHILD GRID LAYERS     :',LA(I)%NUM_CHILD
			DO K = 1,NUM_GRID
			IF (LA(K)%LAYSWITCH.EQ.0 .AND. LA(K)%PARENT.EQ.LA(I)%ID) THEN
			WRITE (*,*) '               CHILD GRID LAYER ID     :',LA(K)%ID
			ENDIF
			ENDDO
			WRITE (*,*) '       BATHYMETRY DATA FILE NAME       :',LA(I)%DEPTH_NAME
			WRITE (*,*) '------------------------------------------------------------'
		 ENDIF
	  ENDDO
	  WRITE (*,*) ' '

!.....SAVE SIMULATION INFORMATION INTO DATA FILE: INPUT_INFO.TXT
	  OPEN (UNIT=23,FILE='SIMULATION_INFO.TXT',STATUS='UNKNOWN')

	  WRITE (23,*) '***********************************************************'
	  WRITE (23,*) '*           INPUT INFORMATION - COMCOT V1.7               *'
	  WRITE (23,*) '***********************************************************'
	  WRITE (23,*) JOB
      WRITE (23,*) '------------------- GENERAL INFORMATION -------------------'
 	  WRITE (23,*) '    TOTAL RUN TIME            (SECOND) :',TEND
	  WRITE (23,*) '    TIME INTERVAL FOR OUTPUT  (SECOND) :',T_INV
	  IF (T_INV .LT. LO%DT) THEN
	     T_INV = LO%DT
		 IPRT = NINT(T_INV/LO%DT)
		 WRITE (23,*) '    OUTPUT INTERVAL < DT, RESET TO', T_INV
		 WRITE (23,*) '   '
      ENDIF
	  WRITE (23,*) '    TIME STEP SIZE            (SECOND) :',LO%DT
	  WRITE (23,*) '    TOTAL STEPS TO RUN         (STEPS) :',KEND
	  WRITE (23,*) '    STEP INTERVAL FOR OUTPUT   (STEPS) :',IPRT
	  IF (OUT_OPT.EQ.0 .OR. OUT_OPT.EQ.2) WRITE (23,*) '    MAX. SURFACE DISPLACEMENT OUTPUT   : ENABLED'
	  IF (OUT_OPT.EQ.1 .OR. OUT_OPT.EQ.2) WRITE (23,*) '    TIME HISTORY RECORD OUTPUT         : ENABLED'
	  IF (OUT_OPT .EQ. 1) WRITE (23,*) '    MAX. SURFACE DISPLACEMENT OUTPUT   : DISABLED'
	  IF (OUT_OPT .EQ. 0) WRITE (23,*) '    TIME HISTORY RECORD OUTPUT         : DISABLED'
	  IF (START_TYPE.EQ.1) THEN
	     WRITE (23,*) '    HOT START FUNCTION                 : ENABLED'
		 WRITE (23,*) '       SIMULATION RESUMES FROM T =     :',START_TIME
		 ISTEP_LEFT = NINT((TEND-START_TIME)/LO%DT)
		 WRITE (23,*) '       TOTAL NUMBER OF TIME STEPS LEFT :',ISTEP_LEFT
      ENDIF
	  WRITE (23,*) '    SHORELINE LOCATED AT DEPTH CONTOUR :',LO%H_LIMIT
	  IF (BC_TYPE.EQ.0) WRITE (23,*) '    BOUNDARY CONDITION                 : RADIATION (OPEN)'
	  IF (BC_TYPE.EQ.1) WRITE (23,*) '    BOUNDARY CONDITION                 : ABSORBING'
	  IF (BC_TYPE.EQ.2) WRITE (23,*) '    BOUNDARY CONDITION                 : SOLID WALL'
	  IF (BC_TYPE.EQ.3) WRITE (23,*) '    BOUNDARY CONDITION                 : FACTS INPUT'
	  IF (ABS(LO%TIDE_LEVEL).GT.GX) WRITE (23,*) '    TIDAL LEVEL CORRECTION (METERS):',LO%TIDE_LEVEL
	  WRITE (23,*) ' '

      WRITE (23,*) '---------------- INITIAL CONDITION INFORMATION ---------------'
	  IF (INI_SURF .EQ. 0 .OR. INI_SURF .EQ. 9) THEN
	      WRITE (23,*) '    #USE BUILT-IN FAULT MODEL'
		  WRITE (23,*) '       TOTAL NUMBER OF FAULT SEGMENTS  :',FAULT_INFO(1)%NUM_FLT
		  DO K = 1,FAULT_INFO(1)%NUM_FLT
		     IF (FAULT_INFO(K)%SWITCH.EQ.0) THEN
			    WRITE (23,*) '     PARAMETERS FOR FAULT SEGMENT      :',K
		        WRITE (23,*) '       FAULT RUPTURE TIME     (SECOND) :',FAULT_INFO(K)%T0
		        WRITE (23,*) '       EPICENTER (LON, LAT)   (DEGREE) :',FAULT_INFO(K)%X0,FAULT_INFO(K)%Y0
		        WRITE (23,*) '       FOCAL DEPTH         (KILOMETER) :',FAULT_INFO(K)%HH/1000.0
		        WRITE (23,*) '       FAULT LENGTH        (KILOMETER) :',FAULT_INFO(K)%L/1000.0
		        WRITE (23,*) '       FAULT WIDTH         (KILOMETER) :',FAULT_INFO(K)%W/1000.0
		        WRITE (23,*) '       STRIKE ANGLE    (THETA, DEGREE) :',FAULT_INFO(K)%TH
		        WRITE (23,*) '       DIP ANGLE       (DELTA, DEGREE) :',FAULT_INFO(K)%DL
		        WRITE (23,*) '       SLIP ANGLE     (LAMBDA, DEGREE) :',FAULT_INFO(K)%RD
		        WRITE (23,*) '       DISLOCATION             (METER) :',FAULT_INFO(K)%D
			 ELSE
			    WRITE (23,*) '    USE DEFORMATION DATA FOR SEGMENT   :',K
		        WRITE (23,*) '       DEFORMATION DATA FILE NAME      :',FAULT_INFO(K)%DEFORM_NAME
			 ENDIF
		  ENDDO
	  ELSE IF (INI_SURF .EQ. 1) THEN
	      WRITE (23,*) '    #USE INITIAL SURFACE DEFORMATION FILE:',FAULT_INFO(1)%DEFORM_NAME
	  ELSE IF (INI_SURF.EQ.2) THEN
	      WRITE (23,*) '    #USE INCIDENT WAVE MODEL'
          IF (WAVE_INFO%MK_TYPE==1) WRITE (23,*) '         INCIDENT SOLITARY WAVE'
		  IF (WAVE_INFO%MK_TYPE==1) WRITE (23,*) '         CHARACTERISTIC WATER DEPTH (M):',WAVE_INFO%DEPTH
	      IF (WAVE_INFO%MK_TYPE==1) WRITE (23,*) '         CHARACTERISTIC WAVE HEIGHT (M):',WAVE_INFO%AMP
          IF (WAVE_INFO%MK_TYPE==2) WRITE (23,*) '         CUSTOMIZED TIMEHISTORY PROFILE:',WAVE_INFO%FORM_NAME
		  IF (WAVE_INFO%MK_TYPE==2) WRITE (23,*) '         TOTAL ENTRIES IN INPUT FILE   :',WAVE_INFO%FORM_LEN
		  IF (WAVE_INFO%MK_TYPE==2) WRITE (23,*) '         ENDING TIME OF WAVE INPUT     :',WAVE_INFO%T(WAVE_INFO%FORM_LEN)
	  ELSE IF (INI_SURF .EQ. 3) THEN
	      WRITE (23,*) '    #USE SUBMARINE LANDSLIDE MODEL'
		  IF (LANDSLIDE_INFO%OPTION .LE. 1) THEN
		  WRITE (23,*) '       OBTAIN SNAPSHOTS OF WATERDEPTH VARIATIONS FROM INPUT DATA FILE'
		  WRITE (23,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (23,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
          WRITE (23,*) '       X_START                         :',LANDSLIDE_INFO%X_START
          WRITE (23,*) '       X_END                           :',LANDSLIDE_INFO%X_END
          WRITE (23,*) '       Y_START                         :',LANDSLIDE_INFO%Y_START
          WRITE (23,*) '       Y_END                           :',LANDSLIDE_INFO%Y_END
		  WRITE (23,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (23,*) '       DATA FILE NAME OF DEPTH VARIATION SNAPSHOTS:',LANDSLIDE_INFO%FILENAME
		  ENDIF
		  IF (LANDSLIDE_INFO%OPTION .EQ. 2) THEN
		  WRITE (23,*) '       CALC. SNAPSHOTS OF WATERDEPTH VARIATIONS FROM BUILT-IN THEORY'
          WRITE (23,*) '       X START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%X_START
          WRITE (23,*) '       X END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%X_END
          WRITE (23,*) '       Y START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%Y_START
          WRITE (23,*) '       Y END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%Y_END
		  WRITE (23,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (23,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(1)
		  WRITE (23,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(LANDSLIDE_INFO%NT)
		  WRITE (23,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (23,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
		  WRITE (23,*) '       LANDSLIDE STARTS AT (X,Y)       :',LANDSLIDE_INFO%XS,LANDSLIDE_INFO%YS
		  WRITE (23,*) '       LANDSLIDE STOPS AT (X,Y)        :',LANDSLIDE_INFO%XE,LANDSLIDE_INFO%YE
          WRITE (23,*) '       TYPICAL SLOPE OF PATH   (DEGREE):',LANDSLIDE_INFO%SLOPE
          WRITE (23,*) '       LENGTH OF SLIDING MASS   (METER):',LANDSLIDE_INFO%A
          WRITE (23,*) '       WIDTH OF SLIDING MASS    (METER):',LANDSLIDE_INFO%B
          WRITE (23,*) '       THICKNESS OF SLIDING MASS(METER):',LANDSLIDE_INFO%THICKNESS
		  ENDIF
	  ELSE IF (INI_SURF .EQ. 4) THEN
	      WRITE (23,*) '    #USE FAULT MODEL + LANDSLIDE MODEL'
		  WRITE (23,*) '       TOTAL NUMBER OF FAULT SEGMENTS  :',FAULT_INFO(1)%NUM_FLT
		  DO K = 1,FAULT_INFO(1)%NUM_FLT
		     IF (FAULT_INFO(K)%SWITCH.EQ.0) THEN
			    WRITE (23,*) '     PARAMETERS FOR FAULT SEGMENT      :',K
		        WRITE (23,*) '       FAULT RUPTURE TIME     (SECOND) :',FAULT_INFO(K)%T0
		        WRITE (23,*) '       EPICENTER (LON, LAT)   (DEGREE) :',FAULT_INFO(K)%X0,FAULT_INFO(K)%Y0
		        WRITE (23,*) '       FOCAL DEPTH         (KILOMETER) :',FAULT_INFO(K)%HH/1000.0
		        WRITE (23,*) '       FAULT LENGTH        (KILOMETER) :',FAULT_INFO(K)%L/1000.0
		        WRITE (23,*) '       FAULT WIDTH         (KILOMETER) :',FAULT_INFO(K)%W/1000.0
		        WRITE (23,*) '       STRIKE ANGLE    (THETA, DEGREE) :',FAULT_INFO(K)%TH
		        WRITE (23,*) '       DIP ANGLE       (DELTA, DEGREE) :',FAULT_INFO(K)%DL
		        WRITE (23,*) '       SLIP ANGLE     (LAMBDA, DEGREE) :',FAULT_INFO(K)%RD
		        WRITE (23,*) '       DISLOCATION             (METER) :',FAULT_INFO(K)%D
			 ELSE
			    WRITE (23,*) '     USE DEFORMATION DATA FOR SEGMENT   :',K
		        WRITE (23,*) '       DEFORMATION DATA FILE NAME      :',FAULT_INFO(K)%DEFORM_NAME
			 ENDIF
		  ENDDO
	      WRITE (23,*) '       PARAMETERS FOR LANDSLIDE MODEL  :'
		  IF (LANDSLIDE_INFO%OPTION .LE. 1) THEN
		  WRITE (23,*) '       OBTAIN SNAPSHOTS OF WATERDEPTH VARIATIONS FROM INPUT DATA FILE'
		  WRITE (23,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (23,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
          WRITE (23,*) '       X_START                         :',LANDSLIDE_INFO%X_START
          WRITE (23,*) '       X_END                           :',LANDSLIDE_INFO%X_END
          WRITE (23,*) '       Y_START                         :',LANDSLIDE_INFO%Y_START
          WRITE (23,*) '       Y_END                           :',LANDSLIDE_INFO%Y_END
		  WRITE (23,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (23,*) '       DATA FILE NAME OF DEPTH VARIATION SNAPSHOTS:',LANDSLIDE_INFO%FILENAME
		  ENDIF
		  IF (LANDSLIDE_INFO%OPTION .EQ. 2) THEN
		  WRITE (23,*) '       CALC. SNAPSHOTS OF WATERDEPTH VARIATIONS FROM BUILT-IN THEORY'
          WRITE (23,*) '       X START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%X_START
          WRITE (23,*) '       X END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%X_END
          WRITE (23,*) '       Y START OF LANDSLIDE AREA       :',LANDSLIDE_INFO%Y_START
          WRITE (23,*) '       Y END OF LANDSLIDE AREA         :',LANDSLIDE_INFO%Y_END
		  WRITE (23,*) '       GRID DIMENSION          (NX,NY) :',LANDSLIDE_INFO%NX,LANDSLIDE_INFO%NY
		  WRITE (23,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(1)
		  WRITE (23,*) '       LANDSLIDE STARTING TIME(SECOND) :',LANDSLIDE_INFO%T(LANDSLIDE_INFO%NT)
		  WRITE (23,*) '       LANDSLIDE DURATION     (SECOND) :',LANDSLIDE_INFO%DURATION
          WRITE (23,*) '       TOTAL NUMBER OF SNAPSHOTS       :',LANDSLIDE_INFO%NT
		  WRITE (23,*) '       LANDSLIDE STARTS AT (X,Y)       :',LANDSLIDE_INFO%XS,LANDSLIDE_INFO%YS
		  WRITE (23,*) '       LANDSLIDE STOPS AT (X,Y)        :',LANDSLIDE_INFO%XE,LANDSLIDE_INFO%YE
          WRITE (23,*) '       TYPICAL SLOPE OF PATH   (DEGREE):',LANDSLIDE_INFO%SLOPE
          WRITE (23,*) '       LENGTH OF SLIDING MASS   (METER):',LANDSLIDE_INFO%A
          WRITE (23,*) '       WIDTH OF SLIDING MASS    (METER):',LANDSLIDE_INFO%B
          WRITE (23,*) '       THICKNESS OF SLIDING MASS(METER):',LANDSLIDE_INFO%THICKNESS
		  ENDIF
	  ENDIF
	  IF (BC_TYPE.EQ.3) THEN
		  WRITE (23,*) '    USE FACTS INPUT BOUNDARY CONDITION!'
		  WRITE (23,*) '       GRID DIMENSION OF FACTS INPUT   :',BCI_INFO%NX,'*',BCI_INFO%NY
		  WRITE (23,*) '       TOTAL SNAPSHOTS OF FACTS INPUT  :',BCI_INFO%NT
		  WRITE (23,*) '       STARTING TIME OF FACTS INPUT    :',BCI_INFO%T(1)
		  WRITE (23,*) '       ENDING TIME OF FACTS INPUT      :',BCI_INFO%T(BCI_INFO%NT)
		  WRITE (23,*) '       DURATION OF FACTS INPUT         :',BCI_INFO%DURATION
		  WRITE (23,*) '       X_START OF FACTS INPUT          :',BCI_INFO%X(1)
		  WRITE (23,*) '       X_END OF FACTS INPUT            :',BCI_INFO%X(BCI_INFO%NX)
		  WRITE (23,*) '       Y_START OF FACTS INPUT          :',BCI_INFO%Y(1)
		  WRITE (23,*) '       Y_END OF FACTS INPUT            :',BCI_INFO%Y(BCI_INFO%NY)
		  WRITE (23,*) '       FILENAMES OF FACTS INPUT        :'
		  WRITE (23,*) BCI_INFO%FNAMEH
		  WRITE (23,*) BCI_INFO%FNAMEU
		  WRITE (23,*) BCI_INFO%FNAMEV
	  ENDIF
	  WRITE (23,*) ' '
	  WRITE (23,*) '--------------- 1ST-LEVEL GRID INFORMATION -----------------'
	  WRITE (23,*) '    #  GRID IDENTIFICATION NUMBER      :',LO%ID
	  IF (LO%LAYCORD .EQ. 0) WRITE (23,*) '       USE SPHERICAL COORDINATE SYSTEM'
      IF (LO%LAYCORD .EQ. 1) WRITE (23,*) '       USE CARTESIAN COORDINATE SYSTEM'
	  IF (LO%LAYGOV .EQ. 0)  WRITE (23,*) '       USE LINEAR SHALLOW WATER EQUATIONS'
      IF (LO%LAYGOV .EQ. 1)  WRITE (23,*) '       USE NONLINEAR SHALLOW WATER EQUATIONS'
	  IF (LO%LAYGOV .EQ. 2)  WRITE (23,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
      IF (LO%LAYGOV .EQ. 3)  WRITE (23,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
	  IF (LO%LAYGOV .EQ. 4)  WRITE (23,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
      IF (LO%LAYGOV .EQ. 5)  WRITE (23,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
	  IF (LO%FRIC_SWITCH .EQ. 1) WRITE (23,*) '       BOTTOM FRICTION                 : DISABLED'
	  IF (LO%FRIC_SWITCH .NE. 1) WRITE (23,*) '       BOTTOM FRICTION                 : ENABLED'
	  IF (LO%FRIC_SWITCH .EQ. 0) WRITE (23,*) '          USE CONSTANT ROUGHNESS COEF. :',LO%FRIC_COEF
	  IF (LO%FRIC_SWITCH .EQ. 2) WRITE (23,*) '          USE VARIABLE ROUGHNESS COEFFICIENTS.'
	  IF (LO%FLUXSWITCH .EQ. 0) WRITE (23,*)  '       VOLUME FLUX OUTPUT              : ENABLED'
	  IF (LO%FLUXSWITCH .EQ. 1) WRITE (23,*)  '       VOLUME FLUX OUTPUT              : DISABLED'
	  IF (LO%FLUXSWITCH .NE. 2) WRITE (23,*)  '       SURFACE DISPLACEMENT OUTPUT     : ENABLED'
	  IF (LO%FLUXSWITCH .EQ. 2) WRITE (23,*)  '       SURFACE DISPLACEMENT OUTPUT     : DISABLED'
	  WRITE (23,*) '       GRID LAYER POSITIONS            :'
	  IF (LO%LAYCORD .EQ. 0) THEN
	     WRITE (23,*) '           X_START            (DEGREE) :',LO%X(1)
	     WRITE (23,*) '           X_END              (DEGREE) :',LO%X(LO%NX)
	     WRITE (23,*) '           Y_START            (DEGREE) :',LO%Y(1)
	     WRITE (23,*) '           Y_END              (DEGREE) :',LO%Y(LO%NY)
		 WRITE (23,*) '       GRID DIMENSION          (NX*NY) :',LO%NX,'*',LO%NY
	     WRITE (23,*) '       GRID SIZE           (DX,MINUTE) :',LO%DX
	     WRITE (23,*) '       GRID SIZE           (DY,MINUTE) :',LO%DY
		 WRITE (23,*) '       TIME STEP SIZE     (DT, SECOND) :',LO%DT
	  ENDIF
	  IF (LO%LAYCORD .EQ. 1) THEN
	     WRITE (23,*) '         X_START               (METER) :',LO%X(1)
	     WRITE (23,*) '         X_END                 (METER) :',LO%X(LO%NX)
		 WRITE (23,*) '         Y_START               (METER) :',LO%Y(1)
		 WRITE (23,*) '         Y_END                 (METER) :',LO%Y(LO%NY)
		 WRITE (23,*) '       GRID DIMENSION          (NX*NY) :',LO%NX,'*',LO%NY
         WRITE (23,*) '       GRID SIZE           (DX, METER) :',LO%DX
	     WRITE (23,*) '       GRID SIZE           (DY, METER) :',LO%DY
	     WRITE (23,*) '       TIME STEP SIZE     (DT, SECOND) :',LO%DT
	  ENDIF
	  WRITE (23,*) '       NUMBER OF CHILD GRID LAYERS     :',LO%NUM_CHILD
	  DO I = 1,NUM_GRID
	  IF (LA(I)%LAYSWITCH.EQ.0 .AND. LA(I)%PARENT.EQ.LO%ID) THEN
	  WRITE (23,*) '               CHILD GRID LAYER ID     :',LA(I)%ID
	  ENDIF
	  ENDDO
	  WRITE (23,*) '       BATHYMETRY DATA FILE NAME       :',LO%DEPTH_NAME

!*	  CALL CR_CHECK (LO,LA)

	  DO I=1,NUM_GRID
	     IF (LA(I)%LAYSWITCH .EQ. 0) THEN
            WRITE (23,*) '--------------- SUB-LEVEL GRID INFORMATION -----------------'
		    WRITE (23,*) '    #  GRID IDENTIFICATION NUMBER      :',LA(I)%ID
	        IF (LA(I)%LAYCORD .EQ. 0) WRITE (23,*) '       USE SPHERICAL COORDINATE SYSTEM'
            IF (LA(I)%LAYCORD .EQ. 1) WRITE (23,*) '       USE CARTESIAN COORDINATE SYSTEM'
	        IF (LA(I)%LAYGOV .EQ. 0)  WRITE (23,*) '       USE LINEAR SHALLOW WATER EQUATIONS'
            IF (LA(I)%LAYGOV .EQ. 1)  WRITE (23,*) '       USE NONLINEAR SHALLOW WATER EQUATIONS'
	        IF (LA(I)%LAYGOV .EQ. 2)  WRITE (23,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
            IF (LA(I)%LAYGOV .EQ. 3)  WRITE (23,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME'
	        IF (LA(I)%LAYGOV .EQ. 4)  WRITE (23,*) '       USE LINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
            IF (LA(I)%LAYGOV .EQ. 5)  WRITE (23,*) '       USE NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME*'
	        IF (LA(I)%FRIC_SWITCH .EQ. 1) WRITE (23,*) '       BOTTOM FRICTION                 : DISABLED'
			IF (LA(I)%FRIC_SWITCH .NE. 1) WRITE (23,*) '       BOTTOM FRICTION                 : ENABLED'
			IF (LA(I)%FRIC_SWITCH .EQ. 0) WRITE (23,*) '          USE CONSTANT ROUGHNESS COEF. :',LO%FRIC_COEF
			IF (LA(I)%FRIC_SWITCH .EQ. 2) WRITE (23,*) '          USE VARIABLE MANNING ROUGHNESS COEFFICIENTS.'
			IF (LA(I)%FLUXSWITCH .EQ. 0) WRITE (23,*)  '       VOLUME FLUX OUTPUT              : ENABLED'
			IF (LA(I)%FLUXSWITCH .EQ. 1) WRITE (23,*)  '       VOLUME FLUX OUTPUT              : DISABLED'
			IF (LA(I)%FLUXSWITCH .NE. 2) WRITE (23,*)  '       SURFACE DISPLACEMENT OUTPUT     : ENABLED'
			IF (LA(I)%FLUXSWITCH .EQ. 2) WRITE (23,*)  '       SURFACE DISPLACEMENT OUTPUT     : DISABLED'
			WRITE (23,*) '       PARENT GRID ID                  :',LA(I)%PARENT
			WRITE (23,*) '       GRID SIZE RATIO                 :',LA(I)%REL_SIZE
			WRITE (23,*) '       TIME STEP SIZE RATIO            :',LA(I)%REL_TIME
			WRITE (23,*) '       POSITIONS IN ITS PARENT LAYER   :'
			IF (LA(I)%LAYCORD .EQ. 0) THEN
	           WRITE (23,*) '           X_START            (DEGREE) :',LA(I)%X(2)
	           WRITE (23,*) '           X_END              (DEGREE) :',LA(I)%X(LA(I)%NX)
	           WRITE (23,*) '           Y_START            (DEGREE) :',LA(I)%Y(2)
	           WRITE (23,*) '           Y_END              (DEGREE) :',LA(I)%Y(LA(I)%NY)
	           WRITE (23,*) '           I_START             (INDEX) :',LA(I)%CORNERS(1)
	           WRITE (23,*) '           I_END               (INDEX) :',LA(I)%CORNERS(2)
	           WRITE (23,*) '           J_START             (INDEX) :',LA(I)%CORNERS(3)
	           WRITE (23,*) '           J_END               (INDEX) :',LA(I)%CORNERS(4)
			   WRITE (23,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
			   WRITE (23,*) '       GRID SIZE           (DX,MINUTE) :',LA(I)%DX
			   WRITE (23,*) '       GRID SIZE           (DY,MINUTE) :',LA(I)%DY
			   WRITE (23,*) '       TIME STEP SIZE     (DT, SECOND) :',LA(I)%DT
			ENDIF
			IF (LA(I)%LAYCORD .EQ. 1) THEN
			   IF (LA(I)%UPZ) THEN
	           WRITE (23,*) '         X_START               (METER) :',LA(I)%X(2)
			   WRITE (23,*) '         X_END                 (METER) :',LA(I)%X(LA(I)%NX)
			   WRITE (23,*) '         Y_START               (METER) :',LA(I)%Y(2)
			   WRITE (23,*) '         Y_END                 (METER) :',LA(I)%Y(LA(I)%NY)
	           WRITE (23,*) '         I_START               (INDEX) :',LA(I)%CORNERS(1)
	           WRITE (23,*) '         I_END                 (INDEX) :',LA(I)%CORNERS(2)
	           WRITE (23,*) '         J_START               (INDEX) :',LA(I)%CORNERS(3)
	           WRITE (23,*) '         J_END                 (INDEX) :',LA(I)%CORNERS(4)
			   WRITE (23,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
		       ELSE
			   WRITE (23,*) '         X_START_S            (DEGREE) :',LA(I)%XT(2)
	           WRITE (23,*) '         X_END_S              (DEGREE) :',LA(I)%XT(LA(I)%NX)
			   WRITE (23,*) '         Y_START_S            (DEGREE) :',LA(I)%YT(2)
	           WRITE (23,*) '         Y_END_S              (DEGREE) :',LA(I)%YT(LA(I)%NY)
	           WRITE (23,*) '         X_START_C             (METER) :',LA(I)%X(2)
	           WRITE (23,*) '         X_END_C               (METER) :',LA(I)%X(LA(I)%NX)
	           WRITE (23,*) '         Y_START_C             (METER) :',LA(I)%Y(2)
	           WRITE (23,*) '         Y_END_C               (METER) :',LA(I)%Y(LA(I)%NY)
			   WRITE (23,*) '       GRID DIMENSION          (NX*NY) :',LA(I)%NX-1,'*',LA(I)%NY-1
		       ENDIF
			   WRITE (23,*) '       GRID SIZE           (DX, METER) :',LA(I)%DX
			   WRITE (23,*) '       GRID SIZE           (DY, METER) :',LA(I)%DY
			   WRITE (23,*) '       TIME STEP SIZE     (DT, SECOND) :',LA(I)%DT
			ENDIF
			WRITE (23,*) '       NUMBER OF CHILD GRID LAYERS     :',LA(I)%NUM_CHILD
			DO K = 1,NUM_GRID
			IF (LA(K)%LAYSWITCH.EQ.0 .AND. LA(K)%PARENT.EQ.LA(I)%ID) THEN
			WRITE (23,*) '               CHILD GRID LAYER ID     :',LA(K)%ID
			ENDIF
			ENDDO
			WRITE (23,*) '       BATHYMETRY DATA FILE NAME       :',LA(I)%DEPTH_NAME
			WRITE (23,*) '------------------------------------------------------------'
		 ENDIF
	  ENDDO
	  WRITE (23,*) ' '
	  CLOSE(23)

      RETURN
	  END



!----------------------------------------------------------------------
      SUBROUTINE WRITE_INI (LO)
! *********  OUTPUT INITIAL CONDITION DATA **************
!.....CREATED ON OCT.29 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO

      OPEN (25,FILE='ini_surface.dat',STATUS='UNKNOWN')
      DO J=1,LO%NY
         WRITE (25,'(15F8.3)') (LO%Z(I,J,1),I=1,LO%NX)
      ENDDO
      CLOSE (25)

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE BATHY_WRITE (LO)
!......................................................................
!DESCRIPTION:
!	  #. WRITE THE BATHYMETRY DATA ONTO HARDDRIVE;
!	  #. THREE DATA FILES WILL BE CREATED: LAYER##.DAT -- X COORDINATES
!        LAYER##_Y.DAT -- Y COORDINATES AND LAYER##.DAT -- WATER DEPTH;
!		 XX REPRESENTS THE 2-DIGIT IDENTIFICATION OF NUMBER OF A LAYER;
!NOTES:
!	  #. CREATED ON OCT 28 2008 (XIAOMING WANG, GNS)
!	  #. LAST REVISE: NOV.21 2008 (XIAOMING WANG)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  CHARACTER(LEN=20) FNAME1,FNAME2,FNAME3,FNAME4,FNAME5
	  CHARACTER(LEN=20) FNAME9
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      IF (LO%ID.EQ.1) THEN
	     IS = 1
		 JS = 1
	  ELSE
	     IS = 2
		 JS = 2
	  ENDIF
	  IE = LO%NX
	  JE = LO%NY

!.....WRITE BATHYMETRY DATA
      WRITE (FNAME1,1) LO%ID
 1    FORMAT('layer',I2.2,'.dat')
      OPEN (25,FILE=FNAME1,STATUS='UNKNOWN')
      DO J = JS,JE
	     DO I = IS,IE
            WRITE (25,'(F12.4)') LO%H(I,J)
		 ENDDO
      ENDDO
      CLOSE (25)

!.....WRITE MODIFIED BATHYMETRY DATA
      IF ( LO%FLUXSWITCH .GE. 9 ) THEN
         WRITE (FNAME9,9) LO%ID
 9       FORMAT('layer',I2.2,'m.dat')
         OPEN (25,FILE=FNAME9,STATUS='UNKNOWN')
         DO J = JS,JE
	        DO I = IS,IE
               WRITE (25,'(F12.4)') LO%H(I,J)
		    ENDDO
         ENDDO
         CLOSE (25)
	  ENDIF


!.....WRITE X COORDINATES
      WRITE (FNAME2,2) LO%ID
 2    format('layer',i2.2,'_x.dat')
      OPEN (25,FILE=FNAME2,STATUS='UNKNOWN')
      DO I = IS,IE
         WRITE (25,'(F17.6)') LO%X(I)
      ENDDO
      CLOSE (25)
!.....WRITE Y COORDINATES
      WRITE (FNAME3,3) LO%ID
 3    format('layer',i2.2,'_y.dat')
      OPEN (25,FILE=FNAME3,STATUS='UNKNOWN')
      DO J = JS,JE
         WRITE (25,'(F17.6)') LO%Y(J)
      ENDDO
      CLOSE (25)

	  IF (.NOT.LO%UPZ) THEN
!.....WRITE LO%XT
      WRITE (FNAME4,4) LO%ID
 4    format('layer',i2.2,'_xs.dat')
      OPEN (25,FILE=FNAME4,STATUS='UNKNOWN')
      DO I = IS,IE
         WRITE (25,'(F17.6)') LO%XT(I)
      ENDDO
      CLOSE (25)
!.....WRITE LO%YT
      WRITE (FNAME5,5) LO%ID
 5    format('layer',i2.2,'_ys.dat')
      OPEN (25,FILE=FNAME5,STATUS='UNKNOWN')
      DO J = JS,JE
         WRITE (25,'(F17.6)') LO%YT(J)
      ENDDO
      CLOSE (25)
	  ENDIF

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE ALL_PRT (K,LO,LA)
! *********  OUTPUT RESULTS IN ASCII FORMAT  **************
! PREPARED BY TOM LOGAN, ARSC (2005)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
!      IMPLICIT NONE
      TYPE (LAYER) :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
      INTEGER      :: K
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      CALL DATA_PRT (K,LO)
      DO I=1,NUM_GRID
	     IF (LA(I)%LAYSWITCH .EQ. 0)  CALL DATA_PRT (K,LA(I))
      ENDDO

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE DATA_PRT (K,LO)
!......................................................................
!DESCRIPTION:
!	  #. OUTPUT FREE SURFACE ELEVATION AND VOLUME FLUXES;
!NOTES:
!	  #. CREATED ON ??? ??, ???? (XIAOMING WANG, CORNELL UNIVERSITY)
!	  #. UPDATED ON MAR 01 2009 (XIAOMING WANG, GNS)
!		 1. ADD TIDAL LEVEL CORRECTION;
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO
	  REAL Z(LO%NX,LO%NY)
	  INTEGER           :: K
	  INTEGER  M,N
      CHARACTER(LEN=20) FNAME
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
                    
      if ( k .GT. 0 ) then
          CALL CUDA_GETZ(LO%Z(:,:,1), LO%ID)
          CALL CUDA_GETMN(LO%M(:,:,1), LO%N(:,:,1), LO%ID)
      end if
      
	  Z = 0.0
	  Z(:,:) = LO%Z(:,:,1)
	  DO I = 1,LO%NX
		 DO J = 1,LO%NY
		    IF (Z(I,J)+LO%H(I,J).LE.GX .AND. LO%H(I,J).LT.ZERO) THEN
			   Z(I,J) = 0.0
			ELSE
			   !APPLY TIDAL LEVEL CORRECTION
	           IF (ABS(LO%TIDE_LEVEL).GT.GX) THEN
			      Z(I,J) = LO%Z(I,J,1) + LO%TIDE_LEVEL
			   ENDIF
			ENDIF
		 ENDDO
	  ENDDO



!     MODIFIED TO REMOVE ADDITIONAL COLUMN AND ROW
      IF (LO%PARENT .LE. 0) THEN
	      IS=1
		  JS=1
		  IE=LO%NX
		  JE=LO%NY
	  ELSE
	      IS=2
		  JS=2
		  IE=LO%NX
		  JE=LO%NY
	  ENDIF
!.....OUTPUT FREE SURFACE ELEVATION
      IF (LO%FLUXSWITCH .LE. 1 .OR. LO%FLUXSWITCH .EQ. 9) THEN
         WRITE (FNAME,1) LO%ID,K
 1       FORMAT('z_',I2.2,'_',I6.6,'.dat')
         OPEN (25,FILE=FNAME,STATUS='UNKNOWN',form="unformatted")
         !DO J = JS,JE
            WRITE (25) (Z(:,:))!(I,J),I=IS,IE)
         !ENDDO
         CLOSE (25)
	  ENDIF
!.....OUTPUT VOLUME FLUX
	  IF (LO%FLUXSWITCH .EQ. 0) THEN
         WRITE (FNAME,2) LO%ID,K
 2       FORMAT('m_',I2.2,'_',I6.6,'.dat')

         OPEN (25,FILE=FNAME,STATUS='UNKNOWN',form="unformatted")
         !DO J = JS,JE
           WRITE (25) (LO%M(:,:,1))!(I,J,1),I=IS,IE)
         !ENDDO
         CLOSE (25)
         WRITE (FNAME,3) LO%ID,K
 3       FORMAT('n_',I2.2,'_',I6.6,'.dat')

         OPEN (25,FILE=FNAME,STATUS='UNKNOWN',form="unformatted")
         !DO J = JS,JE
           WRITE (25) (LO%N(:,:,1))!(I,J,1),I=IS,IE)
         !ENDDO
         CLOSE (25)
	  ENDIF

	  IF (LO%SEDI_SWITCH .EQ. 0) THEN
         WRITE (FNAME,4) LO%ID,K
 4       FORMAT('s_',I2.2,'_',I6.6,'.dat')
         OPEN (25,FILE=FNAME,STATUS='UNKNOWN',form="unformatted")
         !DO J = JS,JE
            WRITE (25) (LO%DH(:,:,2))!(I,J,2),I=IS,IE)
         !ENDDO
         CLOSE (25)
      ENDIF

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE MAX_AMP (LO,LA,TIME,TEND)
!DESCRIPTION:
!	  #. OUTPUT MAX SURFACE ELEVATION/DEPRESSION EVERY HALF HOUR AND
!		 AT THE END OF A SIMULATION;
!NOTE:
!	  #. CREATED ON NOV.21, 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON ???
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  REAL SECS
	  REAL TIME,TEND,TEMP
	  REAL, ALLOCATABLE :: TMP(:,:)
	  INTEGER K,K1,KK, I, J
	  CHARACTER(LEN=40) FNAME1,FNAME2,FNAME3
	  CHARACTER(LEN=40) FNAME4,FNAME5,FNAME6
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
     !!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
      !OBTAIN MAXIMUM ELEVATION AND DEPRESSION FOR THE OUTEST LAYER
      ! DO I = 1,LO%NX
	  !    DO J = 1,LO%NY
		!     IF (LO%H(I,J)+LO%Z(I,J,2) .GT. GX) THEN
		! 		IF (LO%Z(I,J,2).GT.LO%Z_MAX(I,J))					&
		! 							LO%Z_MAX(I,J) = LO%Z(I,J,2)
		! 		! IF (LO%Z(I,J,2).LT.LO%Z_MIN(I,J))					&
		! 		! 					LO%Z_MIN(I,J) = LO%Z(I,J,2)
		! 	ENDIF
		!  ENDDO
      ! ENDDO
	  !OBTAIN MAX ELEVATION AND DEPRESSION FOR THE INNER LAYERS, LA



	  ! DO K = 1,NUM_GRID
		!  IF (LA(K)%LAYSWITCH.EQ.0) THEN
      !       DO I = 1,LA(K)%NX
	  !          DO J = 1,LA(K)%NY
		!           IF (LA(K)%H(I,J)+LA(K)%Z(I,J,2) .GT. GX) THEN
		! 		     IF (LA(K)%Z(I,J,2).GT.LA(K)%Z_MAX(I,J))		&
		! 						LA(K)%Z_MAX(I,J) = LA(K)%Z(I,J,2)
		! 		     IF (LA(K)%Z(I,J,2).LT.LA(K)%Z_MIN(I,J))		&
		! 						LA(K)%Z_MIN(I,J) = LA(K)%Z(I,J,2)
		! 	      ENDIF
		!        ENDDO
      !       ENDDO
		!  ENDIF
	  ! ENDDO
      !!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!

!.....WRITE THE MAXIMUM ELEVATION AND DEPRESSION INTO DATA FILE
!	  EVERY ONE HOUR AND AT THE END OF THE SIMULATION
!.....OUTPUT FREE SURFACE ELEVATION
!     FILE NAME CONVENTION:
!     MAXIMUM ELEVATION:   ETAMAX_LAYER##_####MINS.DAT
!     MAXIMUM DEPRESSION:  ETAMIN_LAYER##_####MINS.DAT
      K1 = FLOOR((TIME-LO%DT)/1800.0)
	  IF (K1 .LE. 0) K1 = 0
	  K2 = FLOOR(TIME/1800.0)
      K = NINT(TIME/60.0)
	  SECS = 60.0*K
      IF (K2.GT.K1 .OR. TIME.GE.TEND-2.0*LO%DT) THEN
         CALL CUDA_GETZMAX(LO%Z_MAX(:,:), LO%ID)
		 ALLOCATE(TMP(LO%NX,LO%NY))
		 DO I = 1,LO%NX
		    DO J = 1,LO%NY
			   IF ((LO%H(I,J)+LO%Z_MAX(I,J)).LE.GX .AND. 			&
											LO%H(I,J).LT.ZERO) THEN
			      TMP(I,J) = 0.0
			   ELSE
			      TMP(I,J) = LO%Z_MAX(I,J)
				  IF (ABS(LO%TIDE_LEVEL).GT.GX) THEN
				     TMP(I,J) = LO%Z_MAX(I,J) + LO%TIDE_LEVEL
				  ENDIF
			   ENDIF
		    ENDDO
		 ENDDO
         WRITE (FNAME1,1) LO%ID,K
 1       FORMAT('zmax_layer',I2.2,'_',I4.4,'mins.dat')
         IF (TIME.GE.TEND-2.0*LO%DT) FNAME1 = 'zmax_layer01.dat'
         OPEN (25,FILE=FNAME1,STATUS='UNKNOWN',form="unformatted")
         !DO J = 1,LO%NY
		    WRITE (25) (TMP(:,:))!(I,J),I=1,LO%NX)
         !ENDDO
         CLOSE (25)
 ! !!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
	! 	 DO I = 1,LO%NX
	! 	    DO J = 1,LO%NY
	! 		   IF ((LO%H(I,J)+LO%Z_MIN(I,J)).LE.GX .AND. 			&
	! 										LO%H(I,J).LT.ZERO) THEN
	! 		      TMP(I,J) = 0.0
	! 		   ELSE
	! 		      TMP(I,J) = LO%Z_MIN(I,J)
	! 			  IF (ABS(LO%TIDE_LEVEL).GT.GX) THEN
	! 			     TMP(I,J) = LO%Z_MIN(I,J) + LO%TIDE_LEVEL
	! 			  ENDIF
	! 		   ENDIF
	! 	    ENDDO
	! 	 ENDDO
 !         WRITE (FNAME2,2) LO%ID,K
 ! 2       FORMAT('zmin_layer',I2.2,'_',I4.4,'hrs.dat')
 !         IF (TIME.GE.TEND-2.0*LO%DT) FNAME2 = 'zmin_layer01.dat'
 !         OPEN (25,FILE=FNAME2,STATUS='UNKNOWN')
 !         DO J = 1,LO%NY
 !            WRITE (25,'(15F9.4)') (TMP(I,J),I=1,LO%NX)
 !         ENDDO
 !         CLOSE (25)
!!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!

		 DEALLOCATE(TMP,STAT = ISTAT)
!!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
 !         DO KK = 1,NUM_GRID
 !            IF (LA(KK)%LAYSWITCH .EQ. 0 ) THEN
	! 		   ALLOCATE(TMP(LA(KK)%NX,LA(KK)%NY))
	! 		   TMP = 0.0
	! 	 	   DO I = 2,LA(KK)%NX
	! 	          DO J = 2,LA(KK)%NY
	! 		         IF ((LA(KK)%H(I,J)+LA(KK)%Z_MAX(I,J)).LE.GX 	&
	! 							.AND. LA(KK)%H(I,J).LT.ZERO) THEN
	! 		            TMP(I,J) = 0.0
	! 		         ELSE
	! 		            TMP(I,J) = LA(KK)%Z_MAX(I,J)
	! 			        IF (ABS(LA(KK)%TIDE_LEVEL).GT.GX) THEN
	! 			           TMP(I,J) = LA(KK)%Z_MAX(I,J)+LA(KK)%TIDE_LEVEL
	! 			        ENDIF
	! 		         ENDIF
	! 	          ENDDO
	! 	       ENDDO
	! 		   IF (TIME.GE.TEND-2.0*LO%DT) THEN
 !                  WRITE (FNAME3,5) LA(KK)%ID
 ! 5                FORMAT('zmax_layer',I2.2,'.dat')
	! 		   ELSE
 !                  WRITE (FNAME3,3) LA(KK)%ID,K
 ! 3                FORMAT('zmax_layer',I2.2,'_',I4.4,'hrs.dat')
 !               ENDIF
 !               OPEN (25,FILE=FNAME3,STATUS='UNKNOWN')
 !               DO J = 2,LA(KK)%NY
 !             	  WRITE (25,'(15F9.4)') (TMP(I,J),I=2,LA(KK)%NX)
 !               ENDDO
 !               CLOSE (25)
 !
	! 		   TMP = 0.0
	! 	 	   DO I = 2,LA(KK)%NX
	! 	          DO J = 2,LA(KK)%NY
	! 		         IF ((LA(KK)%H(I,J)+LA(KK)%Z_MIN(I,J)).LE.GX 	&
	! 							.AND. LA(KK)%H(I,J).LT.ZERO) THEN
	! 		            TMP(I,J) = 0.0
	! 		         ELSE
	! 		            TMP(I,J) = LA(KK)%Z_MIN(I,J)
	! 			        IF (ABS(LA(KK)%TIDE_LEVEL).GT.GX) THEN
	! 			           TMP(I,J) = LA(KK)%Z_MIN(I,J)+LA(KK)%TIDE_LEVEL
	! 			        ENDIF
	! 		         ENDIF
	! 	          ENDDO
	! 	       ENDDO
	! 		   IF (TIME.GE.TEND-2.0*LO%DT) THEN
 !                  WRITE (FNAME4,6) LA(KK)%ID
 ! 6                FORMAT('zmin_layer',I2.2,'.dat')
	! 		   ELSE
 !                  WRITE (FNAME4,4) LA(KK)%ID,K
 ! 4                FORMAT('zmin_layer',I2.2,'_',I4.4,'hrs.dat')
 !               ENDIF
 !               OPEN (25,FILE=FNAME4,STATUS='UNKNOWN')
 !               DO J = 2,LA(KK)%NY
 !               	  WRITE (25,'(15F9.4)') (TMP(I,J),I=2,LA(KK)%NX)
 !               ENDDO
 !               CLOSE (25)
	! 		   DEALLOCATE(TMP,STAT = ISTAT)
 !            ENDIF
	! 	 ENDDO
    !!!!!!!!!!!!!!!!!!!!      COMMEMT ADDED BY TAO     !!!!!!!!!!!!!!!!!!!!!!!!!!
	  ENDIF


	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE GET_TS (DAT,X0,Y0,LO,LA,ID)
!DESCRIPTION:
!	  #. EXTRACT TIME HISTORY RECORD AT SPECIFIED LOCATION;
!NOTE:
!	  #. CREATED ON NOV 25 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN 30 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON FEB08 2009 (XIAOMING WANG, GNS)
!		 1. WATER SURFACE ELEVATION ON DRY LAND WILL TAKE THE VALUE OF
!			LAND ELEVATION, NOW CHANGE IT TO ZERO.
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  REAL DAT,X0,Y0
      REAL L1,L2,L3,L4
	  INTEGER ID,IS,JS,IE,JE,KI,KJ
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!*	  ID = 1
	  DAT = 0.0

!*	  CALL GAUGE_LAYERID (ID,X0,Y0,LO,LA)

	  IF (ID.EQ.1) THEN
         IS = 1
	     JS = 1
	     IE = LO%NX
	     JE = LO%NY

         KI = 0
	     KJ = 0
         DO I = IS,IE-1
		    IF (X0.GE.LO%X(I) .AND. X0.LT.LO%X(I+1)) THEN
		       KI = I
		    END IF
		    IF (X0.LT.LO%X(IS)) KI = IS
		    IF (X0.GE.LO%X(IE)) KI = IE-1
	     END DO
         DO J = JS,JE-1
		    IF (Y0.GE.LO%Y(J) .AND. Y0.LT.LO%Y(J+1)) THEN
		       KJ = J
		    END IF
		    IF (Y0.LT.LO%Y(JS)) KJ = JS
		    IF (Y0.GE.LO%Y(JE)) KJ = JE-1
	     END DO
	     DELTA_X = LO%X(KI+1)-LO%X(KI)
	     DELTA_Y = LO%Y(KJ+1)-LO%Y(KJ)
	     CX = (X0-LO%X(KI))/DELTA_X
	     CY = (Y0-LO%Y(KJ))/DELTA_Y
         CALL GCOMCOT_GET_Z(1, L1, L2, L3, L4, KI, KJ)
         Z1 = L1*(1.0-CX)*(1.0-CY)
	     Z2 = L2*(CX)*(1-CY)
	     Z3 = L3*(1.0-CX)*(CY)
	     Z4 = L4*(CX)*(CY)
         
!        Z1 = LO%Z(KI,KJ,2)*(1.0-CX)*(1.0-CY)
!	     Z2 = LO%Z(KI+1,KJ,2)*(CX)*(1-CY)
!	     Z3 = LO%Z(KI,KJ+1,2)*(1.0-CX)*(CY)
!	     Z4 = LO%Z(KI+1,KJ+1,2)*(CX)*(CY)
	     DAT = Z1+Z2+Z3+Z4
		 IF (ABS(LO%TIDE_LEVEL).GT.GX) DAT = DAT + LO%TIDE_LEVEL
	  ELSE
		 K = ID-1
         IS = 2
	     JS = 2
	     IE = LA(K)%NX
	     JE = LA(K)%NY

         KI = 0
	     KJ = 0
         DO I = IS,IE-1
		    IF (X0.GE.LA(K)%X(I) .AND. X0.LT.LA(K)%X(I+1)) THEN
		       KI = I
		    END IF
		    IF (X0.LT.LA(K)%X(IS)) KI = IS
		    IF (X0.GE.LA(K)%X(IE)) KI = IE-1
	     END DO
         DO J = JS,JE-1
		    IF (Y0.GE.LA(K)%Y(J) .AND. Y0.LT.LA(K)%Y(J+1)) THEN
		    KJ = J
		    END IF
		    IF (Y0.LT.LA(K)%Y(JS)) KJ = JS
		    IF (Y0.GE.LA(K)%Y(JE)) KJ = JE-1
	     END DO
	     DELTA_X = LA(K)%X(KI+1)-LA(K)%X(KI)
	     DELTA_Y = LA(K)%Y(KJ+1)-LA(K)%Y(KJ)
	     CX = (X0-LA(K)%X(KI))/DELTA_X
	     CY = (Y0-LA(K)%Y(KJ))/DELTA_Y

!		 ETALL = LA(K)%Z(KI,KJ,2)
!		 ETALR = LA(K)%Z(KI+1,KJ,2)
!		 ETAUL = LA(K)%Z(KI,KJ+1,2)
!		 ETAUR = LA(K)%Z(KI+1,KJ+1,2)
!		 IF (ETALL+LA(K)%H(KI,KJ) .LT. GX) ETALL = 0.0
!		 IF (ETALR+LA(K)%H(KI+1,KJ) .LT. GX) ETALR = 0.0
!		 IF (ETAUL+LA(K)%H(KI,KJ+1) .LT. GX) ETAUL = 0.0
!		 IF (ETAUR+LA(K)%H(KI+1,KJ+1) .LT. GX) ETAUR = 0.0
!
!        Z1 = ETALL*(1.0-CX)*(1.0-CY)
!	     Z2 = ETALR*CX*(1.0-CY)
!	     Z3 = ETAUL*(1.0-CX)*CY
!	     Z4 = ETAUR*CX*CY
!	     DAT = Z1+Z2+Z3+Z4
!
         CALL GCOMCOT_GET_Z(LA(K)%ID, L1, L2, L3, L4, KI, KJ)
         Z1 = L1*(1.0-CX)*(1.0-CY)
	     Z2 = L2*(CX)*(1-CY)
	     Z3 = L3*(1.0-CX)*(CY)
	     Z4 = L4*(CX)*(CY)

!        Z1 = LA(K)%Z(KI,KJ,2)*(1.0-CX)*(1.0-CY)
!	     Z2 = LA(K)%Z(KI+1,KJ,2)*(CX)*(1-CY)
!	     Z3 = LA(K)%Z(KI,KJ+1,2)*(1.0-CX)*(CY)
!	     Z4 = LA(K)%Z(KI+1,KJ+1,2)*(CX)*(CY)
	     DAT = Z1+Z2+Z3+Z4
		 IF (ABS(LA(K)%TIDE_LEVEL).GT.GX) DAT = DAT + LA(K)%TIDE_LEVEL
	  ENDIF

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE GAUGE_LAYERID (ID,X0,Y0,LO,LA)
!DESCRIPTION:
!	  #. DETERMINE WHICH GRID LAYER THE TIDAL GAGE IS LOCATED IN;
!INPUT:
!	  #. TIDAL GAGE LOCATION: X0, Y0;
!	  #. GRID LAYER INFO: LO, LA;
!OUTPUT:
!	  #. GRID LAYER ID IN WHICH THE TIDAL GAGE IS LOCATED;
!NOTE:
!	  #. CREATED ON JAN 30 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN 30 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  INTEGER ID
	  REAL DAT,X0,Y0
	  INTEGER IS,JS,IE,JE,KI,KJ
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  ID = 1
	  IDTMP = 1
	  LEVEL = 1
	  LTMP = 1

	  IS = 2
	  JS = 2
	  DO K = 1, NUM_GRID
		 IE = LA(K)%NX
		 JE = LA(K)%NY
		 IF (LA(K)%LAYSWITCH.EQ.0) THEN
			IF (X0.GE.LA(K)%X(IS) .AND. X0.LT.LA(K)%X(IE) .AND. 	&
				Y0.GE.LA(K)%Y(JS) .AND. Y0.LT.LA(K)%Y(JE)) THEN
			   IDTMP = LA(K)%ID
			   LTMP = LA(K)%LEVEL
			ENDIF
			IF (LTMP.GT.LEVEL) THEN
			   LEVEL = LTMP
			   ID = IDTMP
		    ENDIF
		 ENDIF
	  ENDDO

!	  WRITE (*,*) ID

	  RETURN
	  END
