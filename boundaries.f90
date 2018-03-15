!*!----------------------------------------------------------------------
!*      SUBROUTINE BC_SETUP (LO,BCI_INFO,WAVE_INFO,BC_TYPE)
!*!.....SET UP VERTICAL WALL BOUNDARY
!*!----------------------------------------------------------------------
!*      USE LAYER_PARAMS
!*	  USE WAVE_PARAMS
!*	  USE BCI_PARAMS
!*	  TYPE (BCI)   :: BCI_INFO
!* 	  TYPE (LAYER) :: LO
!*	  TYPE (WAVE)  :: WAVE_INFO
!*	  INTEGER BC_TYPE
!*
!*      IF (BC_TYPE.EQ.0) CALL OPEN (LO)
!*	  IF (BC_TYPE.EQ.1) CALL SPONGE_LAYER (LO)
!*	  IF (BC_TYPE.EQ.2) CALL BC_WALL (LO,WAVE_INFO)
!*!	  IF (BC_TYPE.EQ.3) CALL BC_INPUT (LO,BCI_INFO,TIME)
!*
!*	  RETURN
!*	  END

!----------------------------------------------------------------------
      SUBROUTINE OPEN (LO)
!.....DEPLOY OPEN BOUNDARY CONDITION (ONLY FOR OUTEST LAYER)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO
      LOGICAL IFIRST,JFIRST,IEND,JEND
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/
!
      I_S = 1
	  J_S = 1
	  I_E = LO%NX
	  J_E = LO%NY
      IFIRST = I_S .EQ. 1
      JFIRST = J_S .EQ. 1
      IEND = I_E .EQ. LO%NX
      JEND = J_E .EQ. LO%NY
!
      IF ( JFIRST ) THEN
         J = 1
         DO I=2,LO%NX-1
            IF (LO%H(I,J) .GT. GX) THEN
               CC = SQRT(GRAV*LO%H(I,J))
               UH = 0.5*(LO%M(I,J,2)+LO%M(I-1,J,2))
               UU = SQRT(UH**2+LO%N(I,J,2)**2)
               ZZ = UU/CC
               ARG = LO%N(I,J,2)
               IF (ARG .GT. ZERO) THEN
                  ZZ = -ZZ
               ENDIF
!			   IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
			   IF (ABS(ZZ) .GT. UB) ZZ=0.0
               LO%Z(I,J,2) = ZZ
            ELSE
               LO%Z(I,J,2) = ZERO
            ENDIF
            ! if ( i .EQ. 715 ) then
            !     WRITE(*,*) "Fortran--------------------->dbdb"
            !     WRITE(*,*) LO%H(I,1), LO%M(I,1,2), LO%N(I,1,2), 1/CC
            !     WRITE(*,*) LO%M(I-1,1,2), ZZ
            ! end if
         END DO
      ENDIF
!
      IF ( JEND ) THEN
         J = LO%NY
         DO I=2,LO%NX-1
            IF (LO%H(I,J) .GT. GX) THEN
               CC = SQRT(GRAV*LO%H(I,J))
               UH = 0.5*(LO%M(I,J,2)+LO%M(I-1,J,2))
               UU = SQRT(UH**2+LO%N(I,J-1,2)**2)
               ZZ = UU/CC
               ARG = LO%N(I,J-1,2)
               IF (ARG .LT. ZERO) THEN
                  ZZ = -ZZ
               ENDIF
!			   IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
			   IF (ABS(ZZ) .GT. UB) ZZ=0.0
               LO%Z(I,J,2) = ZZ
            ELSE
               LO%Z(I,J,2) = ZERO
            ENDIF
         END DO
      ENDIF
!
      IF ( IFIRST ) THEN
         I = 1
         DO J=2,LO%NY-1
            IF (LO%H(I,J) .GT. GX) THEN
               CC = SQRT(GRAV*LO%H(I,J))
               IF (LO%H(I,J-1) .GT. GX) THEN
                  UH = 0.5*(LO%N(I,J,2)+LO%N(I,J-1,2))
               ELSE
                  UH = LO%N(I,J,2)
	           ENDIF
               UU = SQRT(UH**2+LO%M(I,J,2)**2)
               ZZ = UU/CC
               ARG = LO%M(I,J,2)
               IF (ARG .GT. ZERO) THEN
                  ZZ = -ZZ
               ENDIF
!			   IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
			   IF (ABS(ZZ) .GT. UB) ZZ=0.0
               LO%Z(I,J,2) = ZZ
            ELSE
               LO%Z(I,J,2) = ZERO
            ENDIF
         END DO
      ENDIF
!
      IF ( IEND ) THEN
         I = LO%NX
         DO J=2,LO%NY-1
            IF (LO%H(I,J) .GT. GX) THEN
               CC = SQRT(GRAV*LO%H(I,J))
               UH = 0.5*(LO%N(I,J,2)+LO%N(I,J-1,2))
               UU = SQRT(UH**2+LO%M(I-1,J,2)**2)
               ZZ = UU/CC
               ARG = LO%M(I-1,J,2)
               IF (ARG .LT. ZERO) THEN
                  ZZ = -ZZ
               ENDIF
!			   IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
			   IF (ABS(ZZ) .GT. UB) ZZ=0.0
               LO%Z(I,J,2) = ZZ
            ELSE
               LO%Z(I,J,2) = ZERO
            ENDIF
         END DO
      ENDIF
!
      IF ( IFIRST .AND. JFIRST ) THEN
         IF (LO%H(1,1) .GT. GX) THEN
            QX = LO%M(1,1,2)
            QY = LO%N(1,1,2)
            CC = SQRT(GRAV*LO%H(1,1))
            UH = SQRT(QX**2+QY**2)
            ZZ = UH/CC
            IF (QX .GT. ZERO .OR. QY.GT.ZERO) ZZ = -ZZ
!			IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
		    IF (ABS(ZZ) .GT. UB) ZZ=0.0
            LO%Z(1,1,2) = ZZ
         ELSE
            LO%Z(1,1,2) = ZERO
         ENDIF
      ENDIF
!
      IF ( IEND .AND. JFIRST ) THEN
         IF (LO%H(LO%NX,1) .GT. GX) THEN
            QX = LO%M(LO%NX-1,1,2)
            QY = LO%N(LO%NX,1,2)
            CC = SQRT(GRAV*LO%H(LO%NX,1))
            UH = SQRT(QX**2+QY**2)
            ZZ = UH/CC
            IF (QX .LT. ZERO .OR. QY.GT.ZERO) ZZ = -ZZ
!			IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
		    IF (ABS(ZZ) .GT. UB) ZZ=0.0
            LO%Z(LO%NX,1,2) = ZZ
         ELSE
            LO%Z(LO%NX,1,2) = ZERO
         ENDIF
      ENDIF
!
      IF ( IFIRST .AND. JEND ) THEN
         IF (LO%H(1,LO%NY) .GT. GX) THEN
            QX = LO%M(1,LO%NY,2)
            QY = LO%N(1,LO%NY-1,2)
            CC = SQRT(GRAV*LO%H(1,LO%NY))
            UH = SQRT(QX**2+QY**2)
            ZZ = UH/CC
            IF (QX .GT. ZERO .OR. QY.LT.ZERO) ZZ = -ZZ
!			IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
		    IF (ABS(ZZ) .GT. UB) ZZ = 0.0
            LO%Z(1,LO%NY,2) = ZZ
         ELSE
            LO%Z(1,LO%NY,2) = ZERO
         ENDIF
      ENDIF
!
      IF ( IEND .AND. JEND ) THEN
         IF (LO%H(LO%NX,LO%NY) .GT. GX) THEN
            QX = LO%M(LO%NX-1,LO%NY,2)
            QY = LO%N(LO%NX,LO%NY-1,2)
            CC = SQRT(GRAV*LO%H(LO%NX,LO%NY))
            UH =  SQRT(QX**2+QY**2)
			ZZ = UH/CC
			IF (QX.LT.ZERO .OR. QY.LT.ZERO) ZZ = -ZZ
!			IF (ABS(ZZ) .LE. EPS) ZZ = 0.0
			LO%Z(LO%NX,LO%NY,2) = ZZ
!	        LO%Z(LO%NX,LO%NY,2)=0.5*(LO%Z(LO%NX-1,LO%NY,2)			&
!									+ LO%Z(LO%NX,LO%NY-1,2))
		    IF (ABS(LO%Z(LO%NX,LO%NY,2)) .GT. UB)					&
					LO%Z(LO%NX,LO%NY,2)=0.0
         ELSE
            LO%Z(LO%NX,LO%NY,2) = ZERO
         ENDIF
      ENDIF
!
      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE SPONGE_LAYER (LO)
!DESCRIPTION:
!     #. DETERMINE COEFFICIENT ALPHA SO AS TO IMPROVE NUMERICAL DISPERSION
!        REF: WEI AND KIRBY (1995) AND KIRBY ET AL (1998);
!	  #, USED TOGETHER WITH RADIATION/OPEN BOUNDARY CONDITION;
!NOTES:
!	  #. CREATED ON ??? ?? 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!	  #. UPDATED ON FEB 25 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      DO I = 1,LO%NX
	     IP1 = I+1
		 IM1 = I-1
		 IF (IP1.GE.LO%NX) IP1 = LO%NX
		 IF (IM1.LE.1) IM1 = 1
	     DO J = 1,LO%NY
			JP1 = J+1
			JM1 = J-1
			IF (JP1.GE.LO%NY) JP1 = LO%NY
			IF (JM1.LE.1) JM1 = 1
		    IF ( LO%SPONGE_COEFX(I,J).GT.ZERO .OR. 					&
							LO%SPONGE_COEFY(I,J).GT.ZERO ) THEN
   				CX = LO%SPONGE_COEFX(I,J)
				CY = LO%SPONGE_COEFY(I,J)
!				CXY = SQRT(CX**2+CY**2)
!*			    LO%M(I,J,2) = LO%M(I,J,2) + CX*LO%DT				&
!*									*0.5*(LO%M(I,J,1)+LO%M(I,J,2))
!*			    LO%N(I,J,2) = LO%N(I,J,2) + CY*LO%DT				&
!*									*0.5*(LO%N(I,J,1)+LO%N(I,J,2))
			    LO%M(I,J,2) = LO%M(I,J,2) + CX*LO%DT*LO%M(I,J,1)
			    LO%N(I,J,2) = LO%N(I,J,2) + CY*LO%DT*LO%N(I,J,1)
			ENDIF
		 ENDDO
      ENDDO

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE SPONGE_COEF (LO)
!DESCRIPTION:
!     DETERMINE COEFFICIENT ALPHA SO AS TO IMPROVE NUMERICAL DISPERSION
!.....REF: WEI AND KIRBY (1995) AND KIRBY ET AL (1998)
!NOTES:
!	  #. CREATED ON ??? ?? 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!	  #. UPDATED ON FEB 25 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
!	  USE WAVE_PARAMS
!	  USE FAULT_PARAMS
	  TYPE (LAYER) :: LO
!	  TYPE (WAVE)  :: WAVE_INFO
!	  TYPE (FAULT) :: FAULT_INFO
	  REAL DX, DY, DT, G, R_MAX, RX, RY, X_REL, R, WIDTH,H_MEAN
	  REAL COEFX(LO%NX,LO%NY),COEFY(LO%NX,LO%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      DX = LO%DX
	  DY = LO%DY
	  DT = LO%DT
      !WIDTH OF SPONGE LAYER IS 1.5 TIMES THE CHARACTERISTIC WAVE LENGTH

      !TYPICAL WATER LENGTH, USED TO DETERMINE CHARACTERSTIC WAVE NUMBER
	  CALL H_CALC(LO,H_MEAN,H_MAX)

	  DEPTH = H_MEAN

	  IF (DEPTH.LE.GX) DEPTH = GX

!DETERMINE CHARACTERISTIC WAVE LENGTH
	  WAVELENGTH = 20.0*DEPTH
	  WIDTH = 2.0*WAVELENGTH

	  !PARAMETERS USED BY KIRBY ET AL (1998)
	  ALPHA_C = 10.0;
	  ALPHA_MU = ZERO;
	  M = 2.0;
	  T = WAVELENGTH/SQRT(GRAV*DEPTH)
	  OMEGA = 2.0*PI/T

!.....IF LAYER01 ADOPTS SPHERICAL COORD,
!		CONVERT WIDTH (IN METERS) TO ARC MINUTES
	  IF (LO%LAYCORD .EQ. 0) THEN
	     SPONGE_WIDTH = WIDTH/R_EARTH*180.0/PI
	  ELSE
	     SPONGE_WIDTH = WIDTH
	  ENDIF

	  R_MAX = WIDTH !MAXIMUM DISTANCE FOR DAMPING

	  XS = LO%X(1) + SPONGE_WIDTH
	  YS = LO%Y(1) + SPONGE_WIDTH
	  XE = LO%X(LO%NX) - SPONGE_WIDTH
	  YE = LO%Y(LO%NY) - SPONGE_WIDTH

      COEFX = ZERO
	  COEFY = ZERO
      DO I = 1,LO%NX
	     DO J = 1,LO%NY
			X0 = LO%X(I)
			Y0 = LO%Y(J)
            IF (X0 .LE. XS) THEN
			   RX = -(X0 - XS)
			ELSEIF (X0 .GE. XE) THEN
			   RX = X0 - XE
			ELSE
			   RX = ZERO
			ENDIF
            IF (Y0 .LE. YS) THEN
			   RY = -(Y0 - YS)
			ELSEIF (Y0 .GE. YE) THEN
			   RY = Y0 - YE
			ELSE
			   RY = ZERO
			ENDIF
            R = SQRT(RX**2+RY**2)
			X_REL = RX
			Y_REL = RY
			IF (X_REL .LE. ZERO) X_REL = ZERO
			IF (Y_REL .LE. ZERO) Y_REL = ZERO
			IF (X_REL .GE. R_MAX) X_REL = R_MAX
			IF (Y_REL .GE. R_MAX) Y_REL = R_MAX

			COEFX(I,J) = ALPHA_C*OMEGA*(EXP((X_REL/R_MAX)**M)-1.0)	&
										/(EXP(1.0)-1.0)
			COEFY(I,J) = ALPHA_C*OMEGA*(EXP((Y_REL/R_MAX)**M)-1.0)	&
										/(EXP(1.0)-1.0)
			IF (COEFX(I,J).LE.EPS) COEFX(I,J) = 0.0
			IF (COEFY(I,J).LE.EPS) COEFY(I,J) = 0.0
    	 ENDDO
      ENDDO
      LO%SPONGE_COEFX = COEFX
	  LO%SPONGE_COEFY = COEFY

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE H_CALC (LO,H_MEAN,H_MAX)
!DESCRIPTION:
!     DETERMINE COEFFICIENT ALPHA SO AS TO IMPROVE NUMERICAL DISPERSION
!.....REF: WEI AND KIRBY (1995) AND KIRBY ET AL (1998)
!NOTES:
!	  #. CREATED ON ??? ?? 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!	  #. UPDATED ON FEB 25 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  REAL DX, DY, DT, G, R_MAX, RX, RY, X_REL, R, WIDTH
	  REAL H_MEAN,H_SUM
	  INTEGER K
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      K = 0
	  H_MEAN = ZERO
	  H_SUM = ZERO
	  H_MAX = ZERO

      DO I = 1,LO%NX
	     DO J = 1,LO%NY
            IF (LO%H(I,J) .GT. GX) THEN
               H_SUM = H_SUM + LO%H(I,J)
			   IF (H_MAX.LT.LO%H(I,J)) H_MAX = LO%H(I,J)
			   K = K+1
			ENDIF
    	 ENDDO
      ENDDO
      IF (K.GT.0) H_MEAN = H_SUM/K
!	  WRITE(*,*) H_MEAN

	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE BC_WALL (LO,WAVE_INFO)
!.....SET UP VERTICAL WALL BOUNDARY
!Incident direction( 1:top,2:bt,3:lf,4:rt,5:ob )
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE WAVE_PARAMS
	  TYPE (LAYER) :: LO
	  TYPE (WAVE)  :: WAVE_INFO
	  IF (LO%INI_SWITCH.EQ.2 ) THEN
	     IF (WAVE_INFO%INCIDENT.EQ.1) THEN
            LO%H(1:2,:) = -999.0
            LO%H(:,1:2) = -999.0
	        LO%H(LO%NX-1:LO%NX,:) = -999.0
!*	        LO%H(:,LO%NY-1:LO%NY) = -999.0
		 ENDIF
	     IF (WAVE_INFO%INCIDENT.EQ.2) THEN
            LO%H(1:2,:) = -999.0
!*            LO%H(:,1:2) = -999.0
	        LO%H(LO%NX-1:LO%NX,:) = -999.0
	        LO%H(:,LO%NY-1:LO%NY) = -999.0
		 ENDIF
	     IF (WAVE_INFO%INCIDENT.EQ.3) THEN
!*            LO%H(1:2,:) = -999.0
            LO%H(:,1:2) = -999.0
	        LO%H(LO%NX-1:LO%NX,:) = -999.0
	        LO%H(:,LO%NY-1:LO%NY) = -999.0
		 ENDIF
	     IF (WAVE_INFO%INCIDENT.EQ.4) THEN
            LO%H(1:2,:) = -999.0
            LO%H(:,1:2) = -999.0
!*	        LO%H(LO%NX-1:LO%NX,:) = -999.0
	        LO%H(:,LO%NY-1:LO%NY) = -999.0
		 ENDIF
	  ELSE
         LO%H(1:2,:) = -999.0
         LO%H(:,1:2) = -999.0
	     LO%H(LO%NX-1:LO%NX,:) = -999.0
	     LO%H(:,LO%NY-1:LO%NY) = -999.0
	  ENDIF
	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE READ_FACTS (BCI_INFO,LO,SWITCH)
!.....READ FACTS INPUT FILES AS INPUT BOUNDARY CONDITIONS
! SWITCH:
!       0-SURFACE ELEVATION DATA, ####H.ASC
!       1-HORIZONTAL VELOCITY DATA, ####U.ASC
!       2-VERTICAL VELOCITY DATA, ####V.ASC
!.....CREATED ON NOV 11 2008 (XIAOMING WANG, GNS)
!     LAST REVISE: NOV.27, 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE BCI_PARAMS
	  TYPE (LAYER)  :: LO
	  TYPE (BCI)    :: BCI_INFO
      INTEGER NX,NY,NT
	  INTEGER BC_TYPE,POS,SWITCH,STAT
	  REAL H(LO%NX,LO%NY)
	  REAL,ALLOCATABLE :: SNAPSHOT(:,:,:)
	  CHARACTER(LEN=120)   INE,LINE1,LINE2,LINE3
	  CHARACTER(LEN=120)   DUMP,TMP,FNAME
	  CHARACTER(LEN=18)   TEXT
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!	  WRITE(*,*) BCI_INFO%FNAMEH
!	  WRITE(*,*) BCI_INFO%FNAMEU
!	  WRITE(*,*) BCI_INFO%FNAMEV
      FNAME = ''
      IF (SWITCH.EQ.0) FNAME = BCI_INFO%FNAMEH
	  IF (SWITCH.EQ.1) FNAME = BCI_INFO%FNAMEU
	  IF (SWITCH.EQ.2) FNAME = BCI_INFO%FNAMEV
	  WRITE (*,*) FNAME
!
!...../////// READ H,U,V DATA FROM FACTS OUTPUT ////////////
	  OPEN(UNIT=1,FILE=FNAME,STATUS='OLD',IOSTAT=STAT,FORM='FORMATTED')
	  IF (STAT /=0) THEN
         PRINT *,"ERROR:: CAN'T OPEN FACTS DATA",FNAME,"; EXITING."
         STOP
      END IF
	  POS = 0
	  DO WHILE (POS.LE.0)
	     READ (1,*) DUMP
!         WRITE (*,*) DUMP
		 POS = INDEX(DUMP,'GEOMETRY')
		 DUMP = ''
!		 PAUSE
	  ENDDO
      READ (1,*) TEXT,NX,NY,NT
	  WRITE(*,*) NX,NY,NT
	  ALLOCATE(SNAPSHOT(NX,NY,NT))
	  SNAPSHOT = ZERO

      IF (SWITCH.EQ.0) THEN
	     ALLOCATE(BCI_INFO%X(NX))
	     ALLOCATE(BCI_INFO%Y(NY))
	     ALLOCATE(BCI_INFO%T(NT))
	     ALLOCATE(BCI_INFO%SNAPSHOT(NX,NY,NT))
		 ALLOCATE(BCI_INFO%SNAPSHOTU(NX,NY,NT))
		 ALLOCATE(BCI_INFO%SNAPSHOTV(NX,NY,NT))
	     ALLOCATE(BCI_INFO%Z_VERT(LO%NY,NT,2))
	     ALLOCATE(BCI_INFO%U_VERT(LO%NY,NT,2))
	     ALLOCATE(BCI_INFO%V_VERT(LO%NY,NT,2))
	     ALLOCATE(BCI_INFO%Z_HORI(LO%NX,NT,2))
	     ALLOCATE(BCI_INFO%U_HORI(LO%NX,NT,2))
	     ALLOCATE(BCI_INFO%V_HORI(LO%NX,NT,2))
      END IF

	  READ (1,*) DUMP
	  READ (1,*) DUMP
	  READ (1,*) DUMP
	  READ (1,*) (BCI_INFO%X(I),I=1,NX)
	  READ (1,*) DUMP
	  READ (1,*) (BCI_INFO%Y(I),I=1,NY)
	  READ (1,*) DUMP
	  READ (1,*) (BCI_INFO%T(I),I=1,NT)
	  READ (1,*) DUMP

	  DO K = 1,NT
	     DO J = 1,NY
		    READ (1,*) (SNAPSHOT(I,J,K),I=1,NX)
		 ENDDO
	  ENDDO

	  CLOSE(1)

	  IF (SWITCH.EQ.0) THEN
	     BCI_INFO%NX = NX
	     BCI_INFO%NY = NY
	     BCI_INFO%NT = NT
         BCI_INFO%DURATION = BCI_INFO%T(NT)-BCI_INFO%T(1)
!*	     BCI_INFO%T(:) = BCI_INFO%T(:)-BCI_INFO%T(1)
	  ENDIF
!.....CONVERT UNITS FROM CM OR CM/S TO M OR M/S
	  IF (SWITCH.EQ.0) BCI_INFO%SNAPSHOT = SNAPSHOT/100.0
	  IF (SWITCH.EQ.1) BCI_INFO%SNAPSHOTU = SNAPSHOT/100.0
	  IF (SWITCH.EQ.2) BCI_INFO%SNAPSHOTV = SNAPSHOT/100.0

!	  WRITE(*,*) NX,NY,NT

	  DEALLOCATE(SNAPSHOT,STAT = ISTAT)

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE  GET_BC_DATA (BI,LO)
!......................................................................
!DESCRIPTION:
!	  #. MAP THE BOUNDARY CONDITION FROM FACTS TO COMPUTATIONAL GRIDS;
!NOTES:
!	  #. CREATED ON NOV 11 2008 (XIAOMING WANG, GNS)
!     #. LAST REVISE: NOV.20, 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON MAR 10 2009 (XIAOMING WANG, GNS)
!		 1. ADD SUPPORT ON UTM COORDINATES
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE BCI_PARAMS
	  TYPE (LAYER) :: LO
	  TYPE (BCI)   :: BI
	  INTEGER BC_TYPE,POS
	  REAL TEMPX(LO%NX,2),TEMPY(2,LO%NY)
	  REAL SH
!*	  CHARACTER(LEN=120) 	:: LINE,LINE1,LINE2
!*	  CHARACTER(LEN=120) 	:: TMP,FNAMEH,FNAMEU,FNAMEV
!*	  CHARACTER(LEN=18)  TEXT
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  WRITE(*,*) 'READING PARAMETERS FOR FACTS INPUT BOUNDARIES...'
      CALL READ_FACTS (BI,LO,0)
	  CALL READ_FACTS (BI,LO,1)
	  CALL READ_FACTS (BI,LO,2)
	  WRITE(*,*) 'PROCESSING DATA FOR FACTS INPUT BOUNDARIES...'
	  IF (LO%LAYGOV.EQ.0) THEN
	     SH = LO%DX/60.0*0.5
	  ELSE
	     SH = 0.5*LO%DX/(RAD_MIN*R_EARTH*COS(LO%YO*RAD_DEG))/60.0
	  ENDIF

!.....OBTAIN VALUES OF Z,HU,HV ALONG THE BOUNDARIES OF NUMERICAL DOMAIN
      DO K = 1,BI%NT
	  	 WRITE (*,*) '   PROCESSING FACTS SNAPSHOT: ',K,' OUT OF ',BI%NT
!IF LAYER LO USES SPHERICAL COORDINATES //
		 IF (LO%LAYCORD.EQ.0) THEN
!...........OBTAIN VALUES OF Z AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X,LO%Y(1:2),LO%NX,2,				&
						BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
		       BI%Z_HORI(I,K,1) = TEMPX(I,1)
!			   WRITE (*,*) BI%Z_HORI(I,K,1)
		    ENDDO
		    !TOP BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X,LO%Y(LO%NY-1:LO%NY),LO%NX,2,	&
						BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
               BI%Z_HORI(I,K,2) = TEMPX(I,2)
		    ENDDO
		    !LEFT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(1:2),LO%Y,2,LO%NY,				&
						BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%Z_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(LO%NX-1:LO%NX),LO%Y,2,LO%NY,	&
						BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%Z_VERT(J,K,2) = TEMPY(2,J)
		    ENDDO

!...........OBTAIN VALUES OF U AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X+SH,LO%Y(1:2),LO%NX,2,			&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
		       BI%U_HORI(I,K,1) = TEMPX(I,1)
		    ENDDO
		    !TOP BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X+SH,LO%Y(LO%NY-1:LO%NY),LO%NX,	&
							2,BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
               BI%U_HORI(I,K,2) = TEMPX(I,2)
		    ENDDO
		    !LEFT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(1:2)+SH,LO%Y,2,LO%NY,			&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%U_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(LO%NX-1:LO%NX)+SH,				&
							LO%Y,2,LO%NY,BI%SNAPSHOTU(:,:,K),			&
							BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%U_VERT(J,K,2) = TEMPY(1,J)
		    ENDDO

!...........OBTAIN VALUES OF V AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X,LO%Y(1:2)+SH,LO%NX,2,			&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
		       BI%V_HORI(I,K,1) = TEMPX(I,1)
		    ENDDO
		    !TOP BOUNDARY
		    CALL GRID_INTERP (TEMPX,LO%X,LO%Y(LO%NY-1:LO%NY)+SH,LO%NX,	&
							2,BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO I = 1,LO%NX
               BI%V_HORI(I,K,2) = TEMPX(I,1)
		    ENDDO
		    !LEFT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(1:2),LO%Y+SH,2,LO%NY,			&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%V_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    CALL GRID_INTERP (TEMPY,LO%X(LO%NX-1:LO%NX),LO%Y+SH,		&
							2,LO%NY,BI%SNAPSHOTV(:,:,K),				&
							BI%X,BI%Y,BI%NX,BI%NY)
		    DO J = 1,LO%NY
		       BI%V_VERT(J,K,2) = TEMPY(2,J)
		    ENDDO
		 ENDIF
!IF LAYER LO USES UTM COORDINATES ///
		 IF (LO%LAYCORD.EQ.1) THEN
!...........OBTAIN VALUES OF Z AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,1),LO%CXY(I,1,1),				&
										LO%CXY(I,1,2),1,1,				&
							BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%Z_HORI(I,K,1) = TEMPX(I,1)
		    ENDDO
		    !TOP BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,2),LO%CXY(I,LO%NY,1),			&
									LO%CXY(I,LO%NY,2),1,1,				&
							BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%Z_HORI(I,K,2) = TEMPX(I,2)
		    ENDDO
		    !LEFT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(1,J),LO%CXY(1,J,1),				&
									LO%CXY(1,J,2),1,1,					&
							BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%Z_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(2,J),LO%CXY(LO%NX,J,1),			&
									LO%CXY(LO%NX,J,2),1,1,				&
							BI%SNAPSHOT(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%Z_VERT(J,K,2) = TEMPY(2,J)
		    ENDDO

!...........OBTAIN VALUES OF U AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,1),LO%CXY(I,1,1)+SH,			&
										LO%CXY(I,1,2),1,1,				&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%U_HORI(I,K,1) = TEMPX(I,1)
		    ENDDO
		    !TOP BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,2),LO%CXY(I,LO%NY,1)+SH,		&
									LO%CXY(I,LO%NY,2),1,1,				&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%U_HORI(I,K,2) = TEMPX(I,2)
		    ENDDO
		    !LEFT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(1,J),LO%CXY(1,J,1)+SH,			&
									LO%CXY(1,J,2),1,1,					&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%U_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(2,J),LO%CXY(LO%NX-1,J,1)+SH,		&
									LO%CXY(LO%NX-1,J,2),1,1,				&
							BI%SNAPSHOTU(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%U_VERT(J,K,2) = TEMPY(2,J)
		    ENDDO

!...........OBTAIN VALUES OF V AT BOUNDARIES
		    TEMPX = 0.0
		    TEMPY = 0.0
		    !BOTTOM BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,1),LO%CXY(I,1,1),				&
										LO%CXY(I,1,2)+SH,1,1,			&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%V_HORI(I,K,1) = TEMPX(I,1)
		    ENDDO
		    !TOP BOUNDARY
		    DO I = 1,LO%NX
		       CALL GRID_INTERP (TEMPX(I,2),LO%CXY(I,LO%NY-1,1),			&
									LO%CXY(I,LO%NY-1,2)+SH,1,1,			&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%V_HORI(I,K,2) = TEMPX(I,2)
		    ENDDO
		   !LEFT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(1,J),LO%CXY(1,J,1),				&
									LO%CXY(1,J,2)+SH,1,1,				&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%V_VERT(J,K,1) = TEMPY(1,J)
		    ENDDO
		    !RIGHT BOUNDARY
		    DO J = 1,LO%NY
		 	   CALL GRID_INTERP (TEMPY(2,J),LO%CXY(LO%NX,J,1),			&
									LO%CXY(LO%NX,J,2)+SH,1,1,			&
							BI%SNAPSHOTV(:,:,K),BI%X,BI%Y,BI%NX,BI%NY)
		       BI%V_VERT(J,K,2) = TEMPY(2,J)
		    ENDDO
		 ENDIF

!         WRITE (*,*) (BI%Z_VERT(I,K,1),I = 1,LO%NY)
!OBTAIN FLUXES AND ADJUST VALUES ON LAND
		 J = 1	!BOTTOM BOUNDARY
	     DO I = 1,LO%NX
	        IF (LO%H(I,J).LE.GX) THEN
		       BI%Z_HORI(I,K,1) = ZERO
			   BI%U_HORI(I,K,1) = ZERO
			   BI%V_HORI(I,K,1) = ZERO
			ELSE
			   BI%U_HORI(I,K,1) = BI%U_HORI(I,K,1)*LO%HP(I,J)
			   BI%V_HORI(I,K,1) = BI%V_HORI(I,K,1)*LO%HQ(I,J)
			ENDIF
		 ENDDO
		 J = LO%NY	!TOP BOUNDARY
		 DO I = 1,LO%NX
	        IF (LO%H(I,J).LE.GX) THEN
		       BI%Z_HORI(I,K,2) = ZERO
			   BI%U_HORI(I,K,2) = ZERO
			   BI%V_HORI(I,K,2) = ZERO
			ELSE
			   BI%U_HORI(I,K,2) = BI%U_HORI(I,K,2)*LO%HP(I,J)
			   BI%V_HORI(I,K,2) = BI%V_HORI(I,K,2)*LO%HQ(I,J)
			ENDIF
	     ENDDO

		 I = 1	!LEFT BOUNDARY
	     DO J = 1,LO%NY
	        IF (LO%H(I,J).LE.GX) THEN
		       BI%Z_VERT(J,K,1) = ZERO
			   BI%U_VERT(J,K,1) = ZERO
			   BI%V_VERT(J,K,1) = ZERO
			ELSE
			   BI%U_VERT(J,K,1) = BI%U_VERT(J,K,1)*LO%HP(I,J)
			   BI%V_VERT(J,K,1) = BI%V_VERT(J,K,1)*LO%HQ(I,J)
			ENDIF
		 ENDDO
		 I = LO%NX	!RIGHT BOUNDARY
		 DO J = 1,LO%NY
	        IF (LO%H(I,J).LE.GX) THEN
		       BI%Z_VERT(J,K,2) = ZERO
			   BI%U_VERT(J,K,2) = ZERO
			   BI%V_VERT(J,K,2) = ZERO
			ELSE
			   BI%U_VERT(J,K,2) = BI%U_VERT(J,K,2)*LO%HP(I,J)
			   BI%V_VERT(J,K,2) = BI%V_VERT(J,K,2)*LO%HQ(I,J)
			ENDIF
	     ENDDO

      ENDDO

!	  WRITE (*,*) LO%CXY(1,1,1),LO%CXY(1,1,2)
!	  WRITE (*,*) LO%CXY(LO%NX,LO%NY,1),LO%CXY(LO%NX,LO%NY,2)
!.....DEALLOCATE VARIABLE NO LONGER BEING USED.
      DEALLOCATE(BI%SNAPSHOT,BI%SNAPSHOTU,BI%SNAPSHOTV,STAT=ISTAT)

	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE BC_INPUT (BI,LO,TIME)
!......................................................................
!DESCRIPTION:
!	  #. MAP THE BOUNDARY CONDITION FROM FACTS TO COMPUTATIONAL GRIDS;
!NOTES:
!	  #. CREATED ON NOV 11 2008 (XIAOMING WANG, GNS)
!     #. LAST REVISE: NOV.24, 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON MAR 10 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE BCI_PARAMS
	  TYPE (LAYER) :: LO
	  TYPE (BCI)   :: BI
	  REAL T(BI%NT),TMPX(LO%NX),TMPY(LO%NY)
	  REAL TIME,C1,C2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      T(:) = BI%T(:)

!.....OBTAIN SURFACE ELEVATION AT T = TIME FROM FACTS DATA
	  KT = 1
      DO K = 1,BI%NT-1
         IF (TIME-0.5*LO%DT.GE.BI%T(K) .AND.						&
				TIME-0.5*LO%DT.LT.BI%T(K+1)) KT = K
	  ENDDO
	  C1 = (TIME-0.5*LO%DT-BI%T(KT))/(BI%T(KT+1)-BI%T(KT))
	  C2 = 1.0 - C1

      DO I = 1,LO%NX
		 LO%Z(I,1,2) = C1*BI%Z_HORI(I,KT+1,1) + C2*BI%Z_HORI(I,KT,1)
	     LO%Z(I,LO%NY,2) = C1*BI%Z_HORI(I,KT+1,2) + C2*BI%Z_HORI(I,KT,2)
!		 DZDT = (BI%Z_HORI(I,KT+1,1)-BI%Z_HORI(I,KT,1))/(T(KT+1)-T(KT))
!	     LO%Z(I,1,2) = LO%Z(I,1,2) + LO%DT*DZDT
!		 DZDT = (BI%Z_HORI(I,KT+1,2)-BI%Z_HORI(I,KT,2))/(T(KT+1)-T(KT))
!	     LO%Z(I,LO%NY,2) = LO%Z(I,LO%NY,2) + LO%DT*DZDT
	  ENDDO
	  DO J = 1,LO%NY
	     LO%Z(1,J,2) = C1*BI%Z_VERT(J,KT+1,1) + C2*BI%Z_VERT(J,KT,1)
	     LO%Z(LO%NX,J,2) = C1*BI%Z_VERT(J,KT+1,2) + C2*BI%Z_VERT(J,KT,2)
!	     DZDT = (BI%Z_VERT(J,KT+1,1)-BI%Z_VERT(J,KT,1))/(T(KT+1)-T(KT))
!	     LO%Z(1,J,2) = LO%Z(1,J,2) + LO%DT*DZDT
!	     DZDT = (BI%Z_VERT(J,KT+1,2)-BI%Z_VERT(J,KT,2))/(T(KT+1)-T(KT))
!	     LO%Z(LO%NX,J,2) = LO%Z(LO%NX,J,2) + LO%DT*DZDT
	  ENDDO

!.....OBTAIN VOLUME FLUX AT T = TIME
	  KT = 1
	  DO K = 1,BI%NT-1
         IF (TIME.GE.BI%T(K) .AND. TIME.LT.BI%T(K+1)) KT = K
	  ENDDO
	  C1 = (TIME-BI%T(KT))/(BI%T(KT+1)-BI%T(KT))
	  C2 = 1.0 - C1

	  DO I = 1,LO%NX
		 LO%M(I,1,2) = C1*BI%U_HORI(I,KT+1,1) + C2*BI%U_HORI(I,KT,1)
	     LO%M(I,LO%NY,2) = C1*BI%U_HORI(I,KT+1,2) + C2*BI%U_HORI(I,KT,2)
!		 DPDT = (BI%U_HORI(I,KT+1,1)-BI%U_HORI(I,KT,1))/(T(KT+1)-T(KT))
!	     LO%M(I,1,1) = LO%M(I,1,1) + LO%DT*DPDT
!		 DPDT = (BI%U_HORI(I,KT+1,2)-BI%U_HORI(I,KT,2))/(T(KT+1)-T(KT))
!	     LO%M(I,LO%NY,1) = LO%M(I,LO%NY,1) + LO%DT*DPDT
	  ENDDO
	  DO J = 1,LO%NY
	     LO%M(1,J,2) = C1*BI%U_VERT(J,KT+1,1) + C2*BI%U_VERT(J,KT,1)
	     LO%M(LO%NX-1,J,2) = C1*BI%U_VERT(J,KT+1,2) + C2*BI%U_VERT(J,KT,2)
!	     DPDT = (BI%U_VERT(J,KT+1,1)-BI%U_VERT(J,KT,1))/(T(KT+1)-T(KT))
!	     LO%M(1,J,1) = LO%M(1,J,1) + LO%DT*DPDT
!	     DPDT = (BI%U_VERT(J,KT+1,2)-BI%U_VERT(J,KT,2))/(T(KT+1)-T(KT))
!	     LO%M(LO%NX-1,J,1) = LO%M(LO%NX-1,J,1) + LO%DT*DPDT
	  ENDDO
	  DO I = 1,LO%NX
		 LO%N(I,1,2) = C1*BI%V_HORI(I,KT+1,1) + C2*BI%V_HORI(I,KT,1)
	     LO%N(I,LO%NY-1,2) = C1*BI%V_HORI(I,KT+1,2) + C2*BI%V_HORI(I,KT,2)
!		 DQDT = (BI%V_HORI(I,KT+1,1)-BI%V_HORI(I,KT,1))/(T(KT+1)-T(KT))
!	     LO%N(I,1,1) = LO%N(I,1,1) + LO%DT*DQDT
!         DQDT = (BI%V_HORI(I,KT+1,2)-BI%V_HORI(I,KT,2))/(T(KT+1)-T(KT))
!	     LO%N(I,LO%NY-1,1) = LO%N(I,LO%NY-1,1) + LO%DT*DQDT
	  ENDDO
	  DO J = 1,LO%NY
	     LO%N(1,J,2) = C1*BI%V_VERT(J,KT+1,1) + C2*BI%V_VERT(J,KT,1)
	     LO%N(LO%NX,J,2) = C1*BI%V_VERT(J,KT+1,2) + C2*BI%V_VERT(J,KT,2)
!	     DQDT = (BI%V_VERT(J,KT+1,1)-BI%V_VERT(J,KT,1))/(T(KT+1)-T(KT))
!	     LO%N(1,J,1) = LO%N(1,J,1) + LO%DT*DQDT
!	     DQDT = (BI%V_VERT(J,KT+1,2)-BI%V_VERT(J,KT,2))/(T(KT+1)-T(KT))
!	     LO%N(LO%NX,J,1) = LO%N(LO%NX,J,1) + LO%DT*DQDT
	  ENDDO

	  RETURN
	  END
