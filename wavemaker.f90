
!----------------------------------------------------------------------
      SUBROUTINE WAVE_MAKER (TIME,LO,WV)
!.....SEND IN INCIDENT WAVES THROUGH A BOUNDARY
!.....CREATED BY XIAOMING WANG (MAR 15, 2004)
!     UPDATED BY XIAOMING WANG (SEP 17 2006)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      USE WAVE_PARAMS
!      IMPLICIT NONE
      TYPE (LAYER)  :: LO
      TYPE (WAVE)   :: WV
	  REAL ETA,FLUX,TIME,X,DX,SOUTH_LAT,DX_RAD,DD,CC,TIME_LAG,T
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA OSIXTY/0.016666666667/, BIG/-999./
      DATA RAD/0.01745329252/
!      DATA RR/6.37E6/           !      RADIUS OF EARTH
!      DATA PI/3.14159265358979/


	  IF (LO%LAYCORD .EQ. 0) THEN
	     SOUTH_LAT = LO%SOUTH_LAT*RAD_DEG
	     DX_RAD = LO%DX*RAD_MIN
         DX = R_EARTH*COS(SOUTH_LAT)*DX_RAD   !CONVERT RAD TO METERS
      ELSE
		 DX = LO%DX
      ENDIF
	  CC = SQRT(GRAV*(WV%DEPTH+WV%AMP))

      IF (WV%MK_TYPE .EQ. 5) THEN
!	     TIME = TIME-19*LO%DT
         DO I = 1,LO%NX
		    X = I*LO%DX - 50.0
			IF (X.LE.0.0) THEN
			   CALL SOLIT (LO,ETA,FLUX,X,TIME-LO%DT,WV)
		       DO J = 1,LO%NY
				  LO%Z(I,J,1) = ETA
				  LO%M(I,J,1) = FLUX
			   ENDDO
			   CALL SOLIT (LO,ETA,FLUX,X,TIME,WV)
			   DO J = 1,LO%NY
				  LO%Z(I,J,2) = ETA
				  LO%M(I,J,2) = FLUX
			   ENDDO
            ENDIF
		 ENDDO
	  ENDIF

	  X = 0.0
!      CALL SOLIT (ETA,FLUX,X,TIME,WV)
      IF (WV%MK_TYPE .LE. 2) THEN ! FOR SOLITARY OR TIMEHISTORY INPUT
	     CALL SOLIT (LO,ETA,FLUX,X,TIME,WV)
	     SELECT CASE (WV%INCIDENT)
	     CASE (1)     !FROM TOP BOUNDARY
	        DO I = 1,LO%NX
	           LO%Z(I,LO%NY,2) = ETA !+LO%Z(I,LO%NY,1)
	           LO%N(I,LO%NY,2) = -FLUX !+LO%N(I,LO%NY,1)
	        ENDDO
	     CASE (2)    !FROM BOTTOM BOUNDARY
	        DO I = 1,LO%NX
	           LO%Z(I,1,2) = ETA !+LO%Z(I,1,1)
	           LO%N(I,1,2) = FLUX !+LO%N(I,1,1)
			ENDDO
	     CASE (3)   !FROM LEFT BOUNDARY
	        DO J = 1,LO%NY
	           LO%Z(1,J,2) = ETA !+LO%Z(1,J,1)
	           LO%M(1,J,2) = FLUX !+LO%M(1,J,1)
			ENDDO
	     CASE (4)   !FROM RIGHT BOUNDARY
	        DO J = 1,LO%NY
	           LO%Z(LO%NX,J,2) = ETA !+LO%Z(LO%NX,J,1)
	           LO%M(LO%NX,J,2) = -FLUX !+LO%M(LO%NX,J,1)
			ENDDO
         CASE (5)  !OBLIQUE INCIDENT WAVE
			!PROPAGATE TO UPPER-RIGHT
		    IF (WV%ANG.GE.0.0 .AND. WV%ANG.LT.90.0) THEN
			   ANG = WV%ANG*RAD_DEG
			   SN = SIN(ANG)
			   CS = COS(ANG)
               DO J=1,LO%NY
                  DIS = DBLE(J-1)*DX*CS
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(1,J,2) = ETA
				  LO%M(1,J,2) = FLUX*SN
				  LO%N(1,J,2) = FLUX*CS
			   ENDDO
			   DO I=1,LO%NX
                  DIS = DBLE(I-1)*DX*SN
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(I,1,2) = ETA
				  LO%M(I,1,2) = FLUX*SN
				  LO%N(I,1,2) = FLUX*CS
			   ENDDO
			ENDIF
			!PROPAGATE TO LOWER-RIGHT
		    IF (WV%ANG.GE.90.0 .AND. WV%ANG.LT.180.0) THEN
			   ANG = (180.0-WV%ANG)*RAD_DEG
			   CS = COS(ANG)
			   SN = SIN(ANG)
               DO J=1,LO%NY
                  DIS = DBLE(LO%NY-J)*DX*CS
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(1,J,2) = ETA
				  LO%M(1,J,2) = FLUX*SN
				  LO%N(1,J,2) = -FLUX*CS
			   ENDDO
			   DO I=1,LO%NX
                  DIS = DBLE(I-1)*DX*SN
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(I,LO%NY,2) = ETA
				  LO%M(I,LO%NY,2) = FLUX*SN
				  LO%N(I,LO%NY,2) = -FLUX*CS
			   ENDDO
			ENDIF
			!PROPAGATE TO LOWER-LEFT
		    IF (WV%ANG.GE.180.0 .AND. WV%ANG.LT.270.0) THEN
			   ANG = (270.0-WV%ANG)*RAD_DEG
			   CS = COS(ANG)
			   SN = SIN(ANG)
               DO J=1,LO%NY
                  DIS = DBLE(LO%NY-J)*DX*SN
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(LO%NX,J,2) = ETA
				  LO%M(LO%NX,J,2) = -FLUX*CS
				  LO%N(LO%NX,J,2) = -FLUX*SN
			   ENDDO
			   DO I=1,LO%NX
                  DIS = DBLE(LO%NX-I)*DX*CS
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(I,LO%NY,2) = ETA
				  LO%M(I,LO%NY,2) = -FLUX*CS
				  LO%N(I,LO%NY,2) = -FLUX*SN
			   ENDDO
			ENDIF
			!PROPAGATE TO UPPER-LEFT
		    IF (WV%ANG.GE.270.0 .AND. WV%ANG.LT.360.0) THEN
			   ANG = (360.0-WV%ANG)*RAD_DEG
			   CS = COS(ANG)
			   SN = SIN(ANG)
               DO J=1,LO%NY
                  DIS = DBLE(J-1)*DX*CS
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(LO%NX,J,2) = ETA
				  LO%M(LO%NX,J,2) = -FLUX*SN
				  LO%N(LO%NX,J,2) = FLUX*CS
			   ENDDO
			   DO I=1,LO%NX
                  DIS = DBLE(LO%NX-I)*DX*SN
				  TIME_LAG = DIS/CC
				  T = TIME - TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
				  LO%Z(I,1,2) = ETA
				  LO%M(I,1,2) = -FLUX*SN
				  LO%N(I,1,2) = FLUX*CS
			   ENDDO
			ENDIF
         END SELECT
	  ENDIF
	  IF (WV%MK_TYPE .EQ. 3) THEN ! FOCUSING SOLITARY WAVE
		 SELECT CASE (WV%INCIDENT)
			CASE (1)  !WAVE FROM TOP BOUNARY
!*		       D0 = DBLE(LO%NY-1)*DX-WV%POINT(2)
		       D0 = LO%Y(LO%NY)-WV%POINT(2)
			   DO I = 1,LO%NX
			      DIS_VERT = D0
!*				  DIS_HORI = DBLE(I-1)*DX-WV%POINT(1)
				  DIS_HORI = LO%X(I)-WV%POINT(1)
			      DD = SQRT(DIS_VERT**2+DIS_HORI**2)-D0
                  CC = SQRT(GRAV*(WV%DEPTH+WV%AMP))
                  TIME_LAG = DD/CC
				  T = TIME + TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
	              LO%Z(I,LO%NY,2) = ETA !+LO%Z(1,J,1)
	              LO%N(I,LO%NY,2) = -FLUX*DIS_VERT/(DD+D0) !+LO%M(1,J,1)
				  LO%M(I,LO%NY,2) = -FLUX*DIS_HORI/(DD+D0)
			   ENDDO
		    CASE (2)  ! WAVE FROM BOTTOM BOUNDARY
		       D0 = WV%POINT(2)-LO%Y(1)
			   DO I = 1,LO%NX
			      DIS_VERT = D0
!*				  DIS_HORI = DBLE(I-1)*DX-WV%POINT(1)
				  DIS_HORI =  LO%X(I)-WV%POINT(1)
			      DD = SQRT(DIS_HORI**2+DIS_VERT**2)-D0
                  !CC = SQRT(9.807*(WV%DEPTH+WV%AMP))
                  TIME_LAG = DD/CC
				  T = TIME + TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
	              LO%Z(I,1,2) = ETA !+LO%Z(1,J,1)
	              LO%N(I,1,2) = FLUX*DIS_VERT/(DD+D0) !+LO%M(1,J,1)
                  LO%M(I,1,2) = -FLUX*DIS_HORI/(DD+D0)
			   ENDDO
		    CASE (3)  !WAVE FROM LEFT BOUNDARY
!*		       D0 = WV%POINT(1)
		       D0 = WV%POINT(1) - LO%X(1)
			   DO J = 1,LO%NY
			      DIS_HORI = D0
!*				  DIS_VERT = DBLE(J-1)*DX-WV%POINT(2)
				  DIS_VERT = LO%Y(J)-WV%POINT(2)
			      DD = SQRT(DIS_VERT**2+DIS_HORI**2)-D0
                  !CC = SQRT(9.807*(WV%DEPTH+WV%AMP))
                  TIME_LAG = DD/CC
				  T = TIME + TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
	              LO%Z(1,J,2) = ETA !+LO%Z(1,J,1)
	              LO%M(1,J,2) = FLUX*DIS_HORI/(DD+D0) !+LO%M(1,J,1)
				  LO%N(1,J,2) = -FLUX*DIS_VERT/(DD+D0)
			   ENDDO
			CASE (4) !WAVE FROM RIGHT BOUNDARY
!*		       D0 = DBLE(LO%NX-1)*DX-WV%POINT(1)
		       D0 = LO%X(LO%NX)-WV%POINT(1)
			   DO J = 1,LO%NY
			      DIS_HORI = D0
!*				  DIS_VERT = DBLE(J-1)*DX-WV%POINT(2)
				  DIS_VERT = LO%Y(J)-WV%POINT(2)
			      DD = SQRT(DIS_VERT**2+DIS_HORI**2)-D0
                  !CC = SQRT(9.807*(WV%DEPTH+WV%AMP))
                  TIME_LAG = DD/CC
				  T = TIME + TIME_LAG
				  CALL SOLIT(LO,ETA,FLUX,X,T,WV)
	              LO%Z(LO%NX,J,2) = ETA !+LO%Z(1,J,1)
	              LO%M(LO%NX,J,2) = -FLUX*DIS_HORI/(DD+D0) !+LO%M(1,J,1)
				  LO%N(LO%NX,J,2) = -FLUX*DIS_VERT/(DD+D0)
			   ENDDO
         END SELECT
      ENDIF
      ! WAVE MAKER STOP WORKING ON LAND
      DO J=1,LO%NY
	     IF (LO%H(1,J) .LE. 0.0) THEN
	        LO%Z(1,J,2) = 0.0
		    LO%M(1,J,2) = 0.0
		    LO%N(1,J,2) = 0.0
		 ENDIF
		 IF (LO%H(LO%NX,J) .LE. 0.0) THEN
	        LO%Z(LO%NX,J,2) = 0.0
		    LO%M(LO%NX,J,2) = 0.0
		    LO%N(LO%NX,J,2) = 0.0
		 ENDIF
	  ENDDO
      DO I=1,LO%NX
	     IF (LO%H(I,1) .LE. 0.0) THEN
	        LO%Z(I,1,2) = 0.0
		    LO%M(I,1,2) = 0.0
		    LO%N(I,1,2) = 0.0
		 ENDIF
		 IF (LO%H(I,LO%NY) .LE. 0.0) THEN
	        LO%Z(I,LO%NY,2) = 0.0
		    LO%M(I,LO%NY,2) = 0.0
		    LO%N(I,LO%NY,2) = 0.0
		 ENDIF
	  ENDDO

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE SOLIT (LO,ETA,FLUX,X,TIME,WAVE_INFO)
!     CREATED ON SEP 17, 2006 (XIAOMING WANG)
!     ADDITIONAL PASSING PARAMETER, LO, IS ADDED ON MAR 18 2008
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      USE WAVE_PARAMS
	  TYPE (LAYER):: LO
	  TYPE (WAVE) :: WAVE_INFO
      REAL ETA, FLUX, TIME, X
	  REAL TIMELAG,THI,CE,WAVEPERI,WLENGTH
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!      PI = 3.14159265358979
	  G = GRAV
	  IF (LO%LAYCORD .EQ. 0) THEN
	     DX = LO%DX*RAD_MIN*R_EARTH
	  ELSE
	     DX = LO%DX
	  ENDIF
	  DT = LO%DT

      IF (WAVE_INFO%MK_TYPE.EQ.0) THEN
	     A0 = WAVE_INFO%AMP
		 C0 = SQRT(GRAV*(WAVE_INFO%AMP+WAVE_INFO%DEPTH))
		 PERI = 1.30
		 WLENGTH = C0*PERI
		 WNUM = 2.*PI/WLENGTH
		 OMIGA = 2.*PI/PERI

         ETA = A0*COS(WNUM*X-OMIGA*TIME+PI/2.0)
		 ETAP = A0*COS(WNUM*(X+LO%DX/2.)-OMIGA*(TIME+DT/2.0)+PI/2.0)
!	     CE = SQRT(9.807*(WAVE_INFO%DEPTH+ETA))
		 FLUX = ETAP*C0
      ENDIF

      IF (WAVE_INFO%MK_TYPE.EQ.5) THEN
	     A0 = WAVE_INFO%AMP
		 C0 = SQRT(GRAV*(WAVE_INFO%AMP+WAVE_INFO%DEPTH))
		 PERI = 1.30
		 WLENGTH = C0*PERI
		 WNUM = 2.0*PI/WLENGTH
		 OMIGA = 2.0*PI/PERI

         ETA = A0*COS(WNUM*X-OMIGA*TIME+PI/2.0)
		 ETAP = A0*COS(WNUM*(X+LO%DX/2.)-OMIGA*(TIME+DT/2.0)+PI/2.0)
!	     CE = SQRT(9.807*(WAVE_INFO%DEPTH+ETA))
		 FLUX = ETAP*C0
      ENDIF

	  IF (WAVE_INFO%MK_TYPE.EQ.1 .OR. WAVE_INFO%MK_TYPE.EQ.3) THEN
         THI = SQRT(3.0*WAVE_INFO%AMP/(4.0*WAVE_INFO%DEPTH**3))
         CE = SQRT(GRAV*(WAVE_INFO%DEPTH+WAVE_INFO%AMP))
         WAVEPERI = 2.0/(THI*CE)*(3.8+WAVE_INFO%AMP/WAVE_INFO%DEPTH)
         WLENGTH  = 2.0*2.12*WAVE_INFO%DEPTH						&
							/SQRT(WAVE_INFO%AMP/WAVE_INFO%DEPTH)
         TIMELAG  = 0.6 * WAVEPERI
         ETA = WAVE_INFO%AMP/COSH(SQRT(3.0/4.0*WAVE_INFO%AMP		&
						/WAVE_INFO%DEPTH**3)*CE*(TIMELAG-TIME))**2
!         PRINT *,ETA
		 FLUX = ETA*CE  !/(ETA+WAVE_INFO%DEPTH)
	  ELSEIF (WAVE_INFO%MK_TYPE.EQ.2) THEN
         DO I=1,WAVE_INFO%FORM_LEN-1
		  IF (TIME.LT.WAVE_INFO%T(1)) THEN
		     ETA = 0.0
			 FLUX = 0.0
          ELSEIF (TIME.GE.WAVE_INFO%T(I) .AND.						&
					TIME.LT.WAVE_INFO%T(I+1)) THEN
	         ETA=(WAVE_INFO%FSE(I+1)-WAVE_INFO%FSE(I))				&
					/(WAVE_INFO%T(I+1) - WAVE_INFO%T(I))			&
					*(TIME-WAVE_INFO%T(I))+WAVE_INFO%FSE(I)
             CE = SQRT(GRAV*(WAVE_INFO%DEPTH+ETA))
			 FLUX = ETA*CE  !/(ETA+WAVE_INFO%DEPTH)
		  ELSEIF (TIME.GE.WAVE_INFO%T(WAVE_INFO%FORM_LEN)) THEN
             ETA = 0.0
			 FLUX = 0.0
          ENDIF
         ENDDO
      ENDIF

	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE READ_WAVE (WAVE_INFO)
!.....CUSTOMIZED INPUT WAVE PROFILE
!     ONLY USED WHEN WAVE TYPE OPTION IS 2.
!----------------------------------------------------------------------
	  USE WAVE_PARAMS
	  TYPE (WAVE) :: WAVE_INFO
	  REAL H_MAX,SOUTH_LAT,DX,CR
	  REAL TEMP1,TEMP2
	  INTEGER COUNT
      CHARACTER(LEN=80) FNAME
      INTEGER :: RSTAT = 0

	  TEMP1 = 0.0
	  TEMP2 = 0.0
      IF (WAVE_INFO%MK_TYPE==2) THEN
         OPEN(UNIT=20,FILE=WAVE_INFO%FORM_NAME,STATUS='OLD',IOSTAT=ISTAT)
         IF (ISTAT /=0) THEN
            PRINT *,"ERROR:: CAN'T OPEN TIME HISTORY DATA; EXITING."
            STOP
         END IF
		 COUNT = -1
		 DO WHILE (RSTAT == 0)
		    COUNT = COUNT + 1
			READ (20,*,IOSTAT=RSTAT) TEMP1,TEMP2
	     ENDDO
!*		 CLOSE(20)
         WAVE_INFO%FORM_LEN = COUNT
		 ALLOCATE(WAVE_INFO%T(WAVE_INFO%FORM_LEN))
         ALLOCATE(WAVE_INFO%FSE(WAVE_INFO%FORM_LEN))
		 WAVE_INFO%T = 0.0
		 WAVE_INFO%FSE = 0.0
		 REWIND(20)
!*		 OPEN(UNIT=20,FILE=WAVE_INFO%FORM_NAME,STATUS='OLD',IOSTAT=ISTAT)
	     DO I=1,WAVE_INFO%FORM_LEN
	        READ(20,*) WAVE_INFO%T(I),WAVE_INFO%FSE(I)
	     ENDDO
		 CLOSE(20)
	  ENDIF
	  IF (WAVE_INFO%INCIDENT.EQ.5) THEN   ! FOR OBLIQUE WAVE
         WRITE (*,*) '    YOU ARE USING OBLIQUE INCIDENT WAVE. INCIDENT ANGLE IS MEASURED'
		 WRITE (*,*) '    CLOCKWISE FROM THE NORTH (UPWARD), RANGING 0.0 TO 360.'
		 WRITE (*,*) '    PLEASE INPUT INCDIENT ANGLE (IN DEGREES):'
		 READ *, WAVE_INFO%ANG
	  ENDIF
	  IF (WAVE_INFO%MK_TYPE.EQ.3) THEN   ! FOR FOCUSING INCIDENT WAVE
         WRITE (*,*) '    YOU ARE USING FOCUSING INCIDENT WAVE. THE FOCUS IS WHERE'
		 WRITE (*,*) '    THE INCIDENT WAVE FROM A BOUNDARY CONVERGES.'
		 WRITE (*,*) '    PLEASE INPUT X COORD. OF THE FOCUS (IN METERS):'
		 READ *, WAVE_INFO%POINT(1)
		 WRITE (*,*) '    PLEASE INPUT Y COORD. OF THE FOCUS (IN METERS):'
		 READ *, WAVE_INFO%POINT(2)
	  ENDIF


	  RETURN
	  END
