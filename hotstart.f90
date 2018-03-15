
!----------------------------------------------------------------------
      SUBROUTINE HOT_START (LO,LA,START_TYPE,START_STEP)
!     HOT START FUNCTION
!     START_TYPE:
!          = 0: SAVE RESUMING SNAPSHOTS
!          = 1: LOAD RESUMING DATA FROM ONLY FIRST-LEVEL GRIDS (LAYER1)
!          = 2: LOAD RESUMING DATA FROM ALL GRIDS
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER)  :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  INTEGER       :: START_TYPE,START_STEP,SWITCH
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  CALL HOTSTART_DATA (LO,START_TYPE,START_STEP)
	  IF (START_TYPE .NE. 1) THEN
	     DO I = 1,NUM_GRID
	        IF (LA(I)%LAYSWITCH .EQ. 0) CALL HOTSTART_DATA (LA(I),	&
											START_TYPE,START_STEP)
		 ENDDO
      ENDIF

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE HOTSTART_DATA (LO,START_TYPE,START_STEP)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER)  :: LO
	  INTEGER       :: START_TYPE,START_STEP,SWITCH
	  INTEGER IOSTAT
      CHARACTER(LEN=20) FNAME1,FNAME2,FNAME3,FNAME4
	  CHARACTER(LEN=20) FNAME5,FNAME6
	  CHARACTER(LEN=20) FNAME10,FNAME11
      
	  LAYID = LO%ID
	  IF ( LAYID .EQ. 1) THEN
	     IS = 1
	     IE = LO%NX
	     JS = 1
	     JE = LO%NY
	  ELSE
	     IS = 2
		 JS = 2
		 IE = LO%NX
		 JE = LO%NY
	  ENDIF

      WRITE (FNAME1,4) LAYID,START_STEP
 4    FORMAT('z1_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME2,5) LAYID,START_STEP
 5    FORMAT('Z2_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME3,6) LAYID,START_STEP
 6    FORMAT('m1_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME4,7) LAYID,START_STEP
 7    FORMAT('m2_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME5,8) LAYID,START_STEP
 8    FORMAT('n1_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME6,9) LAYID,START_STEP
 9    FORMAT('n2_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME10,10) LAYID,START_STEP
 10   FORMAT('zmax1_',I2.2,'_',I6.6,'.dat')
      WRITE (FNAME11,11) LAYID,START_STEP
 11   FORMAT('zmin1_',I2.2,'_',I6.6,'.dat')
!  ----- OUTPUT DATA FOR FUTURE HOT START -----------
!! OUTPUT WATER SURFACE DISPLACEMENT DATA Z
      IF (START_TYPE .EQ. 0) THEN
         OPEN (25,FILE=FNAME1,STATUS='UNKNOWN')
         DO J = JS,JE
            WRITE (25,'(10F12.5)') (LO%Z(I,J,1),I=IS,IE)
         ENDDO
         CLOSE (25)
!! OUTPUT FLUX DATA IN X DIRECTION
         OPEN (25,FILE=FNAME3,STATUS='UNKNOWN')
         DO J = JS,JE
            WRITE (25,'(10F12.5)') (LO%M(I,J,1),I=IS,IE)
         ENDDO
         CLOSE (25)
!! OUTPUT FLUX DATA IN Y DIRECTION
         OPEN (25,FILE=FNAME5,STATUS='UNKNOWN')
         DO J = JS,JE
           WRITE (25,'(10F12.5)') (LO%N(I,J,1),I=IS,IE)
         ENDDO
         CLOSE (25)
!! OUTPUT MAXIMUM SURFACE ELEVATION DATA
         OPEN (25,FILE=FNAME10,STATUS='UNKNOWN')
         DO J = JS,JE
           WRITE (25,'(10F12.5)') (LO%Z_MAX(I,J),I=IS,IE)
         ENDDO
         CLOSE (25)
!! OUTPUT MAXIMUM SURFACE DEPRESSION DATA
         OPEN (25,FILE=FNAME11,STATUS='UNKNOWN')
         DO J = JS,JE
           WRITE (25,'(10F12.5)') (LO%Z_MIN(I,J),I=IS,IE)
         ENDDO
         CLOSE (25)
	  ELSE
!------- LOAD DATA FOR HOT START ------------------------------------
!! LOAD WATER SURFACE DISPLACEMENT DATA FOR HOT START
         OPEN (UNIT=25,FILE=FNAME1,STATUS='OLD',IOSTAT=ISTAT,		&
												FORM='FORMATTED')
         IF (ISTAT /= 0) THEN
            PRINT *,"WARNING:: HOTSTART CAN'T OPEN Z FILE; USE ZEROS INSTEAD."
			LO%Z(:,:,1) = 0.0
         ELSE
            DO J = JS,JE
               READ (25,'(10F12.5)') (LO%Z(I,J,1),I=IS,IE)
            ENDDO
		 ENDIF
         CLOSE (25)
!! LOAD VOLUME FLUX DATA IN X DIRECTION FOR HOT START
		 OPEN (UNIT=25,FILE=FNAME3,STATUS='OLD',IOSTAT=ISTAT,		&
												FORM='FORMATTED')
         IF (ISTAT /= 0) THEN
            PRINT *,"WARNING:: HOTSTART CAN'T OPEN P FILE; USE ZEROS INSTEAD."
			LO%M(:,:,1) = 0.0
         ELSE
            DO J = JS,JE
               READ (25,'(10F12.5)') (LO%M(I,J,1),I=IS,IE)
            ENDDO
		 ENDIF
         CLOSE (25)
!! LOAD VOLUME FLUX DATA IN Y DIRECTION FOR HOT START
		 OPEN (UNIT=25,FILE=FNAME5,STATUS='OLD',IOSTAT=ISTAT,		&
												FORM='FORMATTED')
         IF (ISTAT /= 0) THEN
            PRINT *,"WARNING:: HOTSTART CAN'T OPEN Q FILE; USE ZEROS INSTEAD."
			LO%N(:,:,1) = 0.0
         ELSE
            DO J = JS,JE
               READ (25,'(10F12.5)') (LO%N(I,J,1),I=IS,IE)
            ENDDO
		 ENDIF
         CLOSE (25)
!! LOAD MAXIMUM SURFACE ELEVATION DATA FOR HOT START
		 OPEN (UNIT=25,FILE=FNAME10,STATUS='OLD',IOSTAT=ISTAT,		&
												FORM='FORMATTED')
         IF (ISTAT /= 0) THEN
            PRINT *,"WARNING:: HOTSTART CAN'T OPEN ZMAX FILE; USE ZEROS INSTEAD."
			LO%Z_MAX(:,:) = 0.0
         ELSE
            DO J = JS,JE
               READ (25,'(10F12.5)') (LO%Z_MAX(I,J),I=IS,IE)
            ENDDO
		 ENDIF
         CLOSE (25)
!! LOAD MAXIMUM SURFACE DEPRESSION DATA FOR HOT START
		 OPEN (UNIT=25,FILE=FNAME11,STATUS='OLD',IOSTAT=ISTAT,		&
												FORM='FORMATTED')
         IF (ISTAT /= 0) THEN
            PRINT *,"WARNING:: HOTSTART CAN'T OPEN ZMIN FILE; USE ZEROS INSTEAD."
			LO%Z_MIN(:,:) = 0.0
         ELSE
            DO J = JS,JE
               READ (25,'(10F12.5)') (LO%Z_MIN(I,J),I=IS,IE)
            ENDDO
		 ENDIF
         CLOSE (25)
	  ENDIF

      RETURN
	  END
