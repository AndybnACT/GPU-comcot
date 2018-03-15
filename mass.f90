!----------------------------------------------------------------------
      SUBROUTINE MASS (LO)
! ....SOLVE CONTINUITY EQUATION 
! OPTION:
!	0 - LINEAR SWE WITHOUT DISPERSION ADJUSTMENT
!	1 - NONLINEAR SWE WITHOUT DISPERSION ADJUSTMENT
!	2 - LINEAR SWE WITH DISPERSION ADJUSTMENT
!	3 - NONLINEAR SWE WITH DISPERSION ADJUSTMENT
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: LO

!.....FOR SPHERICAL COORDINATES
      IF (LO%LAYCORD .EQ. 0) THEN
		 SELECT CASE (LO%LAYGOV)
			CASE (0)
				CALL MASS_S (LO)
			CASE (1)
				CALL CONMASS_S (LO)
			CASE (2)
				CALL MASS_S_D (LO)
			CASE (3)
				CALL MASS_S_D (LO)
			CASE (9)
			    CALL CONMASS_S (LO)
			CASE DEFAULT
				CALL MASS_S (LO)
		 END SELECT
!.....FOR CARTESIAN COORDINATES
	  ELSE
		 SELECT CASE (LO%LAYGOV)
			CASE (0)
				CALL MASS_C (LO)
			CASE (1) 
				CALL CONMASS (LO)
			CASE (2)
				CALL MASS_C_D (LO)
			CASE (3)
				CALL MASS_C_D (LO)
			CASE (5)
				CALL CONMASS (LO)
			CASE (9)
			    CALL CONMASS (LO)
			CASE DEFAULT
				CALL MASS_C (LO)
		 END SELECT
	  ENDIF

	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE MASS_S (L)
! ....SOLVE CONTINUITY EQUATION UNDER SPHERICAL COORD. 
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/
!
      IS = 2
	  JS = 2
	  IF (L%ID .EQ. 1) THEN
	     IE = L%NX-1    !FOR OUTEST LAYER
	     JE = L%NY-1
	  ELSE
	     IE = L%NX      !FOR INNER LAYER
         JE = L%NY
      END IF
      DO J=JS,JE
         JM1 = J-1
         DO I=IS,IE
            IF (L%H(I,J) .GT. GX) THEN
               ZZ = L%Z(I,J,1)-L%R1(I,J)*(L%M(I,J,1)-L%M(I-1,J,1))	&
								-L%R11(I,J)*(L%N(I,J,1)*L%R6(J)		&
								-L%N(I,JM1,1)*L%R6(JM1))
			   IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4)		&
								ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))
			   IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
!*			   IF (ABS(ZZ) .GT. UB) ZZ=SIGN(UB,ZZ)
			   !DEPRESSION CANNOT BE LESS THAN BOTTOM ELEVATION
			   IF ( (ZZ+L%H(I,J)) .LE. EPS ) ZZ = -L%H(I,J)
               L%Z(I,J,2) = ZZ
		    ELSE
		       L%Z(I,J,2) = 0.0
            ENDIF
         END DO
      END DO

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE MASS_C (L)		 
! ....SOLVE CONTINUITY EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!
!     NOTES: ADD SUPPORT FOR DX\=DY (FOR HIGH-LATITUDE, 05/04/2007)
!            RX = L%R
!            RY = L%DT/L%DY
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/

	  RX = L%RX
	  RY = L%RY
!
      IS = 2
      JS = 2
	  IF (L%ID .EQ. 1) THEN  !OUTTEST LAYER				 
	     IE = L%NX -1
	     JE = L%NY -1
	  ELSE				  ! INNER LAYER
	     IE = L%NX
	     JE = L%NY
	  ENDIF

	  IF (L%DIM.EQ.1) THEN
		 L%N(:,:,:) = 0.0
		 JS = NINT(L%NY/2.0)
		 JE = NINT(L%NY/2.0)
	  ENDIF

      DO J=JS,JE
         DO I=IS,IE
            IF (L%H(I,J) .GT. GX) THEN
               ZZ = L%Z(I,J,1) - RX*(L%M(I,J,1)-L%M(I-1,J,1))			&
							- RY*(L%N(I,J,1)-L%N(I,J-1,1))
			   IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4)			&
									ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))  
			   IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
			   !DEPRESSION CANNOT BE LESS THAN BOTTOM ELEVATION
!*			   IF (ABS(ZZ) .GT. UB) ZZ=SIGN(UB,ZZ)  
			   IF ( (ZZ+L%H(I,J)) .LE. EPS ) ZZ = -L%H(I,J)
               L%Z(I,J,2) = ZZ
		    ELSE
               L%Z(I,J,2) = 0.0
            ENDIF
		    IF (L%DIM.EQ.1) L%Z(I,:,2) = L%Z(I,J,2)
         END DO
      END DO

      RETURN
      END



!----------------------------------------------------------------------
      SUBROUTINE CONMASS (L)
!.....SOVLE CONTINUITY EQUATION (NONLINEAR) IN CARTESIAN COORD.
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!.....DD  TOTAL WATER DEPTH
!----------------------------------------------------------------------
      USE LAYER_PARAMS
!	  USE DFLIB
      TYPE (LAYER) 	:: L
	  REAL FM
!	  INTEGER ERR, ERR_TABLE(L%NX,L%NY)
!	  INTEGER(2) STATUS
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      RX = L%RX
	  RY = L%RY

	  IF ( L%ID .EQ. 1 ) THEN  !OUTTEST LAYER !20
		 IS = 2
		 JS = 2
	     IE = L%NX -1
		 JE = L%NY -1
	  ELSE
		 IS = 2
		 JS = 2
	     IE = L%NX
		 JE = L%NY
	  ENDIF

	  IF (L%DIM .EQ. 1) THEN
		 L%N(:,:,:) = 0.0
		 JS = NINT(L%NY/2.0)
		 JE = NINT(L%NY/2.0)
	  ENDIF

!*	  L%DZ = 0.0
      DO J = JS, JE		
         DO I = IS, IE	
            IF (L%H(I,J) .GT. ELMAX) THEN
               ZZZ = L%Z(I,J,1) - RX*(L%M(I,J,1)-L%M(I-1,J,1))		&
			                 - RY*(L%N(I,J,1)-L%N(I,J-1,1))
			   IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4)		&
							ZZZ = ZZZ-(L%HT(I,J,2)-L%HT(I,J,1))   
               IF (ABS(ZZZ) .LT. EPS) ZZZ = 0.0
               DD = ZZZ + L%H(I,J)
               IF (DD .GE. GX) THEN
                  L%DZ(I,J,2) = DD
                  L%Z(I,J,2)  = ZZZ
               ELSE
                  L%DZ(I,J,2) = 0.0
                  L%Z(I,J,2)  = -L%H(I,J)
               ENDIF
            ELSE
               L%DZ(I,J,2) = 0.0
               L%Z(I,J,2)  = -L%H(I,J)
            ENDIF
		    IF (L%DIM.EQ.1) L%Z(I,:,2) = L%Z(I,JS,2)
         ENDDO
	  ENDDO

      RETURN
      END
      
!----------------------------------------------------------------------
      SUBROUTINE CONMASS_S (L)
!.....SOVLE CONTINUITY EQUATION (NONLINEAR) IN SPHERICAL COORD.
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!.....DD  TOTAL WATER DEPTH
!----------------------------------------------------------------------
      USE LAYER_PARAMS
!	  USE DFLIB
      TYPE (LAYER) 	:: L
	  REAL FM
!	  INTEGER ERR, ERR_TABLE(L%NX,L%NY)
!	  INTEGER(2) STATUS
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      IS = 2
	  JS = 2
	  IF (L%ID .EQ. 1) THEN
	     IE = L%NX-1    !FOR OUTEST LAYER
	     JE = L%NY-1
	  ELSE
	     IE = L%NX      !FOR INNER LAYER
         JE = L%NY
      END IF
      DO J=JS,JE
         JM1 = J-1
         DO I=IS,IE
            IF (L%H(I,J) .GE. ELMAX) THEN
               ZZ = L%Z(I,J,1)-L%R1(I,J)*(L%M(I,J,1)-L%M(I-1,J,1))	&
			        - L%R11(I,J)*(L%N(I,J,1)*L%R6(J)-L%N(I,JM1,1)	&
					* L%R6(JM1))
               IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4)		&
									ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))   
			   IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
               DD = ZZ + L%H(I,J)
               IF (DD .GE. GX) THEN
                  L%DZ(I,J,2) = DD
                  L%Z(I,J,2)  = ZZ
               ELSE
                  L%DZ(I,J,2) = 0.0
                  L%Z(I,J,2)  = - L%H(I,J)
               ENDIF
            ELSE
               L%DZ(I,J,2) = 0.0
               L%Z(I,J,2)  = - L%H(I,J)
            ENDIF
         END DO
      END DO

      RETURN
      END

