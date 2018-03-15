!----------------------------------------------------------------------
      SUBROUTINE MOMENT (LO)
! ....SOLVE MOMENTUM EQUATION
!OPTIONS:
!	0 - LINEAR EQUATION WITHOUT DISPERSION ADJUSTMENT
!	1 - NONLINEAR EQUATION WITHOUT DISPERSION ADJUSTMENT
!	2 - LINEAR EQUATION WITH DISPERSION ADJUSTMENT
!	3 - NONLINEAR EQUATION WITH DISPERSION ADJUSTMENT
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: LO
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!.....FOR SPHERICAL COORDINATES
      IF (LO%LAYCORD .EQ. 0) THEN
		 SELECT CASE (LO%LAYGOV)
			CASE (0)
				CALL MOMT_S (LO)
			CASE (1)
				CALL CONMOME_S (LO)
			CASE (2)
				CALL MOMT_S_D (LO)
			CASE (3)
				CALL MOMT_S_D (LO)
			CASE (9)
			    CALL CONMOME_S (LO)
			CASE DEFAULT
				CALL MOMT_S (LO)
		 END SELECT
!.....FOR CARTESIAN COORDINATES
	  ELSE
		 SELECT CASE (LO%LAYGOV)
			CASE (0)
				CALL MOMT_C (LO)
			CASE (1)
				CALL CONMOME (LO)
			CASE (2)
				CALL MOMT_C_D (LO)
			CASE (3)
				CALL MOMT_C_D (LO)
			CASE (5)
				CALL CONMOME (LO)
			CASE (9)
			    CALL CONMOME (LO)
			CASE DEFAULT
				CALL MOMT_C (LO)
		 END SELECT
	  ENDIF

	  RETURN
	  END



!----------------------------------------------------------------------
      SUBROUTINE MOMT_S (L)
! ....SOLVE MOMENTUM EQUATION (LINEAR) IN SPHERICAL COORD.
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL SCM
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA TWLVTH/0.08333333333333/
	  SCM = 0.0
!
      IE = L%NX-1
      JE = L%NY-1
	  IS = 2
	  JS = 2
	  IF (L%ID .EQ. 1) THEN
	     IS = 1
		 JS = 1
      END IF
	  DO I=IS,IE
         IP1 = I+1
         DO J=JS,L%NY
			IF ((L%H(I,J).GT.GX) .AND. (L%H(IP1,J).GT.GX)) THEN
				JM1 = J-1
				JP1 = J+1
				IF (JM1.LE.1) JM1 = 1
				IF (JP1.GE.L%NY) JP1 = L%NY
				TOT_N = L%N(I,J,1) + L%N(IP1,J,1) + L%N(I,JM1,1)	&
												+ L%N(IP1,JM1,1)
				XM = L%M(I,J,1) - L%R2(I,J)*(L%Z(IP1,J,2)			&
								- L%Z(I,J,2)) + L%R3(I,J)*TOT_N
				! IF (L%MODSCM .EQ. 0) THEN
				! 	SCM =  L%R2(I,J)*TWLVTH*((L%Z(IP1,JP1,2)		&
				! 				- 2*L%Z(IP1,J,2)+L%Z(IP1,JM1,2))	&
				! 				- (L%Z(I,JP1,2)-2*L%Z(I,J,2)		&
				! 				+ L%Z(I,JM1,2)))
				! 	XM = XM - SCM
				! ENDIF
				IF (ABS(XM) .LT. EPS) XM = ZERO
				L%M(I,J,2) = XM
			ELSE
			    L%M(I,J,2) = 0.0
			END IF
            ! if ( i .EQ. 448 .AND. j .EQ. 783 ) then
            !     WRITE(*,*) L%R2(I,J), L%R3(I,J),  L%N(I,J,1),  L%N(IP1,J,1), L%N(I,JM1,1), L%N(IP1,JM1,1)
            !     WRITE(*,*) L%M(I,J,1), L%Z(IP1,J,2), L%Z(I,J,2)
            ! end if
         END DO
      END DO
!
      DO J=JS,JE
         JP1 = J+1
         DO I=IS,L%NX
			IF ((L%H(I,J).GT.GX) .AND. (L%H(I,JP1).GT.GX)) THEN
				IM1 = I-1
				IP1 = I+1
				IF (IM1.LE.1) IM1 = 1
				IF (IP1.GE.L%NX) IP1 = L%NX
				TOT_M = L%M(IM1,J,1) + L%M(IM1,JP1,1) + L%M(I,J,1)	&
													+ L%M(I,JP1,1)
				XN = L%N(I,J,1) - L%R4(I,J)*(L%Z(I,JP1,2)			&
								- L%Z(I,J,2))-L%R5(I,J)*TOT_M
				! IF (L%MODSCM .EQ. 0) THEN
				! 	SCM = L%R4(I,J)*TWLVTH*((L%Z(IP1,JP1,2)			&
				! 				- 2*L%Z(I,JP1,2)+L%Z(IM1,JP1,2))	&
			    !                 - (L%Z(IP1,J,2)-2*L%Z(I,J,2)		&
				! 				+ L%Z(IM1,J,2)))
				! 	XN = XN - SCM
				! ENDIF
				IF (ABS(XN) .LT. EPS) XN = ZERO
				L%N(I,J,2) = XN
			ELSE
			    L%N(I,J,2) = 0.0
			END IF
            ! if ( i .EQ. 450 .AND. j .EQ. 783 ) then
            !     WRITE(*,*)L%R4(I,J),L%R5(I,J),L%M(I,J,1),L%M(I,JP1,1),L%M(IM1,JP1,1),L%M(IM1,J,1)
            !     WRITE(*,*)L%N(I,J,1),L%Z(I,JP1,2),L%Z(I,J,2)
            ! end if
         END DO
      END DO
!
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE MOMT_C (L)
!.....SOLVE MOMENTUM EQUATION (LINEAR) IN CARTESIAN COORD.
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL SCM, GRX,GRY,FF,FM
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA TWLVTH/0.08333333333333/

	  SCM = 0.0
	  CONST = 0.5*L%DT*GRAV
	  FM = L%FRIC_COEF
	  FF = 0.0

	  RX = L%RX
	  RY = L%RY
	  GRX = L%GRX
	  GRY = L%GRY
!
      IE = L%NX-1
      JE = L%NY-1
	  IF (L%ID .EQ. 1) THEN
	     IS = 1
		 JS = 1
	  ELSE
	     IS = 2
		 JS = 2
	  END IF

	  IF (L%DIM.EQ.1) THEN
		 L%N(:,:,:) = 0.0
		 JS = NINT(L%NY/2.0)
		 JE = NINT(L%NY/2.0)
	  ENDIF

      DO I=IS,IE
         IP1 = I+1
         DO J=JS,JE+1
			IF ((L%H(I,J).GT.GX) .AND. (L%H(IP1,J).GT.GX)) THEN
				JM1 = J-1
				JP1 = J+1
				IF (JM1 .LT. 1) JM1 = 1
				IF (JP1 .GT. L%NY ) JP1 = L%NY
				HM = L%HP(I,J)+0.5*(L%Z(I,J,2)+L%Z(IP1,J,2))
				XM = L%M(I,J,1)-L%GRX*HM*(L%Z(IP1,J,2)-L%Z(I,J,2))
!...........USE BOTTOM FRICTION
				IF (L%FRIC_SWITCH.NE.1) THEN
					IF (L%FRIC_SWITCH .EQ. 2) FM = L%FRIC_VCOEF(I,J)
					XQQ = 0.25*(L%N(I,J,1) + L%N(IP1,J,1)			&
								+ L%N(I,JM1,1) + L%N(IP1,JM1,1))
!					DF = 0.5*(L%H(I,J)+L%H(IP1,J))
					IF (HM.GE.1.0E-3) THEN
						FF = CONST*FM*FM*SQRT(L%M(I,J,1)**2			&
										+ XQQ**2)/HM**(2.3333333)
					ELSE
						FF = 0.0
					ENDIF
					XM = XM - FF*L%M(I,J,1)
				ENDIF
				IF (L%MODSCM .EQ. 0) THEN
					SCM = L%GRX*TWLVTH*HM*((L%Z(IP1,JP1,2)			&
								- 2*L%Z(IP1,J,2)+L%Z(IP1,JM1,2))	&
								- (L%Z(I,JP1,2)-2*L%Z(I,J,2)		&
								+ L%Z(I,JM1,2)))
					XM = XM - SCM
				ENDIF
				IF (L%FRIC_SWITCH.NE.1) XM = XM/(1.+FF)
				IF (ABS(XM) .LT. EPS) XM = ZERO
				L%M(I,J,2) = XM
			ELSE
				XM = 0.0
			    L%M(I,J,2) = XM
			END IF
			IF (L%DIM.EQ.1) L%M(I,:,2) = XM
		 END DO
      END DO
!
	  IF (L%DIM .NE. 1) THEN
      DO J=JS,JE
         JP1 = J+1
         DO I=IS,IE+1
			IF ((L%H(I,J).GT.GX) .AND. (L%H(I,JP1).GT.GX)) THEN
				IM1 = I-1
				IP1 = I+1
				IF (IM1 .LT. 1) IM1 = 1
				IF (IP1 .GE. L%NX) IP1 = L%NX
				HN = L%HQ(I,J)+0.5*(L%Z(I,J,2)+L%Z(I,JP1,2))
				XN = L%N(I,J,1)-L%GRY*HN*(L%Z(I,JP1,2)-L%Z(I,J,2))
!...........USE BOTTOM FRICTION
				IF (L%FRIC_SWITCH.NE.1) THEN
					IF (L%FRIC_SWITCH .EQ. 2) FM = L%FRIC_VCOEF(I,J)
					XPP = 0.25*(L%M(I,J,1) + L%M(I,JP1,1)			&
								+ L%M(IM1,J,1) + L%M(IM1,JP1,1))
!					DF = 0.5*(L%H(I,J)+L%H(I,JP1))
					IF (HN.GE.1.0E-5) THEN
						FF = CONST*FM*FM*SQRT(L%N(I,J,1)**2			&
										+ XPP**2)/HN**(2.3333333)
					ELSE
						FF = 0.0
					ENDIF
					XN = XN - FF*L%N(I,J,1)
				ENDIF

				IF (L%MODSCM .EQ. 0) THEN
					SCM = L%GRY*TWLVTH*HN*((L%Z(IP1,JP1,2)			&
								- 2.*L%Z(I,JP1,2)+L%Z(IM1,JP1,2))	&
			                    - (L%Z(IP1,J,2)-2.*L%Z(I,J,2)		&
								+ L%Z(IM1,J,2)))
					XN = XN - SCM
				ENDIF

				IF (L%FRIC_SWITCH.NE.1) XN = XN/(1.0+FF)
				IF (ABS(XN) .LT. EPS) XN = ZERO
				L%N(I,J,2) = XN
			ELSE
			    L%N(I,J,2) = 0.0
			END IF
         END DO
      END DO
	  ENDIF
!
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE CONMOME (L)
!....SOVLE NONLINEAR MOMENTUM EQUATION, CARTESIAN COORD.
!....DP,DQ: TOTAL WATER DEPTH AT DISCHARGE POINT
!....HP,HQ: STILL WATER DEPTH AT DISCHARGE POINT
!....HZ   : WATER DEPTH (+: WATER; -: LAND)
!....Z    : FREE SURFACE ELEVATION
!    P,Q  : VOLUME FLUX IN X AND Y DIRECTION;
!    NX,NY: DIMENSION OF DOMAIN IN X AND Y DIRECTION
!    LAST REVISE: NOV.17 2008 BY XIAOMING WANG
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
      REAL DP(L%NX,L%NY,2),DQ(L%NX,L%NY,2),DZ(L%NX,L%NY,2)
	  REAL HP(L%NX,L%NY),HQ(L%NX,L%NY)
      REAL HZ(L%NX,L%NY),P(L%NX,L%NY,2),Q(L%NX,L%NY,2),Z(L%NX,L%NY,2)
!	  REAL T_BREAK(L%NX,L%NY)
	  REAL RX, RY, DT, GRX, GRY,SCM,NU
	  INTEGER IX, JY, IFRIC, LIN_CHK, MOD_SCM
	  INTEGER LAYID,SPAN,SPANX,SPANY
	  INTEGER NOFLUX   !FLAG TO DETERMINE IF FLUX IS CALCULATED
!	  INTEGER MASK(L%NX,L%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  SCM = 0.0
	  ADVX = 0.0
	  ADVY = 0.0
	  NOFLUX = 0
!	  L%MASK = 0

	  IDIM = L%DIM

	  IF (IDIM.EQ.1) L%N(:,:,:) = 0.0
	  IF (L%LAYGOV.EQ.9) L%LINCHK = 0

	  LAYID=L%ID
	  IX = L%NX
	  JY = L%NY
      HZ = L%H
	  P = L%M
	  Q = L%N
	  Z = L%Z
	  DZ = L%DZ
	  DT = L%DT
	  RX = L%RX
	  RY = L%RY
      GRX = L%GRX
	  GRY = L%GRY

	  IFRIC = L%FRIC_SWITCH
	  FM = L%FRIC_COEF
	  LIN_CHK = L%LINCHK
	  MOD_SCM = L%MODSCM

	  HP = L%HP
	  HQ = L%HQ

!.....FOR FRICTIONAL EFFECTS
      IF (IFRIC .NE. 1) THEN
         CONST = 0.5*DT*GRAV
      ELSE
         CONST = 0.0
      ENDIF

!CALCULATE TOTAL WATER DEPTH AND TOTAL WATER DEPTH AT DISCHARGE POINT
      DO I=1,IX  !-1
	  	 IP1 = I+1   !!
		 IF (IP1 .GE. IX) IP1=IX  !!
         DO J=1,JY  !-1
			JP1 = J+1  !!
			IF (JP1 .GE. JY) JP1=JY   !!
			DP1 = 0.25*(DZ(I,J,2)+DZ(I,J,1)+DZ(IP1,J,2)+DZ(IP1,J,1))
			DP2 = 0.50*(DZ(I,J,2)+DZ(IP1,J,2))
			DQ1 = 0.25*(DZ(I,J,2)+DZ(I,J,1)+DZ(I,JP1,2)+DZ(I,JP1,1))
			DQ2 = 0.50*(DZ(I,J,2)+DZ(I,JP1,2))

			IF (DP1 .LT. GX) DP1 = 0.0
			IF (DP2 .LT. GX) DP2 = 0.0
			IF (DQ1 .LT. GX) DQ1 = 0.0
			IF (DQ2 .LT. GX) DQ2 = 0.0

			DP(I,J,1) = DP1
			DP(I,J,2) = DP2
			DQ(I,J,1) = DQ1
			DQ(I,J,2) = DQ2
         END DO
      END DO
	  !WIDTH OF BUFFER GRIDS BEFORE IMPLEMENTING NONLINEAR CALCULATION
	  SPAN = 10
	  ISPANS = SPAN  !
	  ISPANE = IX - SPAN
	  JSPANS = SPAN
	  JSPANE = JY - SPAN

!.....SOLVE MOMENTUM EQUATION IN X DIRECTION
!!3  MODIFY THE DOMAIN OF COMPUTATION

      IS = 2
	  JS = 2
	  IE = IX-1
!	  IE = IX
	  JE = JY
	  IF (LAYID .EQ. 1) THEN
	      IS = 1
		  JS = 1
		  IE = IX-1
		  JE = JY-1
	  ENDIF

	  IF (IDIM.EQ.1) JS = FLOOR(L%NY/2.0)
	  IF (IDIM.EQ.1) JE = FLOOR(L%NY/2.0)


      DO I = IS, IE
	     IP2 = I+2
		 IP1 = I+1
		 IM1 = I-1
		 IF (IM1.LE.1) IM1 = 1
		 IF (IP1 .GE. IX) IP1 = IX
		 IF (IP2 .GE. IX) IP2 = IX
         DO J = JS, JE
			SCM = 0.0
			NOFLUX = 0  ! 0: CALC. FLUXES; 1: DON'T CALC. FLUXES
!			L%MASK(I,J) = 0
			JP1 = J+1
			JM1 = J-1
			IF (JM1.LE.1) JM1 = 1
			IF (JP1 .GE. JY) JP1 = JY
!..CALCULATE X-DIRECTION LINEAR TERMS
	        IF (HZ(I,J).LE.ELMAX .OR. HP(I,J) .LE. ELMAX) THEN
				P(I,J,2) = 0.0
				NOFLUX = 1
            ELSE
!DETERMINE MOVING BOUNDARY (CHANGE 0.0 TO GX)
				IF (DZ(I,J,2).GT.GX .AND. DZ(IP1,J,2).GT.GX) THEN
					DD = DP(I,J,2)
					DF = DP(I,J,1)
				ELSEIF (DZ(I,J,2).GT.GX .AND. DZ(IP1,J,2).LE.GX		&
						.AND. HZ(IP1,J)+Z(I,J,2).GT.GX) THEN
!*					DD = HZ(IP1,J) + Z(I,J,2)
!					DD = HP(I,J) + 0.5*(Z(I,J,2)+Z(IP1,J,2))
					DD = 0.5*DZ(I,J,2)
					DF = DD
!					L%MASK(I,J) = 1
				ELSEIF (DZ(I,J,2).LE.GX .AND. DZ(IP1,J,2).GT.GX		&
						.AND. HZ(I,J)+Z(IP1,J,2) .GT. GX) THEN
!*					DD = HZ(I,J) + Z(IP1,J,2)
!					DD = HP(I,J) + 0.5*(Z(I,J,2)+Z(IP1,J,2))
					DD = 0.5*DZ(IP1,J,2)
					DF = DD
!					L%MASK(I,J) = 1
				ELSE
					P(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF

				IF (DD .LT. GX) THEN
					P(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF
			ENDIF


!---------->COMPUTE X FLUX COMPONENT
            IF (NOFLUX .NE. 1) THEN  !--<<<<X<<<<
!..ESTABLISH A LOWER BOUND
			IF (DF .LT. 1.0E-3) THEN
				DF = 1.0E-3
			ENDIF

			XQQ = 0.25*(Q(I,J,1) + Q(IP1,J,1)						&
								+ Q(I,JM1,1) + Q(IP1,JM1,1)) !7!

!..FIXED FRICTIONAL COEFFICIENT - MANNING
			!USING VARIABLE MANNING COEF.
			IF (IFRIC .EQ. 2) FM = L%FRIC_VCOEF(I,J)
			IF (IFRIC .NE. 1) THEN
!*				XQQ = 0.25*(Q(I,J,1) + Q(IP1,J,1)					&
!*								+ Q(I,JM1,1) + Q(IP1,JM1,1)) !7!
				FF = CONST*FM*FM*SQRT(P(I,J,1)**2					&
								+ XQQ**2)/DF**(2.3333333)
			ELSE
				FF = 0.0
			ENDIF

!.....>>COMPUTE LINEAR PART IN NONLINEAR MOMENTUM EQUATION<<....
			XP = P(I,J,1) - GRX*DD*(Z(IP1,J,2)-Z(I,J,2))
!!7!->  MODIFIED OR OLD SCHEME FOR MOMENTUM EQUATION
			IF (MOD_SCM.EQ.0 .AND. J.LT.JY) THEN
				SCM = GRX*DD/12.									&
						*((Z(IP1,JP1,2)-2*Z(IP1,J,2)+Z(IP1,JM1,2))	&
						-(Z(I,JP1,2)-2*Z(I,J,2)+Z(I,JM1,2)))
				XP = XP - SCM
			ENDIF

			IF (IFRIC .NE. 1) XP = XP - FF*P(I,J,1)
!!7! <-

!.....>>COMPUTE CONVECTION TERMS IN NONLINEAR MOMENTUM EQUATION<<...
			IF (LIN_CHK.EQ.1 .AND. DP(I,J,1).GE.GX .AND.			&
				(I.GT.ISPANS .AND. I.LT.ISPANE .AND.J.GT.JSPANS		&
				.AND. J.LT.JSPANE)) THEN
!			IF (LIN_CHK.EQ.1 .AND. MASK(I,J).EQ.1) THEN

!*				XQQ = 0.25*(Q(I,J,1) + Q(IP1,J,1)					&
!*								+ Q(I,JM1,1) + Q(IP1,JM1,1)) !7!
				ADVX = 0.0
				ADVY = 0.0
!..UPWIND SCHEME FOR X-DIRECTION VOLUME FLUX COMPONENT
				IF (P(I,J,1) .LT. 0.0) THEN   !CHANGE .LE. TO .LT.
					IF (DP(IP1,J,1) .LT. GX .OR.					&
										DZ(IP1,J,2) .LT. GX) THEN
						ADVX = RX*(-P(I,J,1)**2/DP(I,J,1))
					ELSE
						ADVX = RX*(P(IP1,J,1)**2/DP(IP1,J,1)		&
										-P(I,J,1)**2/DP(I,J,1))
					ENDIF
!*					IF (DP(IP1,J,1).GE.GX .AND. DP(I,J,1).GE.GX) THEN
!*						XP = XP - RX*(P(IP1,J,1)**2/DP(IP1,J,1)		&
!*								- P(I,J,1)**2/DP(I,J,1))
!*					ENDIF
				ELSE
					IF (DP(IM1,J,1) .LT. GX .OR.					&
											DZ(I,J,2) .LT. GX) THEN
						ADVX = RX*(P(I,J,1)**2/DP(I,J,1))
					ELSE
						ADVX = RX*(P(I,J,1)**2/DP(I,J,1)			&
										-P(IM1,J,1)**2/DP(IM1,J,1))
					ENDIF
!*					IF (DP(I,J,1).GE.GX .AND. DP(IM1,J,1).GT.GX) THEN
!*						ADVX = RX*(P(I,J,1)**2/DP(I,J,1)			&
!*								- P(IM1,J,1)**2/DP(IM1,J,1))
!*					ENDIF

				ENDIF

!..UPWIND SCHEME FOR Y-DIRECTION VOLUME FLUX COMPONENT
				IF (XQQ .LT. 0.0) THEN
!*					XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)			&
!*									+ Q(I,J,1) + Q(IP1,J,1))
					IF (DZ(I,JP1,2) .LT. GX .OR.					&
										DZ(IP1,JP1,2) .LT. GX) THEN
						ADVY = RY*(-P(I,J,1)*XQQ/DP(I,J,1))
					ELSEIF (DP(I,JP1,1) .LT. GX) THEN
						ADVY = RY*(-P(I,J,1)*XQQ/DP(I,J,1))
					ELSE
						XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)		&
									+ Q(I,J,1) + Q(IP1,J,1))
						ADVY = RY*(P(I,JP1,1)*XQE/DP(I,JP1,1)		&
										-P(I,J,1)*XQQ/DP(I,J,1))
					ENDIF
!*					IF (DP(I,JP1,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)		&
!*									+ Q(I,J,1) + Q(IP1,J,1))
!*						ADVY = RY*(P(I,JP1,1)*XQE/DP(I,JP1,1)		&
!*								- P(I,J,1)*XQQ/DP(I,J,1))
!*					ENDIF
				ELSE
!*					XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)			&
!*									+ Q(I,J-2,1) + Q(IP1,J-2,1))

					IF (DZ(I,JM1,2) .LT. GX .OR.					&
										DZ(IP1,JM1,2) .LT. GX) THEN
						ADVY = RY*(P(I,J,1)*XQQ/DP(I,J,1))
					ELSEIF (DP(I,JM1,1) .LT. GX) THEN
						ADVY = RY*(P(I,J,1)*XQQ/DP(I,J,1))
					ELSE
						XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)		&
									+ Q(I,J-2,1) + Q(IP1,J-2,1))
						ADVY = RY*(P(I,J,1)*XQQ/DP(I,J,1)			&
									-P(I,JM1,1)*XQE/DP(I,JM1,1))
					ENDIF
!*					IF (DP(I,JM1,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)		&
!*									+ Q(I,J-2,1) + Q(IP1,J-2,1))
!*						ADVY = RY*(P(I,J,1)*XQQ/DP(I,J,1)			&
!*								- P(I,JM1,1)*XQE/DP(I,JM1,1))
!*					ENDIF
				ENDIF
				CALL WAVE_BREAKING (L,XP,XQ,ADVX,ADVY,I,J,1)
				XP = XP - ADVX - ADVY
			ENDIF
		    IF (IFRIC.NE.1) XP = XP/(1.0+FF)

!........LIMIT THE DISCHARGE
			IF (ABS(XP) .LT. EPS) XP = 0.0
			PQ_LIMIT = V_LIMIT*DD
			IF (XP .GT. PQ_LIMIT) XP = PQ_LIMIT
			IF (XP .LT. -PQ_LIMIT) XP =-PQ_LIMIT

			P(I,J,2) = XP
			IF (IDIM.EQ.1) P(I,:,2) = XP

			ELSE   ! IF NOFLUX=1
			   P(I,J,2) = 0.0
			ENDIF  !--<<<<X<<<<
	     ENDDO
	  ENDDO
       !7!

!.................................................................
!.....SOLVE MOMENTUM EQUATION IN Y DIRECTION
      IF (IDIM.NE.1) THEN
      IS = 2
	  JS = 2
	  IE = IX
	  JE = JY-1
!	  JE = JY
	  IF (LAYID .EQ. 1) THEN
	      IS = 1
		  JS = 1
		  IE = IX-1
		  JE = JY-1
	  ENDIF
      DO I = IS, IE          !3
	  	 IP1 = I+1
		 IM1 = I-1
		 IF (IM1.LE.1) IM1 = 1
		 IF (IP1 .GE. IX) IP1 = IX
		 DO J = JS, JE        !3
			SCM = 0.0
			NOFLUX = 0    ! FLAG: 0-CALCULATE Y FLUX; 1-DON'T CALCULATE Y FLUX
!			L%MASK(I,J) = 0
			JP2 = J+2
			JP1 = J+1
			JM1 = J-1
			IF (JM1.LE.1) JM1 = 1
			IF (JP1 .GE. JY) JP1 = JY
			IF (JP2 .GE. JY) JP2 = JY
!..GUESS THE RUNUP NEVER REACH AT HZ(I,J)=ELMAX
			IF (HZ(I,J).LE.ELMAX .OR. HQ(I,J).LE.ELMAX) THEN
				Q(I,J,2) = 0.0
				NOFLUX = 1
			ELSE
!DETERMINE MOVING BOUNDARY SCHEME
				IF (DZ(I,J,2).GT.GX .AND. DZ(I,JP1,2).GT.GX) THEN
					DD = DQ(I,J,2)
					DF = DQ(I,J,1)
				ELSEIF (DZ(I,J,2).GT.GX .AND. DZ(I,JP1,2).LE.GX		&
						.AND. HZ(I,JP1) + Z(I,J,2) .GT. GX) THEN
!*     				DD = HZ(I,JP1) + Z(I,J,2)
!					DD = HQ(I,J) + 0.5*(Z(I,J,2)+Z(I,JP1,2))
					DD = 0.5*DZ(I,J,2)
					DF = DD
				ELSEIF (DZ(I,J,2).LE.GX .AND. DZ(I,JP1,2).GT.GX		&
						.AND. HZ(I,J)+Z(I,JP1,2) .GT. GX) THEN
!*     				DD = HZ(I,J) + Z(I,JP1,2)
!					DD = HQ(I,J) + 0.5*(Z(I,J,2)+Z(I,JP1,2))
					DD = 0.5*DZ(I,JP1,2)
					DF = DD
				ELSE
	 				Q(I,J,2) = 0.0
           			NOFLUX = 1
				ENDIF

				IF (DD .LT. GX) THEN
					Q(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF
			ENDIF

!...........COMPUTE Y FLUX WHEN NOFLUX.NE.1
            IF (NOFLUX .NE. 1) THEN    !--<<<<Y<<<<

!..CALCULATE Y-DIRECTION LINEAR TERMS
!..ESTABLISH A LOWER BOUND
			IF (DF .LT. 1.0E-3) THEN
				DF = 1.0E-3
			ENDIF
			XPP = 0.25*(P(I,J,1) + P(I,JP1,1)						&
								+ P(IM1,J,1) + P(IM1,JP1,1))  !7!
!..FIXED FRICTIONAL COEFFICIENT - MANNING
			IF (IFRIC .EQ. 2) FM = L%FRIC_VCOEF(I,J) !VARIABLE MANNING'S COEF.
			IF (IFRIC .NE. 1) THEN
!*				XPP = 0.25*(P(I,J,1) + P(I,JP1,1)					&
!*								+ P(IM1,J,1) + P(IM1,JP1,1))  !7!
				FF = CONST*FM*FM*SQRT(Q(I,J,1)**2					&
										+XPP**2)/DF**2.3333333
			ELSE
				FF = 0.0
			ENDIF

!...........SOLVE LINEAR VERSION OF MOMENTUM EQUATION IN Y DIRECTION
			XQ = Q(I,J,1) - GRY*DD*(Z(I,JP1,2)-Z(I,J,2))
!!7!->  MODIFIED OR OLD SCHEME FOR MOMENTUM EQUATION
			IF (MOD_SCM .EQ. 0 .AND. I.LT.IX) THEN
				SCM = GRY  * DD / 12.								&
					*((Z(IP1,JP1,2)-2*Z(I,JP1,2)+Z(IM1,JP1,2))		&
					-(Z(IP1,J,2)-2*Z(I,J,2)+Z(IM1,J,2)))
				XQ = XQ - SCM
			ENDIF

			IF (IFRIC .NE. 1) XP = XP - FF*Q(I,J,1)

			IF (LIN_CHK.EQ.1 .AND. DQ(I,J,1).GE.GX .AND.			&
				(I.GT.ISPANS .AND. I.LT.ISPANE .AND. J.GT.JSPANS	&
				.AND. J.LT.JSPANE)) THEN
!			IF (LIN_CHK.EQ.1 .AND. MASK(I,J).EQ.1) THEN

!*				XPP = 0.25*(P(I,J,1) + P(I,JP1,1)					&
!*								+ P(IM1,J,1) + P(IM1,JP1,1))  !7!
				ADVX = 0.0
				ADVY = 0.0
!..UPWIND SCHEME FOR Y-DIRECTION VOLUME FLUX COMPONENT
				IF (Q(I,J,1) .LT. 0.0) THEN
					IF (DQ(I,JP1,1) .LT. GX .OR.					&
										DZ(I,JP1,2) .LT. GX) THEN
						ADVY = RY*(-Q(I,J,1)**2/DQ(I,J,1))
					ELSE
						ADVY = RY*(Q(I,JP1,1)**2/DQ(I,JP1,1)		&
										-Q(I,J,1)**2/DQ(I,J,1))
					ENDIF
!*					IF (DQ(I,JP1,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						ADVY = RY*(Q(I,JP1,1)**2/DQ(I,JP1,1)		&
!*								- Q(I,J,1)**2/DQ(I,J,1))
!*					ENDIF
				ELSE
					IF (DQ(I,JM1,1) .LT. GX .OR.					&
											DZ(I,J,2) .LT. GX) THEN
						ADVY = RY*(Q(I,J,1)**2/DQ(I,J,1))
					ELSE
						ADVY = RY*(Q(I,J,1)**2/DQ(I,J,1)			&
									-Q(I,JM1,1)**2/DQ(I,JM1,1))
					ENDIF
!*					IF (DQ(I,JM1,1).GT.GX) THEN
!*					XQ = XQ - RY*(Q(I,J,1)**2/DQ(I,J,1)				&
!*							- Q(I,JM1,1)**2/DQ(I,JM1,1))
!*					ENDIF
				ENDIF

!..UPWIND SCHEME FOR X-DIRECTION VOLUME FLUX COMPONENT
				IF (XPP .LT. 0.0) THEN
!*					XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)			&
!*									+ P(I,J,1) + P(I,JP1,1))

					IF (DZ(IP1,J,2) .LT. GX .OR.					&
										DZ(IP1,JP1,2) .LT. GX) THEN
						ADVX = RX*(-Q(I,J,1)*XPP/DQ(I,J,1))
					ELSEIF (DQ(IP1,J,1) .LT. GX) THEN
						ADVX = RX*(-Q(I,J,1)*XPP/DQ(I,J,1))
					ELSE
						XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)		&
									+ P(I,J,1) + P(I,JP1,1))
						ADVX = RX*(Q(IP1,J,1)*XPE/DQ(IP1,J,1)		&
										-Q(I,J,1)*XPP/DQ(I,J,1))
					ENDIF
!*					IF (DQ(IP1,J,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)		&
!*									+ P(I,J,1) + P(I,JP1,1))
!*						XQ = XQ - RX*(Q(IP1,J,1)*XPE/DQ(IP1,J,1)	&
!*								- Q(I,J,1)*XPP/DQ(I,J,1))
!*					ENDIF
				ELSE
!*					XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)			&
!*									+ P(I-2,J,1) + P(I-2,JP1,1))

					IF (DZ(IM1,J,2) .LT. GX .OR.					&
										DZ(IM1,JP1,2) .LT. GX) THEN
						ADVX = RX*(Q(I,J,1)*XPP/DQ(I,J,1))
					ELSEIF (DQ(I-1,J,1) .LT. GX) THEN
						ADVX = RX*(Q(I,J,1)*XPP/DQ(I,J,1))
					ELSE
						XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)		&
									+ P(I-2,J,1) + P(I-2,JP1,1))
						ADVX = RX*(Q(I,J,1)*XPP/DQ(I,J,1)			&
									-Q(IM1,J,1)*XPE/DQ(IM1,J,1))
					ENDIF
!*					IF (DQ(IM1,J,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)		&
!*									+ P(I-2,J,1) + P(I-2,JP1,1))
!*						ADVX = RX*(Q(I,J,1)*XPP/DQ(I,J,1)
!*								- Q(IM1,J,1)*XPE/DQ(IM1,J,1))
!*					ENDIF
				ENDIF
				CALL WAVE_BREAKING (L,XP,XQ,ADVX,ADVY,I,J,2)
				XQ = XQ - ADVX - ADVY
			ENDIF
		    IF (IFRIC .NE. 1) XQ = XQ/(1.0+FF)
!...........LIMIT THE DISCHARGE
			IF (ABS(XQ) .LT. EPS) XQ = 0.0
			PQ_LIMIT = V_LIMIT*DD
	        IF (XQ .GT. PQ_LIMIT) XQ = PQ_LIMIT
			IF (XQ .LT. -PQ_LIMIT) XQ =-PQ_LIMIT

			Q(I,J,2) = XQ

			ELSE  ! IF NOFLUX = 1
			   Q(I,J,2) = 0.0
            ENDIF   !--<<<<Y<<<<
		  ENDDO
      ENDDO !7!

	  ENDIF

      L%M = P
	  L%N = Q

!*	  CALL WAVE_BREAKING (L)

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE CONMOME_S (L)
!....SOVLE NONLINEAR MOMENTUM EQUATION, SPHERICAL COORD.
!....DP,DQ: TOTAL WATER DEPTH AT DISCHARGE POINT
!....HP,HQ: STILL WATER DEPTH AT DISCHARGE POINT
!....HZ   : WATER DEPTH (+: WATER; -: LAND)
!....Z    : FREE SURFACE ELEVATION
!    P,Q  : VOLUME FLUX IN X AND Y DIRECTION;
!    IX,JY: DIMENSION OF DOMAIN IN X AND Y DIRECTION
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
      REAL DP(L%NX,L%NY,2),DQ(L%NX,L%NY,2),DZ(L%NX,L%NY,2)
	  REAL HP(L%NX,L%NY),HQ(L%NX,L%NY)
      REAL HZ(L%NX,L%NY),P(L%NX,L%NY,2),Q(L%NX,L%NY,2),Z(L%NX,L%NY,2)
	  REAL RX, RY, DT, GRX, GRY
	  REAL ADVX,ADVY
	  INTEGER IX, JY, IFRIC, LIN_CHK, MOD_SCM
	  INTEGER LAYID,SPANX,SPANY
	  INTEGER NOFLUX   !FLAG TO DETERMINE IF FLUX IS CALCULATED
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  NOFLUX = 0

	  IDIM = 2
	  IF (IDIM.EQ.1) L%N(:,:,:) = 0.0

	  IF (L%LAYGOV.EQ.9) L%LINCHK = 0

	  LAYID=L%ID
	  IX = L%NX
	  JY = L%NY
      HZ = L%H
	  P = L%M
	  Q = L%N
	  Z = L%Z
	  DZ = L%DZ
	  DT = L%DT
	  RX = L%RX
	  RY = L%RY
      GRX = L%GRX
	  GRY = L%GRY

	  IFRIC = L%FRIC_SWITCH
	  FM = L%FRIC_COEF
	  LIN_CHK = L%LINCHK
	  MOD_SCM = L%MODSCM

	  HP = L%HP
	  HQ = L%HQ

	  !WIDTH OF BUFFER GRIDS BEFORE IMPLEMENTING NONLINEAR CALCULATION
	  NN = 10
	  ISPANS = NN
	  ISPANE = L%NX-NN
	  JSPANS = NN
	  JSPANE = L%NY-NN

!..FOR FRICTIONAL EFFECTS
      IF (IFRIC .NE. 1) THEN
         CONST = 0.5*DT*GRAV
      ELSE
         CONST = 0.0
      ENDIF

!..CALCULATE TOTAL WATER DEPTH AND TOTAL WATER DEPTH AT DISCHARGE POINT
      DO I=1,IX  !-1
         DO J=1,JY  !-1
			IP1 = I+1   !!
			IF (IP1 .GE. IX) IP1=IX  !!
			JP1 = J+1  !!
			IF (JP1 .GE. JY) JP1=JY   !!
			DP1 = 0.25*(DZ(I,J,2)+DZ(I,J,1)+DZ(IP1,J,2)+DZ(IP1,J,1))
			DP2 = 0.50*(DZ(I,J,2)+DZ(IP1,J,2))
			DQ1 = 0.25*(DZ(I,J,2)+DZ(I,J,1)+DZ(I,JP1,2)+DZ(I,JP1,1))
			DQ2 = 0.50*(DZ(I,J,2)+DZ(I,JP1,2))

			IF (DP1 .LT. GX) DP1 = 0.0
			IF (DP2 .LT. GX) DP2 = 0.0
			IF (DQ1 .LT. GX) DQ1 = 0.0
			IF (DQ2 .LT. GX) DQ2 = 0.0

			DP(I,J,1) = DP1
			DP(I,J,2) = DP2
			DQ(I,J,1) = DQ1
			DQ(I,J,2) = DQ2
         END DO
      END DO
!.....SOLVE MOMENTUM EQUATION IN X DIRECTION
!!3  MODIFY THE DOMAIN OF COMPUTATION
      IS = 2
	  JS = 2
	  IE = IX-1
!	  IE = IX
	  JE = JY
	  IF (LAYID .EQ. 1) THEN
	      IS = 1
		  JS = 1
		  IE = IX-1
		  JE = JY-1
	  ENDIF

	  IF (IDIM.EQ.1) JS = FLOOR(L%NY/2.0)
	  IF (IDIM.EQ.1) JE = FLOOR(L%NY/2.0)


      DO I = IS, IE
		 IP1 = I+1
		 IM1 = I-1
		 IF (IM1.LE.1) IM1 = 1
		 IF (IP1 .GE. IX) IP1 = IX
         DO J = JS, JE
			NOFLUX = 0  !0: CALC FLUXES; 1: DON'T CALC FLUXES
			JP1 = J+1
			JM1 = J-1
			IF (JM1.LE.1) JM1 = 1
			IF (JP1 .GE. JY) JP1 = JY
!..CALCULATE X-DIRECTION LINEAR TERMS
	        IF (HZ(I,J).LE.ELMAX .OR. HP(I,J) .LE. ELMAX) THEN
				P(I,J,2) = 0.0
				NOFLUX = 1
            ELSE
!..MOVING BOUNDARY
				IF (DZ(I,J,2).GT.GX .AND. DZ(IP1,J,2).GT.GX) THEN
					DD = DP(I,J,2)
					DF = DP(I,J,1)
				ELSEIF (DZ(I,J,2).GT.GX .AND. DZ(IP1,J,2).LE.GX	&
						.AND. HZ(IP1,J)+Z(I,J,2).GT.GX) THEN
!*					DD = HZ(IP1,J) + Z(I,J,2)
					DD = HP(I,J) + 0.5*(Z(I,J,2)+Z(IP1,J,2))
					DF = DD
				ELSEIF (DZ(I,J,2).LE.GX .AND. DZ(IP1,J,2).GT.GX	&
						.AND. HZ(I,J)+Z(IP1,J,2) .GT. GX) THEN
!*					DD = HZ(I,J) + Z(IP1,J,2)
					DD = HP(I,J) + 0.5*(Z(I,J,2)+Z(IP1,J,2))
					DF = DD
				ELSE
					P(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF

				IF (DD .LT. GX) THEN
					P(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF
			ENDIF


!---------->COMPUTE X FLUX COMPONENT
            IF (NOFLUX .NE. 1) THEN  !--<<<<X<<<<
!..ESTABLISH A LOWER BOUND
			IF (DF .LT. 1.0E-3) THEN
				DF = 1.0E-3
			ENDIF

!..FIXED FRICTIONAL COEFFICIENT - MANNING
			!USING VARIABLE MANNING COEF.
			XQQ = 0.25*(Q(I,J,1)+Q(IP1,J,1)+Q(I,JM1,1)+Q(IP1,JM1,1))
			IF (IFRIC .EQ. 2) FM = L%FRIC_VCOEF(I,J)
			IF (IFRIC .NE. 1) THEN
!!!				XQQ = 0.25*(Q(I,J,1)+Q(IP1,J,1)						&
!!!									+Q(I,JM1,1)+Q(IP1,JM1,1)) !7!
				FF = CONST*FM*FM*SQRT(P(I,J,1)**2					&
										+XQQ**2)/DF**(2.3333333)
			ELSE
				FF = 0.0
			ENDIF

!.....>>COMPUTE LINEAR PART IN NONLINEAR MOMENTUM EQUATION<<....
!!7!->  MODIFIED OR OLD SCHEME FOR MOMENTUM EQUATION
			IF (MOD_SCM.EQ.0 .AND. J.LT.JY) THEN
				XP = P(I,J,1) - L%R2(I,J)*DD*(Z(IP1,J,2)-Z(I,J,2))	&
					- L%R2(I,J)  * DD / 12.							&
					* ((Z(IP1,JP1,2)-2*Z(IP1,J,2)+Z(IP1,JM1,2))		&
					-(Z(I,JP1,2)-2*Z(I,J,2)+Z(I,JM1,2)))
			ELSE
				XP = P(I,J,1) - L%R2(I,J)*DD*(Z(IP1,J,2)-Z(I,J,2))
			ENDIF
			IF (IFRIC.NE.1) XP = XP - FF*P(I,J,1)
!!!!!!!!COMPUTE CORIOLIS EFFECT
            XP = XP + 4.0*L%R3(I,J)*XQQ
!!7! <-

!.....>>COMPUTE CONVECTION TERMS IN NONLINEAR MOMENTUM EQUATION<<...
			IF (LIN_CHK.EQ.1 .AND. DP(I,J,1).GE.GX .AND.			&
				(I.GT.ISPANS .AND. I.LT.ISPANE .AND.				&
				J.GT.JSPANS .AND. J.LT.JSPANE) ) THEN

!!!				XQQ = 0.25*(Q(I,J,1)+Q(IP1,J,1)						&
!!!									+Q(I,JM1,1)+Q(IP1,JM1,1)) !7!
				ADVX = 0.0
				ADVY = 0.0
!..UPWIND SCHEME FOR X-DIRECTION VOLUME FLUX COMPONENT
				IF (P(I,J,1) .LT. 0.0) THEN   !CHANGE .LE. TO .LT.
					IF (DP(IP1,J,1) .LT. GX .OR.					&
										DZ(IP1,J,2) .LT. GX) THEN
						ADVX = L%R21(I,J)*(-P(I,J,1)**2/DP(I,J,1))
					ELSE
						ADVX = L%R21(I,J)						&
								*(P(IP1,J,1)**2/DP(IP1,J,1)			&
								-P(I,J,1)**2/DP(I,J,1))
					ENDIF
!*					IF (DP(IP1,J,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						ADVX = L%R21(I,J)*(P(IP1,J,1)**2			&
!*								/ DP(IP1,J,1)-P(I,J,1)**2/DP(I,J,1))
!*					ENDIF
				ELSE
					IF (DP(IM1,J,1) .LT. GX .OR.					&
											DZ(I,J,2) .LT. GX) THEN
						ADVX = L%R21(I,J)*(P(I,J,1)**2/DP(I,J,1))
					ELSE
						ADVX = L%R21(I,J)*(P(I,J,1)**2/DP(I,J,1)	&
										-P(IM1,J,1)**2/DP(IM1,J,1))
					ENDIF
!*					IF (DP(IM1,J,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						ADVX = L%R21(I,J)*(P(I,J,1)**2/DP(I,J,1)	&
!*								- P(IM1,J,1)**2/DP(IM1,J,1))
!*					ENDIF
				ENDIF

!..UPWIND SCHEME FOR Y-DIRECTION VOLUME FLUX COMPONENT
				IF (XQQ .LT. 0.0) THEN
!*					XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)			&
!*									+ Q(I,J,1) + Q(IP1,J,1))

					IF (DZ(I,JP1,2) .LT. GX .OR.					&
										DZ(IP1,JP1,2) .LT. GX) THEN
						ADVY = L%R0(I,J)*(-P(I,J,1)*XQQ/DP(I,J,1))
					ELSEIF (DP(I,JP1,1) .LT. GX) THEN
						ADVY = L%R0(I,J)*(-P(I,J,1)*XQQ/DP(I,J,1))
					ELSE
						XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)		&
									+ Q(I,J,1) + Q(IP1,J,1))
						ADVY = L%R0(I,J)							&
								*(P(I,JP1,1)*XQE/DP(I,JP1,1)		&
								-P(I,J,1)*XQQ/DP(I,J,1))
					ENDIF
!*					IF (DP(I,JP1,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						XQE = 0.25*(Q(I,JP1,1) + Q(IP1,JP1,1)		&
!*									+ Q(I,J,1) + Q(IP1,J,1))
!*						ADVY = L%R0(I,J)*(P(I,JP1,1)*XQE			&
!*								/ DP(I,JP1,1)-P(I,J,1)*XQQ/DP(I,J,1))
!*					ENDIF
				ELSE
!*					XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)			&
!*									+ Q(I,J-2,1) + Q(IP1,J-2,1))

					IF (DZ(I,JM1,2) .LT. GX .OR.					&
										DZ(IP1,JM1,2) .LT. GX) THEN
						ADVY = L%R0(I,J)*(P(I,J,1)*XQQ/DP(I,J,1))
					ELSEIF (DP(I,JM1,1) .LT. GX) THEN
						ADVY = L%R0(I,J)*(P(I,J,1)*XQQ/DP(I,J,1))
					ELSE
						XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)		&
									+ Q(I,J-2,1) + Q(IP1,J-2,1))
						ADVY = L%R0(I,J)*(P(I,J,1)*XQQ/DP(I,J,1)	&
										-P(I,JM1,1)*XQE/DP(I,JM1,1))
					ENDIF
!*					IF (DP(I,JM1,1).GT.GX .AND. DP(I,J,1).GT.GX) THEN
!*						XQE = 0.25*(Q(I,JM1,1) + Q(IP1,JM1,1)		&
!*									+ Q(I,J-2,1) + Q(IP1,J-2,1))
!*						ADVY = L%R0(I,J)*(P(I,J,1)*XQQ/DP(I,J,1)	&
!*								- P(I,JM1,1)*XQE/DP(I,JM1,1))
!*					ENDIF
				ENDIF
				CALL WAVE_BREAKING (L,XP,XQ,ADVX,ADVY,I,J,1)
				XP = XP - ADVX - ADVY
			ENDIF
		    IF (IFRIC.NE.1) XP = XP/(1.0+FF)

!........LIMIT THE DISCHARGE
			IF (ABS(XP) .LT. EPS) XP = 0.0
			PQ_LIMIT = V_LIMIT*DD
			IF (XP .GT. PQ_LIMIT) XP = PQ_LIMIT
			IF (XP .LT. -PQ_LIMIT) XP =-PQ_LIMIT

			P(I,J,2) = XP
			IF (IDIM.EQ.1) P(I,:,2) = XP

			ELSE
			P(I,J,2) = 0
			ENDIF  !--<<<<X<<<<
	     ENDDO
	  ENDDO
       !7!

!.................................................................
!.....SOLVE MOMENTUM EQUATION IN Y DIRECTION
      IF (IDIM.NE.1) THEN
      IS = 2
	  JS = 2
	  IE = IX
	  JE = JY-1
!	  JE = JY
	  IF (LAYID .EQ. 1) THEN
	      IS = 1
		  JS = 1
		  IE = IX-1
		  JE = JY-1
	  ENDIF
      DO I = IS, IE          !3
		 IP1 = I+1
	  	 IM1 = I-1
		 IF (IP1 .GE. IX) IP1 = IX
		 IF (IM1.LE.1) IM1 = 1
		 DO J = JS, JE        !3
			NOFLUX = 0    ! FLAG: 0-CALCULATE Y FLUX; 1-DON'T CALCULATE Y FLUX
			JP1 = J+1
			JM1 = J-1
			IF (JP1 .GE. JY) JP1 = JY
			IF (JM1.LE.1) JM1 = 1
!..GUESS THE RUNUP NEVER REACH AT HZ(I,J)=ELMAX
			IF (HZ(I,J).LE.ELMAX .OR. HQ(I,J).LE.ELMAX) THEN
				Q(I,J,2) = 0.0
				NOFLUX = 1
			ELSE
!..MOVING BOUNDARY
				IF (DZ(I,J,2).GT.GX .AND. DZ(I,JP1,2).GT.GX) THEN
					DD = DQ(I,J,2)
					DF = DQ(I,J,1)
				ELSEIF (DZ(I,J,2).GT.GX .AND. DZ(I,JP1,2).LE.GX	&
						.AND. HZ(I,JP1) + Z(I,J,2).GT.GX) THEN
!*     				DD = HZ(I,JP1) + Z(I,J,2)
					DD = HQ(I,J) + 0.5*(Z(I,J,2)+Z(I,JP1,2))
					DF = DD
				ELSEIF (DZ(I,J,2) .LE. GX .AND. DZ(I,JP1,2).GT.GX &
						.AND. HZ(I,J)+Z(I,JP1,2).GT.GX) THEN
!*     				DD = HZ(I,J) + Z(I,JP1,2)
					DD = HQ(I,J) + 0.5*(Z(I,J,2)+Z(I,JP1,2))
					DF = DD
				ELSE
	 				Q(I,J,2) = 0.0
           			NOFLUX = 1
				ENDIF

				IF (DD .LT. GX) THEN
					Q(I,J,2) = 0.0
					NOFLUX = 1
				ENDIF
			ENDIF

!...........COMPUTE Y FLUX WHEN NOFLUX.NE.1
            IF (NOFLUX .NE. 1) THEN    !--<<<<Y<<<<

!..CALCULATE Y-DIRECTION LINEAR TERMS
!..ESTABLISH A LOWER BOUND
			IF (DF .LT. 1.0E-3) THEN
				DF = 1.0E-3
			ENDIF

!..FIXED FRICTIONAL COEFFICIENT - MANNING
            XPP = 0.25*(P(I,J,1)+P(I,JP1,1)+P(IM1,J,1)+P(IM1,JP1,1))
			IF (IFRIC .EQ. 2) FM = L%FRIC_VCOEF(I,J) !VAR MANNING'S N.
			IF (IFRIC .NE. 1) THEN
!!!				XPP = 0.25*(P(I,J,1)+P(I,JP1,1)						&
!!!									+P(IM1,J,1)+P(IM1,JP1,1))  !7!
				FF = CONST*FM*FM*SQRT(Q(I,J,1)**2					&
										+ XPP**2)/DF**2.3333333
			ELSE
				FF = 0.0
			ENDIF

!...........SOLVE LINEAR VERSION OF MOMENTUM EQUATION IN Y DIRECTION
!!7!->  MODIFIED OR OLD SCHEME FOR MOMENTUM EQUATION
			IF (MOD_SCM .EQ. 0 .AND. I.LT.IX) THEN
				XQ = Q(I,J,1) - L%R4(I,J)*DD*(Z(I,JP1,2)-Z(I,J,2))	&
						- L%R4(I,J)  * DD / 12.						&
						* ((Z(IP1,JP1,2)-2*Z(I,JP1,2)+Z(IM1,JP1,2))	&
						- (Z(IP1,J,2)-2*Z(I,J,2)+Z(IM1,J,2)))
			ELSE
				XQ = Q(I,J,1) - L%R4(I,J)*DD*(Z(I,JP1,2)-Z(I,J,2))
			ENDIF
			IF (IFRIC.NE.1) XQ = XQ - FF*Q(I,J,1)

!!! CORIOLIS EFFECT
            XQ = XQ - 4.0*L%R5(I,J)*XPP

			IF (LIN_CHK.EQ.1 .AND. DQ(I,J,1).GE.GX .AND.			&
					(I.GT.ISPANS .AND. I.LT.ISPANE .AND.			&
					J.GT.JSPANS .AND. J.LT.JSPANE) ) THEN

!!!				XPP = 0.25*(P(I,J,1)+P(I,JP1,1)						&
!!!									+P(IM1,J,1)+P(IM1,JP1,1))  !7!
				ADVX = 0.0
				ADVY = 0.0
!..UPWIND SCHEME FOR Y-DIRECTION VOLUME FLUX COMPONENT
				IF (Q(I,J,1) .LT. 0.0) THEN
					IF (DQ(I,JP1,1) .LT. GX .OR.					&
										DZ(I,JP1,2) .LT. GX) THEN
						ADVY = L%R0(I,J)*(-Q(I,J,1)**2/DQ(I,J,1))
					ELSE
						ADVY = L%R0(I,J)							&
								*(Q(I,JP1,1)**2/DQ(I,JP1,1)			&
								-Q(I,J,1)**2/DQ(I,J,1))
					ENDIF
!*					IF (DQ(I,JP1,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						ADVY = L%R0(I,J)*(Q(I,JP1,1)**2			&
!*								/ DQ(I,JP1,1)-Q(I,J,1)**2/DQ(I,J,1))
!*					ENDIF
				ELSE
					IF (DQ(I,JM1,1) .LT. GX .OR.					&
											DZ(I,J,2) .LT. GX) THEN
						ADVY = L%R0(I,J)*(Q(I,J,1)**2/DQ(I,J,1))
					ELSE
						ADVY = L%R0(I,J)*(Q(I,J,1)**2/DQ(I,J,1)	&
										-Q(I,JM1,1)**2/DQ(I,JM1,1))
					ENDIF
!*					IF (DQ(I,JM1,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*					ADVY = L%R0(I,J)*(Q(I,J,1)**2/DQ(I,J,1)		&
!*							- Q(I,JM1,1)**2/DQ(I,JM1,1))
!*					ENDIF
				ENDIF

!..UPWIND SCHEME FOR X-DIRECTION VOLUME FLUX COMPONENT
				IF (XPP .LT. 0.0) THEN
!*					XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)			&
!*									+ P(I,J,1) + P(I,JP1,1))

					IF (DZ(IP1,J,2) .LT. GX .OR.					&
										DZ(IP1,JP1,2) .LT. GX) THEN
						ADVX = L%R22(I,J)*(-Q(I,J,1)*XPP/DQ(I,J,1))
					ELSEIF (DQ(IP1,J,1) .LT. GX) THEN
						ADVX = L%R22(I,J)*(-Q(I,J,1)*XPP/DQ(I,J,1))
					ELSE
						XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)		&
									+ P(I,J,1) + P(I,JP1,1))
						ADVX = L%R22(I,J)						&
								*(Q(IP1,J,1)*XPE/DQ(IP1,J,1)		&
								-Q(I,J,1)*XPP/DQ(I,J,1))
					ENDIF
!*					IF (DQ(IP1,J,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						XPE = 0.25*(P(IP1,J,1) + P(IP1,JP1,1)		&
!*									+ P(I,J,1) + P(I,JP1,1))
!*						ADVX = L%R22(I,J)*(Q(IP1,J,1)*XPE		&
!*								/ DQ(IP1,J,1)-Q(I,J,1)*XPP/DQ(I,J,1))
!*					ENDIF
				ELSE
					XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)			&
									+ P(I-2,J,1) + P(I-2,JP1,1))

					IF (DZ(IM1,J,2) .LT. GX .OR.					&
										DZ(IM1,JP1,2) .LT. GX) THEN
						ADVX = L%R22(I,J)*(Q(I,J,1)*XPP/DQ(I,J,1))
					ELSEIF (DQ(I-1,J,1) .LT. GX) THEN
						ADVX = L%R22(I,J)*(Q(I,J,1)*XPP/DQ(I,J,1))
					ELSE
						XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)		&
									+ P(I-2,J,1) + P(I-2,JP1,1))
						ADVX = L%R22(I,J)						&
								*(Q(I,J,1)*XPP/DQ(I,J,1)			&
								-Q(IM1,J,1)*XPE/DQ(IM1,J,1))
					ENDIF
!*					IF (DQ(IM1,J,1).GT.GX .AND. DQ(I,J,1).GT.GX) THEN
!*						XPE = 0.25*(P(IM1,J,1) + P(IM1,JP1,1)		&
!*									+ P(I-2,J,1) + P(I-2,JP1,1))
!*						XQ = XQ-L%R22(I,J)*(Q(I,J,1)*XPP/DQ(I,J,1)	&
!*								- Q(IM1,J,1)*XPE/DQ(IM1,J,1))
!*					ENDIF
				ENDIF
				CALL WAVE_BREAKING (L,XP,XQ,ADVX,ADVY,I,J,2)
				XQ = XQ - ADVX - ADVY
			ENDIF
		    IF (IFRIC.NE.1) XQ = XQ/(1.0+FF)
!...........LIMIT THE DISCHARGE
			IF (ABS(XQ) .LT. EPS) XQ = 0.0
			PQ_LIMIT = V_LIMIT*DD
	        IF (XQ .GT. PQ_LIMIT) XQ = PQ_LIMIT
			IF (XQ .LT. -PQ_LIMIT) XQ =-PQ_LIMIT

			Q(I,J,2) = XQ

			ELSE
			Q(I,J,2) = 0.0
            ENDIF   !--<<<<Y<<<<
		  ENDDO
      ENDDO !7!

	  ENDIF

      L%M = P
	  L%N = Q

!	  CALL WAVE_BREAKING (L)

      RETURN
      END


!-----------------------------------------------------------------------
      SUBROUTINE WAVE_BREAKING (L,XP,XQ,ADVX,ADVY,I,J,DIRECTION)
!DESCRIPTION:
!	  #. DETECT WAVE BREAKING AND USE ATIFICIAL DAMPING COEF FOR DAMPING;
!NOTES:
!	  #. CREATED ON FEB 24 2009 (XIAOMING WANG, GNS)
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL ADVX,ADVY
	  REAL NU
	  INTEGER DIRECTION
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!.....CHECK WAVE BREAKING CONDITION
	  DT = L%DT
!*	  IF (L%LAYCORD.EQ.1) THEN
!*	     DX = L%DX
!*	     DY = L%DY
!*	  ELSE
!*	     DX = L%DX*RAD_MIN*R_EARTH
!*		 DY = L%DY*RAD_MIN*R_EARTH
!*	  ENDIF
	  DELTA = 6.5
!*	  DO I = 3,L%NX-2
!*		 IP1 = I+1
!*		 IM1 = I-1
!*		 DO J = 3,L%NY-2
			JP1 = J+1
			JM1 = J-1
			C = SQRT(GRAV*L%DZ(I,J,2))
			DZDT_F = 0.08*C
			DZDT_I = 0.65*C
			DZDT = ABS(L%Z(I,J,2)-L%Z(I,J,1))/L%DT
!*		    IF (DZDT .GT. DZDT_I) THEN
!*			   T_B = 8.0*C/GRAV
!*			   !SETUP THRESHOLD TO EVOKE WATER BREAKING MODEL
!*			   DZDT_B = DZDT_I
!*			   L%MASK(I,J) = NINT(T_B/L%DT)+1
!*			ENDIF

!			IF (DZDT .GT. DZDT_F) THEN
!			   COEF = (DZDT_I-DZDT_F)/(DZDT_I-DZDT_F+DZDT-DZDT_F)
!			   IF (DIRECTION.EQ.1) ADVX = COEF*ADVX
!			   IF (DIRECTION.EQ.2) ADVY = COEF*ADVY
!			ENDIF

			   IF (DZDT .GT. DZDT_F) THEN
!				  COEF = (DZDT_I-DZDT_F)/(DZDT_I-DZDT_F+DZDT-DZDT_F)
				  REF = (DZDT-DZDT_F)/(DZDT_I-DZDT_F)
				  A = 1.0
				  COEF = EXP(-A*REF)
!				  IF (DZDT .GT. DZDT_I) THEN
!			      COEF = (DZDT_I-DZDT_F)/(DZDT_I-DZDT_F+DZDT-DZDT_F)*0.5
!				  ENDIF
				  IF (DIRECTION.EQ.1) THEN
				     ADVX = COEF*ADVX
			         ADVY = COEF*ADVY
				  ENDIF
				  IF (DIRECTION.EQ.2) THEN
				     ADVX = COEF*ADVX
			         ADVY = COEF*ADVY
				  ENDIF
			   ENDIF



!			   IF (L%MASK(I,J).GT.0) THEN
!*			   IF (DZDT .GT. DZDT_F) THEN
!*				  L%T_BREAK(I,J) = T_B
!*				  WRITE (*,*) 'WAVE BREAKS AT:',L%X(I),L%Y(J)


!*			   IF (L%T_BREAK(I,J) .GT. ZERO) THEN
!*				  DZDT_B = DZDT_I + (T_B-L%T_BREAK(I,J))/T_B		&
!*													*(DZDT_F-DZDT_I)
!*				  L%T_BREAK(I,J) = L%T_BREAK(I,J)-L%DT
!*			   ENDIF
!*			   IF (DZDT .GT. 2*DZDT_B) THEN
!*			      B = DELTA
!*				  NU = B*DZDT
!*			   ELSEIF (DZDT.LE.2*DZDT_B .AND. DZDT.GT.DZDT_B) THEN
!*				  B = DELTA*(DZDT/DZDT_B-1.0)
!*				  NU = B*DZDT
!*			   ELSE
!*			      B = 0.0
!*				  NU = 0.0
!*			   ENDIF

!*			   IF (NU .GT. ZERO) THEN
!*				  PXX = (L%M(IP1,J,1)-2.0*L%M(I,J,1)+L%M(IM1,J,1))/DX**2
!*				  QXX = (L%N(IP1,J,1)-2.0*L%N(I,J,1)+L%N(IM1,J,1))/DX**2
!*				  PYY = (L%M(I,JP1,1)-2.0*L%M(I,J,1)+L%M(I,JM1,1))/DY**2
!*				  QYY = (L%N(I,JP1,1)-2.0*L%N(I,J,1)+L%N(I,JM1,1))/DY**2
!*				  PXY = ((L%M(I,J,1)-L%M(IM1,J,1))					&
!*							-(L%M(I,JM1,1)-L%M(IM1,JM1,1)))/(DX*DY)
!*				  QXY = ((L%N(I,J,1)-L%N(I,JM1,1))					&
!*							-(L%N(IM1,J,1)-L%N(IM1,JM1,1)))/(DX*DY)
!*				  RBX = NU*(PXX+0.5*(PYY+QXY))
!*				  RBY = NU*(QYY+0.5*(QXX+PXY))
!				  RBX = NU*(PXX)
!				  RBY = NU*(QYY)
!*				  IF (DIRECTION.EQ.1) XP = XP + L%DT*RBX
!*				  IF (DIRECTION.EQ.2) XQ = XQ + L%DT*RBY
!*			   ENDIF
!				  L%MASK(I,J) = L%MASK(I,J) - 1
!			   ENDIF
!			ENDIF
!*		 ENDDO
!*	  ENDDO

	  RETURN
	  END
