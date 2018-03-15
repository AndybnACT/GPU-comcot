!----------------------------------------------------------------------
	  SUBROUTINE ALPHA_CALC (LO,LA)
!......................................................................
!DESCRIPTION:
!	  #. CALL SUBROUTINE ALPHA() TO DETERMINE HIDDEN GRIDS FOR THE 
!		 PURPOSE OF DISPERSION IMPROVEMENT
!	  #. LO%LAYGOV: 
!				2 - LINEAR SWE WITH DISPERSION-IMPROVED SCHEME;
!	  			3 - NONLINEAR SWE WITH DISPERSION-IMPROVED SCHEME;		
!NOTES:
!	  #. CREATED ON JAN 22 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN.27 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (LAYER), DIMENSION(NUM_GRID)  :: LA
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
	  IF (LO%LAYGOV .EQ. 2 .OR. LO%LAYGOV .EQ. 3) THEN
	  WRITE (*,*) '    CALCULATING COEFICIENTS FOR DISPERSION-IMPROVED SCHEME......'      
	  ENDIF
      IF (LO%LAYSWITCH .EQ. 0) THEN
		 IF (LO%LAYGOV .GE. 2) CALL ALPHA (LO) 	     
	  END IF
	  DO I=1,NUM_GRID
         IF (LA(I)%LAYSWITCH .EQ. 0) THEN
			IF (LA(I)%LAYGOV .GE. 2) CALL ALPHA (LA(I)) 	     
		 END IF
	  END DO

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE ALPHA (LO)
!.....REF: YOON, S.B. (2002), JOURNAL OF GEOPHYSICAL RESEARCH, VOL.107,
!	       NO.C10,3140
!.....DETERMINE COEFFICIENT ALPHA SO AS TO IMPROVE NUMERICAL DISPERSION
!     SUBROUTINE IS DESIGNED TO INCLUDE DISPERSION EFFECT IN SWE
!     ASSUME 'SQUARE' GRID CELLS (I.E., DX=DY)
!.....CREATED ON DEC 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!     UPDATED ON DEC 15, 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  REAL DX, DY, DT, DXS, DYS
	  REAL H(LO%NX,LO%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      DX = LO%DX
	  DY = LO%DY
	  DT = LO%DT
	  H = LO%H

      DO I = 1,LO%NX
	     DO J = 1,LO%NY
		    IF (H(I,J) .GT. GX) THEN
			   !ADJUSTED GRID SIZE IN METERS
		       DELTA = SQRT(4.0*H(I,J)**2+GRAV*H(I,J)*(DT**2))  
               IF (LO%LAYCORD .EQ. 0) THEN
				  LO%ALPHA(I,J) = DELTA/(R_EARTH*COS(LO%Y(J)		&
										*RAD_DEG)*LO%DX*RAD_MIN)
			   ELSE
				  LO%ALPHA(I,J) = DELTA/LO%DX
			   ENDIF
			ELSE
			   LO%ALPHA(I,J) = 1.0
			ENDIF

!SETUP LIMIT ON ALPHA
			A_LIMIT = 10.0
			IF (LO%ALPHA(I,J) .GT. A_LIMIT) LO%ALPHA(I,J) = A_LIMIT

!CALCULATE THE SIZE OF HIDDEN GRIDS
			IF (LO%LAYCORD.EQ.0) THEN
!*			   DYS = LO%DEL_Y(J)*RAD_MIN*LO%ALPHA(I,J)
!*			   DXS = LO%DX*RAD_MIN*LO%ALPHA(I,J)
			   DXS = LO%DX*LO%ALPHA(I,J)
			   DYS = LO%DEL_Y(J)*LO%ALPHA(I,J)
			ELSE
			   DXS = LO%DX*LO%ALPHA(I,J)
			   DYS = LO%DY*LO%ALPHA(I,J)
			ENDIF

!CALCULATE COEF. OF 4-PTS LAGRANGE CUBIC INTERPOLATION
			A = LO%ALPHA(I,J)
			SHIFT = 0.5*(A+1.0)
            N = FLOOR(SHIFT)
			BETA = SHIFT - N
			B1 = 1.0-BETA
			B2 = BETA

			DX = LO%DX
			XF = DX+BETA*DX
			XB = DX+(1.0-BETA)*DX
			X1 = 0.0*DX
			X2 = 1.0*DX
			X3 = 2.0*DX
			X4 = 3.0*DX


			LO%CF(I,J,1) = ((XF-X2)*(XF-X3)*(XF-X4))				&
									/((X1-X2)*(X1-X3)*(X1-X4))
			LO%CF(I,J,2) = ((XF-X1)*(XF-X3)*(XF-X4))				&
									/((X2-X1)*(X2-X3)*(X2-X4))
			LO%CF(I,J,3) = ((XF-X1)*(XF-X2)*(XF-X4))				&
									/((X3-X1)*(X3-X2)*(X3-X4))
			LO%CF(I,J,4) = ((XF-X1)*(XF-X2)*(XF-X3))				&
									/((X4-X1)*(X4-X2)*(X4-X3))

			LO%CB(I,J,1) = ((XB-X2)*(XB-X3)*(XB-X4))				&
									/((X1-X2)*(X1-X3)*(X1-X4))
			LO%CB(I,J,2) = ((XB-X1)*(XB-X3)*(XB-X4))				&
									/((X2-X1)*(X2-X3)*(X2-X4))
			LO%CB(I,J,3) = ((XB-X1)*(XB-X2)*(XB-X4))				&
									/((X3-X1)*(X3-X2)*(X3-X4))
			LO%CB(I,J,4) = ((XB-X1)*(XB-X2)*(XB-X3))				&
									/((X4-X1)*(X4-X2)*(X4-X3))

!CREATE MASK FOR DISPERSION CALCULATION:
!	1 - WITH DISPERSION-IMPROVED SCHEME; 
!	0 - WITHOUT DISPERSION-IMPROVED SCHEME;
			LO%MASK(:,:) = 0
			N_LIMIT = NINT(A_LIMIT)
			IF (I.GT.N_LIMIT .AND. I.LT.LO%NX-N_LIMIT) THEN
			   IF (J.GT.N_LIMIT .AND. J.LT.LO%NY-N_LIMIT) THEN
				  LO%MASK(I,J) = 1
			   ENDIF
			ENDIF
			N_LIMIT = NINT(2.0*A_LIMIT)
			IF (I.GT.N_LIMIT .AND. I.LT.LO%NX-N_LIMIT) THEN
			   IF (J.GT.N_LIMIT .AND. J.LT.LO%NY-N_LIMIT) THEN
				  IF (LO%H(I,J) .LT. GX) THEN
					 DO KI = I-N_LIMIT,I+N_LIMIT
						DO KJ = J-N_LIMIT,J+N_LIMIT
						   LO%MASK(KI,KJ) = 0
						ENDDO
					 ENDDO
				  ENDIF
			   ENDIF
			ENDIF

!......................................................................
!CALCULATE COEF FOR YOON'S METHOD (2002)
			A = LO%ALPHA(I,J)
			IF (A .GT. 3.0) A = 3.0

			IF (LO%LAYCORD .EQ. 0) THEN
		       LO%A1X(I,J) = (A**2-1)/(24.*LO%DX*RAD_MIN)
			   LO%A2X(I,J) = 3.*(9.-A**2)/(24.*LO%DX*RAD_MIN)
			   LO%A1Y(I,J) = (A**2-1)/(24.*LO%DEL_Y(J)*RAD_MIN)
			   LO%A2Y(I,J) = 3.*(9.-A**2)/(24.*LO%DEL_Y(J)*RAD_MIN)
			ELSE
		       LO%A1X(I,J) = (A**2-1)/(24.*LO%DX)
			   LO%A2X(I,J) = 3.*(9.-A**2)/(24.*LO%DX)
			   LO%A1Y(I,J) = (A**2-1)/(24.*LO%DEL_Y(J))
			   LO%A2Y(I,J) = 3.*(9.-A**2)/(24.*LO%DEL_Y(J))
            ENDIF
		 ENDDO
      ENDDO


	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE UVW (U,U0,V,V0,WX,WY,WX0,WY0,LO)
!.....REF: YOON, S.B. (2002), JOURNAL OF GEOPHYSICAL RESEARCH, VOL.107,
!		   NO.C10,3140
!.....DETERMINE COEFFICIENT ALPHA SO AS TO IMPROVE NUMERICAL DISPERSION
!     SUBROUTINE DESIGNED TO INCLUDED DISPERSION IN SWE
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER) :: LO
	  REAL DX, DY, DT, DXS, DYS, G
	  REAL H(LO%NX,LO%NY)
	  REAL U(LO%NX,LO%NY),U0(LO%NX,LO%NY)  !0 MEANS VALUE AT TIME STEP N-1
	  REAL V(LO%NX,LO%NY),V0(LO%NX,LO%NY)
	  REAL WX0(LO%NX,LO%NY),WY0(LO%NX,LO%NY)  !PQ AT CELL CORNERS
	  REAL WX(LO%NX,LO%NY),WY(LO%NX,LO%NY) !PQ AT CELL EDGES 
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
	  DATA RAD/0.01745329252/

	  G = GRAV
      DX = LO%DX
	  DY = LO%DY
	  DT = LO%DT
	  H = LO%H

      DO I = 1,LO%NX
	     IP1 = I + 1
		 IM1 = I - 1
		 IF(IP1.GE.LO%NX) IP1 = LO%NX
		 IF(IM1.LE.1) IM1 = 1
	     DO J = 1,LO%NY
		    JP1 = J + 1
		    JM1 = J - 1
		    IF(JP1.GE.LO%NY) JP1 = LO%NY
		    IF(JM1.LE.1) JM1 = 1
			IF (LO%HP(I,J).GT.0.0 .AND. LO%HP(IP1,J).GT.0.0) THEN
			   IF (LO%HQ(I,J).GT.0.0 .AND. LO%HQ(I,JP1).GT.0.0) THEN
			      HP = LO%HP(I,J) 
				  HQ = LO%HQ(I,J) 
				  HPQ = (LO%HP(I,J) + LO%HP(I,JP1) + LO%HQ(I,J)		&
												+ LO%HQ(IP1,J))/4.0

			      U(I,J) = LO%M(I,J,1)**2/HP
				  U0(I,J) = LO%M0(I,J)**2/HP

                  V(I,J) = LO%N(I,J,1)**2/HQ
				  V0(I,J) = LO%N0(I,J)**2/HQ

				  WX(I,J) = (LO%N(I,J,1)+LO%N(IP1,J,1)				&
								+LO%N(I,JM1,1)+LO%N(IP1,JM1,1))		&
								*(LO%M(I,J,1))/(4.0*HP)
				  WY(I,J) = (LO%M(I,J,1)+LO%M(IM1,J,1)				&
								+LO%M(I,JP1,1)+LO%M(IM1,JP1,1))		&
								*(LO%N(I,J,1))/(4.0*HQ)

                  WX0(I,J) = (LO%N0(I,J)+LO%N0(IP1,J)+LO%N0(I,JM1)	&
							+LO%N0(IP1,JM1))*(LO%M0(I,J))/(4.0*HP)
                  WY0(I,J) = (LO%M0(I,J)+LO%M0(IM1,J)+LO%M0(I,JP1)	&
							+LO%M0(IM1,JP1))*(LO%N0(I,J))/(4.0*HQ)

			   ENDIF
            ENDIF
		 ENDDO
      ENDDO

	  RETURN
	  END


!----------------------------------------------------------------------
      SUBROUTINE MASS_S_D (L)
! ....SOLVE CONTINUITY EQUATION IN SPHERICAL COORD. 
!.....ADOPT DISPERSION-IMPROVED NUMERICAL SCHEME
!.....CREATED ON DEC 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!     UPDATED ON DEC.15 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
      REAL ALPHA, DPDX, DQDY, A1, A2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/
	  DATA OSIXTY/0.016666666667/
!
      IS = 2
	  JS = 2
!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
      IDIM = 2   
!.....INTERPOLATION METHOD. 
!	  1 - YOON 2002; 
!	  2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2   

	  IF (L%ID .EQ. 1) THEN
	    IE = L%NX-1    !FOR OUTEST LAYER
	    JE = L%NY-1
	  ELSE
	    IE = L%NX      !FOR INNER LAYER
        JE = L%NY
      END IF
      DO J=JS,JE
        JM1 = J-1
		JM2 = J-2
		JP1 = J+1
		JP2 = J+2
		IF (JM1.LE.1) JM1 = 1
		IF (JM2.LE.1) JM2 = 1
		IF (JP1.GE.L%NY) JP1 = L%NY
		IF (JP2.GE.L%NY) JP2 = L%NY
        DO I=IS,IE
		  IM1 = I-1
		  IM2 = I-2
		  IP1 = I+1
		  IP2 = I+2
		  IF (IM1.LE.1) IM1 = 1
		  IF (IM2.LE.1) IM2 = 1
		  IF (IP1.GE.L%NX) IP1 = L%NX
		  IF (IP2.GE.L%NX) IP2 = L%NX
          IF (L%H(I,J) .GT. GX) THEN
            A = L%ALPHA(I,J)
			SHIFT = 0.5*(A+1.0)
            NN = FLOOR(SHIFT)

		    DXS = L%DX*RAD_MIN*L%ALPHA(I,J)
			DYS = L%DEL_Y(J)*RAD_MIN*L%ALPHA(I,J)

		    IF (L%MASK(I,J) .EQ. 1) THEN
			   IF (INTP.EQ.1) THEN
		          ALPHA = L%ALPHA(I,J)
		          A1X = L%A1X(I,J)   !(ALPHA**2-1)/(24.*L%DX)
			      A1Y = L%A1Y(I,J)
			      A2X = L%A2X(I,J)   !3.*(9.-ALPHA**2)/(24.*L%DX)
			      A2Y = L%A2Y(I,J)
			      DPDX = A1X*(L%M(IP1,J,1)-L%M(IM2,J,1))			&
							+A2X*(L%M(I,J,1)-L%M(IM1,J,1))
                  DQDY = A1Y*(L%N(I,JP1,1)*L%R6(JP1)-L%N(I,JM2,1)	&
							*L%R6(J-2))+A2Y*(L%N(I,J,1)*L%R6(J)		&
							-L%N(I,JM1,1)*L%R6(JM1))
			      ZZ = L%Z(I,J,1)-L%R1(I,J)*DXS*DPDX				&
									-L%R11(I,J)*DYS*DQDY
			   ENDIF
			   IF  (INTP.EQ.2) THEN
                  PF = L%CF(I,J,1)*L%M(IM1+NN,J,1)					&
						+ L%CF(I,J,2)*L%M(I+NN,J,1)					&
						+ L%CF(I,J,3)*L%M(IP1+NN,J,1)				&
						+ L%CF(I,J,4)*L%M(IP2+NN,J,1)
                  PB = L%CB(I,J,1)*L%M(IM1-NN,J,1)					&
						+ L%CB(I,J,2)*L%M(I-NN,J,1)					&
						+ L%CB(I,J,3)*L%M(IP1-NN,J,1)				&
						+ L%CB(I,J,4)*L%M(IP2-NN,J,1)
				  QF = L%CF(I,J,1)*L%N(I,JM1+NN,1)*L%R6(JM1+NN)		&
						+ L%CF(I,J,2)*L%N(I,J+NN,1)*L%R6(J+NN)		&
						+ L%CF(I,J,3)*L%N(I,JP1+NN,1)*L%R6(JP1+NN)	&
						+ L%CF(I,J,4)*L%N(I,JP2+NN,1)*L%R6(JP2+NN)
				  QB = L%CB(I,J,1)*L%N(I,JM1-NN,1)*L%R6(JM1-NN)		&
						+ L%CB(I,J,2)*L%N(I,J-NN,1)*L%R6(J-NN)		&
						+ L%CB(I,J,3)*L%N(I,JP1-NN,1)*L%R6(JP1-NN)	&
						+ L%CB(I,J,4)*L%N(I,JP2-NN,1)*L%R6(JP2-NN)
			      DPDX = (PF-PB)
			      DQDY = (QF-QB)
!*				  ZZ = L%Z(I,J,1)-L%R1(I,J)*DPDX-L%R0(I,J)*DQDY
				  ZZ = L%Z(I,J,1)-L%R1(I,J)*DPDX-L%R11(I,J)*DQDY

			   ENDIF
			ELSE
               ZZ = L%Z(I,J,1)-L%R1(I,J)*(L%M(I,J,1)-L%M(IM1,J,1))	&
							-L%R11(I,J)*(L%N(I,J,1)*L%R6(J)			&
										-L%N(I,JM1,1)*L%R6(JM1))
			ENDIF

            IF (L%INI_SWITCH .EQ. 3 .OR. L%INI_SWITCH.EQ.4) THEN
			   ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))
			ENDIF
			IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
!			IF (ABS(ZZ) .GT. UB) ZZ=SIGN(UB,ZZ) 
			!DEPRESSION CANNOT BE LESS THAN BOTTOM ELEVATION
			IF ( (ZZ+L%H(I,J)) .LE. ZERO ) ZZ = -L%H(I,J) 
            L%Z(I,J,2) = ZZ
		  ELSE
		    L%Z(I,J,2) = 0.0
          ENDIF
       END DO
      END DO

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE MASS_C_D (L)		 
! ....SOLVE CONTINUITY EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!
!NOTES: 
!	  #. ADD SUPPORT FOR DX\=DY (FOR HIGH-LATITUDE, 05/04/2007)
!            RX = L%R
!            RY = L%DT/L%DY

!     SUBROUTINE DESIGNED TO INCCLUDE DISPERSION IN SWE
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL ALPHA, DPDX, DQDY, A1, A2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/

	  RX = L%RX
	  RY = L%RY
	  GRT = L%DT*GRAV
	  
!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
      IDIM = 2   
!.....INTERPOLATION METHOD. 
!		1 - YOON2002; 
!		2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2   

!
      IS = 2
      JS = 2
	  IF(L%ID .EQ. 1)THEN  !OUTTEST LAYER				 
	    IE = L%NX -1
	    JE = L%NY -1
	  ELSE				  ! INNER LAYER
	    IE = L%NX
	    JE = L%NY
	  ENDIF

      IF (IDIM.EQ.2) THEN
      DO J=JS,JE
        JM1 = J-1
		JM2 = J-2
		JP1 = J+1
		JP2 = J+2
		IF (JM1.LE.1) JM1 = 1
		IF (JM2.LE.1) JM2 = 1
		IF (JP1.GE.L%NY) JP1 = L%NY
		IF (JP2.GE.L%NY) JP2 = L%NY
        DO I=IS,IE
		  IM1 = I-1
		  IM2 = I-2
		  IP1 = I+1
		  IP2 = I+2
		  IF (IM1.LE.1) IM1 = 1
		  IF (IM2.LE.1) IM2 = 1
		  IF (IP1.GE.L%NX) IP1 = L%NX
		  IF (IP2.GE.L%NX) IP2 = L%NX

          IF (L%H(I,J) .GT. GX) THEN

            A = L%ALPHA(I,J)
			SHIFT = 0.5*(A+1.0)
            NN = FLOOR(SHIFT)

		    DXS = L%DX*L%ALPHA(I,J)
			DYS = L%DY*L%ALPHA(I,J)

		    IF (L%MASK(I,J).EQ.1) THEN
			   IF (INTP.EQ.1) THEN
		          A1X = L%A1X(I,J)   !(ALPHA**2-1)/(24.*L%DX)
			      A1Y = L%A1Y(I,J)
			      A2X = L%A2X(I,J)   !3.*(9.-ALPHA**2)/(24.*L%DX)
			      A2Y = L%A2Y(I,J)
			      DPDX = A1X*(L%M(IP1,J,1)-L%M(IM2,J,1))			&
							+A2X*(L%M(I,J,1)-L%M(IM1,J,1))
                  DQDY = A1Y*(L%N(I,JP1,1)-L%N(I,JM2,1))			&
							+A2Y*(L%N(I,J,1)-L%N(I,JM1,1))
			      ZZ = L%Z(I,J,1)-L%DT*(DPDX+DQDY)
			   ENDIF
			   IF  (INTP.EQ.2) THEN
                  PF = L%CF(I,J,1)*L%M(IM1+NN,J,1)					&
						+ L%CF(I,J,2)*L%M(I+NN,J,1)					&
						+ L%CF(I,J,3)*L%M(IP1+NN,J,1)				&
						+ L%CF(I,J,4)*L%M(IP2+NN,J,1)
                  PB = L%CB(I,J,1)*L%M(IM1-NN,J,1)					&
						+ L%CB(I,J,2)*L%M(I-NN,J,1)					&
						+ L%CB(I,J,3)*L%M(IP1-NN,J,1)				&
						+ L%CB(I,J,4)*L%M(IP2-NN,J,1)
				  QF = L%CF(I,J,1)*L%N(I,JM1+NN,1)					&
						+ L%CF(I,J,2)*L%N(I,J+NN,1)					&
						+ L%CF(I,J,3)*L%N(I,JP1+NN,1)				&
						+ L%CF(I,J,4)*L%N(I,JP2+NN,1)
				  QB = L%CB(I,J,1)*L%N(I,JM1-NN,1)					&
						+ L%CB(I,J,2)*L%N(I,J-NN,1)					&
						+ L%CB(I,J,3)*L%N(I,JP1-NN,1)				&
						+ L%CB(I,J,4)*L%N(I,JP2-NN,1)
				  DP = (PF-PB)
			      DQ = (QF-QB)
			      ZZ = L%Z(I,J,1) - (L%DT/DXS*DP+L%DT/DYS*DQ)
			   ENDIF
			ELSE
               ZZ = L%Z(I,J,1) - RX*(L%M(I,J,1)-L%M(IM1,J,1))		&
								- RY*(L%N(I,J,1)-L%N(I,JM1,1))
			ENDIF
			IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4) THEN
			   ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))  
			ENDIF
			IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
!*			IF (ABS(ZZ) .GT. UB) ZZ=SIGN(UB,ZZ) 
			!DEPRESSION CANNOT BE LESS THAN BOTTOM ELEVATION 
			IF ( (ZZ+L%H(I,J)) .LE. ZERO ) ZZ = -L%H(I,J) 
            L%Z(I,J,2) = ZZ
		  ELSE
            L%Z(I,J,2) = 0.0
          ENDIF
        END DO
      END DO
      ENDIF

!.....CODE FOR 1-D STUDY     
      IF (IDIM.EQ.1) THEN
	     L%N(:,:,:) = 0.0
	     J = FLOOR(L%NY/2.)
	     DO I=IS,IE
		    IM1 = I-1
			IM2 = I-2
			IP1 = I+1
			IP2 = I+2
			IF (IM1.LE.1) IM1 = 1
			IF (IM2.LE.1) IM2 = 1
			IF (IP1.GE.L%NX) IP1 = L%NX
			IF (IP2.GE.L%NX) IP2 = L%NX

            IF (L%H(I,J) .GT. GX) THEN
               A = L%ALPHA(I,J)
			   SHIFT = 0.5*(A+1.0)
               NN = FLOOR(SHIFT)

	           DXS = L%DX*L%ALPHA(I,J)
			   BETA = SHIFT - NN

			   IF (L%MASK(I,J).EQ.1) THEN
			      IF (INTP.EQ.1) THEN
		             A1X = L%A1X(I,J)   !(ALPHA**2-1)/(24.*L%DX)
			         A2X = L%A2X(I,J)   !3.*(9.-ALPHA**2)/(24.*L%DX)
			         DPDX = A1X*(L%M(IP1,J,1)-L%M(IM2,J,1))			&
							+A2X*(L%M(I,J,1)-L%M(IM1,J,1))
			         ZZ = L%Z(I,J,1)-L%DT*DPDX
				  ENDIF
				  IF (INTP.EQ.2) THEN
			         IM2 = IM1-1
                     PF = L%CF(I,J,1)*L%M(IM1+NN,J,1)				&
							+ L%CF(I,J,2)*L%M(I+NN,J,1)				&
							+ L%CF(I,J,3)*L%M(IP1+NN,J,1)			&
							+ L%CF(I,J,4)*L%M(IP2+NN,J,1)
                     PB = L%CB(I,J,1)*L%M(IM1-NN,J,1)				&
							+ L%CB(I,J,2)*L%M(I-NN,J,1)				&
							+ L%CB(I,J,3)*L%M(IP1-NN,J,1)			&
							+ L%CB(I,J,4)*L%M(IP2-NN,J,1)
			         DP = (PF-PB)
			         ZZ = L%Z(I,J,1)-L%DT/DXS*DP
			      ENDIF
			   ELSE
                  ZZ = L%Z(I,J,1)-RX*(L%M(I,J,1)-L%M(IM1,J,1))
			   ENDIF
			   IF (L%INI_SWITCH.EQ.3 .OR. L%INI_SWITCH.EQ.4) THEN
			      ZZ = ZZ-(L%HT(I,J,2)-L%HT(I,J,1))  
			   ENDIF
			   IF (ABS(ZZ) .LT. EPS) ZZ = ZERO
!			   IF (ABS(ZZ) .GT. UB) ZZ=SIGN(UB,ZZ)  
			   !DEPRESSION CANNOT BE LESS THAN BOTTOM ELEVATION
			   IF ( (ZZ+L%H(I,J)) .LE. ZERO ) ZZ = -L%H(I,J) 
               L%Z(I,:,2) = ZZ
		    ELSE
               L%Z(I,:,2) = 0.0
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE MOMT_S_D (L)
! ....SOLVE MOMENTUM EQUATION (LINEAR) IN SPHERICAL COORD.
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL MODSCM
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA TWLVTH/0.08333333333333/
	  DATA RAD/0.01745329252/
	  MODSCM = 0.0
!
      GRT = GRAV*L%DT

!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
      IDIM = 2   
!.....INTERPOLATION METHOD. 
!		1 - YOON2002; 
!		2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2  

      IXM1 = L%NX-1
      JYM1 = L%NY-1
	  IS = 2
	  JS = 2
	  IF (L%ID .EQ. 1) THEN
	    IS = 1
		JS = 1
      END IF
      DO I=IS,IXM1
        IP1 = I+1
        DO J=JS,L%NY
          IF ((L%H(I,J).GT.GX) .AND. (L%H(IP1,J).GT.GX)) THEN
            JM1 = J-1
			JM2 = J-2
			JP1 = J+1
			JP2 = J+2
			IF (JM1 .LT. 1) JM1 = 1
			IF (JM2 .LT. 1) JM2 = 1
            IF (JP1 .GT. L%NY ) JP1 = L%NY
            IF (JP2 .GT. L%NY ) JP2 = L%NY

			HM = L%HP(I,J)
            A = L%ALPHA(I,J)
			SHIFT = 0.5*(A+1.0)
            NN = FLOOR(SHIFT)

			DXS = L%DX*RAD_MIN*L%ALPHA(I,J)
			DYS = L%DEL_Y(J)*RAD_MIN*L%ALPHA(I,J)

            IPN = I+NN
			IMN = I-NN

            NS = NINT(A)
			IF (NS.LT.1) NS = 1
			JPNN = J+NS
			JMNN = J-NS


			IP1PN = IP1+NN
			IM1PN = IM1+NN
			IP1MN = IP1-NN
			IM1MN = IM1-NN
			IP2PN = IP2+NN
			IP2MN = IP2-NN


		    IF (L%MASK(I,J) .EQ. 1) THEN
			   IF (INTP.EQ.1) THEN
		          ALPHA = L%ALPHA(I,J)
		          A1X = L%A1X(I,J) !(ALPHA**2-1)/(24.*L%DX)
			      A2X = L%A2X(I,J) !3.*(9.-ALPHA**2)/(24.*L%DX)
			      DZDX = A1X*(L%Z(I+2,J,2)-L%Z(I-1,J,2))			&
									+ A2X*(L%Z(I+1,J,2)-L%Z(I,J,2))
                  TOT_N = L%N(I,J,1) + L%N(IP1,J,1) + L%N(I,JM1,1)	&
									+ L%N(IP1,JM1,1)
				  R2 = L%R2(I,J)*DXS
                  XM = L%M(I,J,1)-R2*DZDX + L%R3(I,J)*TOT_N				  
			      IF (L%MODSCM .EQ. 0) THEN
			         IF (ALPHA .LT. 1.732) THEN
			            GAMMA = ALPHA**2
			            DZFDX = A1X*(L%Z(I+2,J+1,2)-L%Z(I-1,J+1,2))	&
								+A2X*(L%Z(I+1,J+1,2)-L%Z(I,J+1,2))
				        DZBDX = A1X*(L%Z(I+2,J-1,2)-L%Z(I-1,J-1,2))	&
								+A2X*(L%Z(I+1,J-1,2)-L%Z(I,J-1,2))
			         ELSE
			            GAMMA = ALPHA**2/4.0
			            DZFDX = A1X*(L%Z(I+2,J+2,2)-L%Z(I-1,J+2,2))	&
								+A2X*(L%Z(I+1,J+2,2)-L%Z(I,J+2,2))
				        DZBDX = A1X*(L%Z(I+2,J-2,2)-L%Z(I-1,J-2,2))	&
								+A2X*(L%Z(I+1,J-2,2)-L%Z(I,J-2,2))				     
				     ENDIF
			         SCM = GAMMA*R2*TWLVTH					&
											*(DZFDX-2*DZDX+DZBDX)
					 XM = XM - SCM
				  ENDIF
               ENDIF
			   IF (INTP.EQ.2) THEN !4-PTS LAGRANGE INTERPOLATION
                  ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)					&
						+ L%CF(I,J,2)*L%Z(IPN,J,2)					&
						+ L%CF(I,J,3)*L%Z(IP1PN,J,2)				&
						+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				  ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)					&
						+ L%CB(I,J,2)*L%Z(IMN,J,2)					&
						+ L%CB(I,J,3)*L%Z(IP1MN,J,2)				&
						+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
				  DZMDX = (ZMF-ZMB)
                  TOT_N = L%N(I,J,1)+L%N(IP1,J,1)+L%N(I,JM1,1)		&
													+L%N(IP1,JM1,1)
                  XM = L%M(I,J,1)-L%R2(I,J)*DZMDX + L%R3(I,J)*TOT_N		
				  IF (L%MODSCM .EQ. 0) THEN
                     ZUF = L%CF(I,J,1)*L%Z(IM1PN,JPNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JPNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JPNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JPNN,2)
				     ZUB = L%CB(I,J,1)*L%Z(IM1MN,JPNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JPNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JPNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JPNN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IM1PN,JMNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JMNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JMNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JMNN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IM1MN,JMNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JMNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JMNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JMNN,2)
				     DZUDX = (ZUF-ZUB)
				     DZLDX = (ZLF-ZLB)
				     ZADD = DZUDX-2.0*DZMDX+DZLDX
					 GAMMA = L%ALPHA(I,J)**2/(NS**2)
					 SCM = GAMMA*L%R2(I,J)*TWLVTH*ZADD
					 XM = XM - SCM
				  ENDIF
			   ENDIF
			ELSE
		       TOT_N = L%N(I,J,1) + L%N(IP1,J,1) + L%N(I,JM1,1)		&
												+ L%N(IP1,JM1,1)
               XM = L%M(I,J,1)-L%R2(I,J)*(L%Z(IP1,J,2)-L%Z(I,J,2))	&
											+L%R3(I,J)*TOT_N
			   IF (L%MODSCM .EQ. 0) THEN
			      SCM = L%R2(I,J)*TWLVTH*((L%Z(IP1,JP1,2)			&
							-2*L%Z(IP1,J,2)+L%Z(IP1,JM1,2))			&
							-(L%Z(I,JP1,2)-2*L%Z(I,J,2)				&
							+L%Z(I,JM1,2)))
                  XM = XM - SCM
               ENDIF
            ENDIF

            IF (ABS(XM) .LT. EPS) XM = ZERO
            L%M(I,J,2) = XM
          ELSE
		    L%M(I,J,2) = 0.0
		  ENDIF
        END DO
      END DO
!
      DO J=JS,JYM1
        JP1 = J+1
        DO I=IS,L%NX
          IF ((L%H(I,J).GT.ZERO) .AND. (L%H(I,JP1).GT.ZERO)) THEN
            IP1 = I+1
		    IP2 = I+2
		    IM1 = I-1
		    IM2 = I-2
		    IF (IP1 .GE. L%NX) IP1 = L%NX
		    IF (IP2 .GE. L%NX) IP2 = L%NX
		    IF (IM1 .LE. 1) IM1 = 1
		    IF (IM2 .LE. 1) IM2 = 1
              
			HN = L%HQ(I,J)

            A = L%ALPHA(I,J)
			SHIFT = 0.5*(A+1.0)
            NN = FLOOR(SHIFT)

			DXS = L%DX*RAD_MIN*L%ALPHA(I,J)
			DYS = L%DEL_Y(J)*RAD_MIN*L%ALPHA(I,J)

            JPN = J+NN
			JMN = J-NN
			JP1PN = JP1+NN
			JM1PN = JM1+NN
			JP1MN = JP1-NN
			JM1MN = JM1-NN
			JP2PN = JP2+NN
			JP2MN = JP2-NN

            NS = NINT(A)
			IF (NS.LT.1) NS = 1
		    IPNN = I+NS
			IMNN = I-NS

            GAMMA = L%ALPHA(I,J)**2
			SCM = 0.0

		    IF (L%MASK(I,J) .EQ. 1) THEN
			   IF (INTP.EQ.1) THEN
		          ALPHA = L%ALPHA(I,J)
		          A1Y = L%A1Y(I,J) !(ALPHA**2-1)/(24.*L%DY)
			      A2Y = L%A2Y(I,J) !3.*(9.-ALPHA**2)/(24.*L%DY)
			      DZDY = A1Y*(L%Z(I,JP2,2)-L%Z(I,J-1,2))			&
						+A2Y*(L%Z(I,JP1,2)-L%Z(I,J,2))
                  TOT_M = L%M(IM1,J,1)+L%M(IM1,JP1,1)+L%M(I,J,1)	&
													+L%M(I,JP1,1)
				  R4 = L%R4(I,J)*DYS
                  XN = L%N(I,J,1)-R4*DZDY -L%R5(I,J)*TOT_M
			      IF (L%MODSCM .EQ. 0) THEN
			         IF (ALPHA .LT. 1.732) THEN
			            GAMMA = ALPHA**2
			            DZUDY = A1Y*(L%Z(IP1,JP2,2)-L%Z(IP1,JM1,2))	&
								+A2Y*(L%Z(IP1,JP1,2)-L%Z(IP1,J,2))
				        DZLDY = A1Y*(L%Z(IM1,JM2,2)-L%Z(IM1,JM1,2))	&
								+A2Y*(L%Z(IM1,JP1,2)-L%Z(IM1,J,2))
			         ELSE
			            GAMMA = ALPHA**2/4.0
			            DZUDY = A1Y*(L%Z(IP2,JP2,2)-L%Z(IP2,JM1,2))	&
								+A2Y*(L%Z(IP2,JP1,2)-L%Z(IP2,J,2))
				        DZLDY = A1Y*(L%Z(IM2,JP2,2)-L%Z(IM2,JM1,2))	&
								+A2Y*(L%Z(IM2,JP1,2)-L%Z(IM2,J,2))				     
				     ENDIF
			         MODSCM = GAMMA*R4*TWLVTH				&
											*(DZUDY-2*DZDY+DZLDY)
					 XN = XN - MODSCM
                  ENDIF
			   ENDIF
			   IF (INTP.EQ.2) THEN
                  ZMF = L%CF(I,J,1)*L%Z(I,JM1PN,2)					&
						+ L%CF(I,J,2)*L%Z(I,JPN,2)					&
						+ L%CF(I,J,3)*L%Z(I,JP1PN,2)				&
						+ L%CF(I,J,4)*L%Z(I,JP2PN,2)
				  ZMB = L%CB(I,J,1)*L%Z(I,JM1MN,2)					&
						+ L%CB(I,J,2)*L%Z(I,JMN,2)					&
						+ L%CB(I,J,3)*L%Z(I,JP1MN,2)				&
						+ L%CB(I,J,4)*L%Z(I,JP2MN,2)
				  DZMDY = (ZMF-ZMB)
                  TOT_M = L%M(IM1,J,1)+L%M(IM1,JP1,1)+L%M(I,J,1)	&
													+L%M(I,JP1,1)
                  XN = L%N(I,J,1)-L%R4(I,J)*DZMDY-L%R5(I,J)*TOT_M
				  IF (L%MODSCM.EQ.0) THEN
                     ZRF = L%CF(I,J,1)*L%Z(IPNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IPNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IPNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IPNN,JP2PN,2)
				     ZRB = L%CB(I,J,1)*L%Z(IPNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IPNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IPNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IPNN,JP2MN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IMNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IMNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IMNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IMNN,JP2PN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IMNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IMNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IMNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IMNN,JP2MN,2)
                     
				     DZRDY = (ZRF-ZRB)
				     DZLDY = (ZLF-ZLB)
				     ZADD = DZRDY-2.0*DZMDY+DZLDY
					 GAMMA = L%ALPHA(I,J)**2/NS**2

					 SCM = GAMMA*L%R4(I,J)*TWLVTH*ZADD
					 XN = XN - SCM
				  ENDIF
			   ENDIF
			ELSE
               TOT_M = L%M(IM1,J,1)+L%M(IM1,JP1,1)+L%M(I,J,1)		&
												+L%M(I,JP1,1)
               XN = L%N(I,J,1)-L%R4(I,J)*(L%Z(I,JP1,2)-L%Z(I,J,2))	&
												-L%R5(I,J)*TOT_M
			   IF (L%MODSCM .EQ. 0) THEN
			      SCM = L%R4(I,J)*TWLVTH*((L%Z(IP1,JP1,2)			&
							-2*L%Z(I,JP1,2)+L%Z(IM1,JP1,2))			&
							-(L%Z(IP1,J,2)-2*L%Z(I,J,2)				&
							+L%Z(IM1,J,2)))	
				  XN = XN - SCM
               ENDIF
            ENDIF

            IF (ABS(XN) .LT. EPS) XN = ZERO
            L%N(I,J,2) = XN
          ELSE
		    L%N(I,J,2) = 0.0
		  ENDIF
        END DO
      END DO
!
      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE MOMT_C_D (L)
!.....SOLVE MOMENTUM EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!     SUBROUTINE DESIGNED TO INCCLUDE DISPERSION IN SWE
!     OPTION L%LAYGOV ==
!                   2, LINEAR EQUATION WITH DISPERSION CORRECTION
!                   3, NONLINEAR EQUATION WITH DISPERSION CORRECTION
!.....CREATED ON DEC 2007 (XIAOMING WANG, CORNELL UNIVERSITY)
!     UPDATED ON DEC 15, 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL SCM, GRX,GRY,FF,FM
	  REAL U(L%NX,L%NY),U0(L%NX,L%NY)  !0 MEANS VALUE AT TIME STEP N-1
	  REAL V(L%NX,L%NY),V0(L%NX,L%NY)
	  REAL W(L%NX,L%NY),W0(L%NX,L%NY)
	  REAL WX0(L%NX,L%NY),WY0(L%NX,L%NY)  !PQ AT CELL CORNERS
	  REAL WX(L%NX,L%NY),WY(L%NX,L%NY) !PQ AT CELL EDGES
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA TWLVTH/0.08333333333333/

	  SCM = 0.0
	  CONST = 0.5*L%DT*GRAV
	  FM = L%FRIC_COEF
	  FF = 0.0
      
	  DX = L%DX
	  DY = L%DY
	  DT = L%DT
	  GRT = GRAV*L%DT
	  RX = L%RX
	  RY = L%RY
	  GRX = L%GRX
	  GRY = L%GRY

!.....DIMENSIONALITY OF COMPUTATION. 1 - 1D; 2 - 2D
	  IDIM = 2
	  IF (IDIM .EQ. 1) L%N = 0.0
!.....INTERPOLATION METHOD. 
!		0 - USING FIXED GRIDS; 
!		1 - LINEAR; 
!		2 - LAGRANGE 4-POINTS (CUBIC)
	  INTP = 2
!.....CONTROL IF NONLINEAR TERMS USE HIDDEN GRIDS. 
!		0 - NO NONLINEAR;
!		1 - HIDDEN GRIDS, 
!		2 - FIXED GRIDS;
!		3 - FIXED GRIDS (V. COEF.)  
	  NL = 2  
!.....CONTROL IF 5TH-ORDER DERIVATIVES ARE SOLVED. 
!		0 - NO; 
!		1 - YES;	    
	  IT4 = 0    

      IF (L%LAYGOV.EQ.2) NL = 0
	  IF (L%LAYGOV.EQ.3) NL = 2

      IF (NL.NE.0) CALL UVW (U,U0,V,V0,WX,WY,WX0,WY0,L)
!
      IXM1 = L%NX-1
      JYM1 = L%NY-1
	  IF (L%ID .EQ. 1) THEN
	    IS = 1
		JS = 1
	  ELSE
	    IS = 2
		JS = 2
	  END IF

	  IF (IDIM.EQ.2) THEN
	  ! MOMENTUM IN X DIRECTION
      DO I=IS,IXM1
        IP1 = I+1
		IP2 = I+2
		IM1 = I-1
		IM2 = I-2
		IF (IP1 .GE. L%NX) IP1 = L%NX
		IF (IP2 .GE. L%NX) IP2 = L%NX
		IF (IM1 .LE. 1) IM1 = 1
		IF (IM2 .LE. 1) IM2 = 1
        DO J=JS,L%NY
          IF ((L%H(I,J).GT.GX) .AND. (L%H(IP1,J).GT.GX)) THEN
            JM1 = J-1
			JM2 = J-2
			JP1 = J+1
			JP2 = J+2
			IF (JM1 .LT. 1) JM1 = 1
			IF (JM2 .LT. 1) JM2 = 1
            IF (JP1 .GT. L%NY ) JP1 = L%NY
            IF (JP2 .GT. L%NY ) JP2 = L%NY
			HM = L%HP(I,J)+0.5*(L%Z(I,J,2)+L%Z(IP1,J,2))
            !USE BOTTOM FRICTION
			IF (L%FRIC_SWITCH.NE.1) THEN
			   IF (L%FRIC_SWITCH .EQ. 2) FM = L%FRIC_VCOEF(I,J)
			   XQQ = 0.25*(L%N(I,J,1)+L%N(IP1,J,1)+L%N(I,JM1,1)		&
												+L%N(IP1,JM1,1))
!			   DF = 0.5*(L%H(I,J)+L%H(IP1,J))
			   IF (HM.GE.1.0E-3) THEN
			      FF = CONST*FM*FM*SQRT(L%M(I,J,1)**2+XQQ**2)		&
									/HM**(2.3333333)
			   ELSE
			      FF = 0.0
			   ENDIF
            ENDIF !END OF APPLYING BOTTOM FRICTION

            ALPHA = L%ALPHA(I,J)
			SHIFT = 0.5*(ALPHA+1.0)
            NN = FLOOR(SHIFT)

			DXS = L%DX*ALPHA
			DYS = L%DY*ALPHA

            IF (NL.EQ.1) THEN
			   BETA = SHIFT - NN
			   B1 = (1.0-BETA)/DXS
			   B2 = BETA/DXS
			ENDIF

            IPN = I+NN
			IMN = I-NN

            NS = NINT(A)
			IF (NS.LT.1) NS = 1
			JPNN = J+NS
			JMNN = J-NS

			IP1PN = IP1+NN
			IM1PN = IM1+NN
			IP1MN = IP1-NN
			IM1MN = IM1-NN
			IP2PN = IP2+NN
			IP2MN = IP2-NN

			ISPANS = MAX(NN+5,10)
!*            ISPANS = 10
			ISPANE = L%NX-ISPANS
			JSPANS = ISPANS
			JSPANE = L%NY-JSPANS

		    IF (L%MASK(I,J) .EQ. 1) THEN
			   !CALCULATE DISPERSION CORRECTION
			   IF (INTP.EQ.0) THEN   !FIXED POINT MULTIPLIED BY COEF.
				  ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)					&
						+ L%CF(I,J,2)*L%Z(IPN,J,2)					&
						+ L%CF(I,J,3)*L%Z(IP1PN,J,2)				&
						+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				  ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)					&
						+ L%CB(I,J,2)*L%Z(IMN,J,2)					&
						+ L%CB(I,J,3)*L%Z(IP1MN,J,2)				&
						+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
                  DZMDX = (ZMF-ZMB)
                  XM = (1.-FF)*L%M(I,J,1)-GRX/ALPHA*HM*DZMDX
				  IF (L%MODSCM .EQ. 0) THEN
					 !TERM TO CORRECT DIAGNAL DIFFERENTIATION
                     ZADD = ((L%Z(IP1,JP1,2)						&
								-2*L%Z(IP1,J,2)+L%Z(IP1,JM1,2))		&
								-(L%Z(I,JP1,2)-2*L%Z(I,J,2)			&
								+L%Z(I,JM1,2))) 
			         GAMMA = L%ALPHA(I,J)**2
					 SCM = GAMMA*GRX*TWLVTH*HM*ZADD
 			         XM = XM - SCM
				  ENDIF
			   ENDIF
			   IF (INTP.EQ.1) THEN  !YOON (2002) INTERPOLATION
		          A1X = L%A1X(I,J) !(ALPHA**2-1)/(24.*L%DX)
			      A2X = L%A2X(I,J) !3.*(9.-ALPHA**2)/(24.*L%DX)
			      DZMDX = A1X*(L%Z(I+2,J,2)-L%Z(I-1,J,2))			&
									+ A2X*(L%Z(I+1,J,2)-L%Z(I,J,2))
				  XM = (1.-FF)*L%M(I,J,1)-GRT*HM*DZMDX
				  IF (L%MODSCM .EQ. 0) THEN
			         IF (ALPHA .LT. 1.732) THEN
			            GAMMA = ALPHA**2
			            DZFDX = A1X*(L%Z(I+2,J+1,2)-L%Z(I-1,J+1,2))	&
								+A2X*(L%Z(I+1,J+1,2)-L%Z(I,J+1,2))
				        DZBDX = A1X*(L%Z(I+2,J-1,2)-L%Z(I-1,J-1,2))	&
								+A2X*(L%Z(I+1,J-1,2)-L%Z(I,J-1,2))
			         ELSE
			            GAMMA = ALPHA**2/4.0
			            DZFDX = A1X*(L%Z(I+2,J+2,2)-L%Z(I-1,J+2,2))	&
								+A2X*(L%Z(I+1,J+2,2)-L%Z(I,J+2,2))
				        DZBDX = A1X*(L%Z(I+2,J-2,2)-L%Z(I-1,J-2,2))	&
								+A2X*(L%Z(I+1,J-2,2)-L%Z(I,J-2,2))				     
				     ENDIF
			         SCM = GAMMA*GRT*HM*TWLVTH*(DZFDX-2*DZDX+DZBDX)
					 XM = XM - SCM
				  ENDIF
               ENDIF
			   IF (INTP.EQ.2) THEN !4-PTS LAGRANGE INTERPOLATION
                  ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)					&
						+ L%CF(I,J,2)*L%Z(IPN,J,2)					&
						+ L%CF(I,J,3)*L%Z(IP1PN,J,2)				&
						+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				  ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)					&
						+ L%CB(I,J,2)*L%Z(IMN,J,2)					&
						+ L%CB(I,J,3)*L%Z(IP1MN,J,2)				&
						+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
				  DZMDX = (ZMF-ZMB)
				  XM = (1.-FF)*L%M(I,J,1)-GRX/ALPHA*HM*DZMDX
				  IF (L%MODSCM .EQ. 0) THEN
                     ZUF = L%CF(I,J,1)*L%Z(IM1PN,JPNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JPNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JPNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JPNN,2)
				     ZUB = L%CB(I,J,1)*L%Z(IM1MN,JPNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JPNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JPNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JPNN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IM1PN,JMNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JMNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JMNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JMNN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IM1MN,JMNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JMNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JMNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JMNN,2)
				     DZUDX = (ZUF-ZUB)
				     DZLDX = (ZLF-ZLB)
				     ZADD = DZUDX-2.0*DZMDX+DZLDX
					 GAMMA = L%ALPHA(I,J)**2/(NS**2)
					 SCM = GAMMA*GRX/ALPHA*TWLVTH*HM*ZADD
 			         XM = XM - SCM
				  ENDIF
			   ENDIF
			   IF (INTP.EQ.3) THEN !4-PTS LAGRANGE INTERPOLATION TO O(DX^4)
                  ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)					&
						+ L%CF(I,J,2)*L%Z(IPN,J,2)					&
						+ L%CF(I,J,3)*L%Z(IP1PN,J,2)				&
						+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				  ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)					&
						+ L%CB(I,J,2)*L%Z(IMN,J,2)					&
						+ L%CB(I,J,3)*L%Z(IP1MN,J,2)				&
						+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
				  DZMDX = (ZMF-ZMB)
				  XM = (1.-FF)*L%M(I,J,1)-GRX/ALPHA*HM*DZMDX
				  IF (L%MODSCM .EQ. 0) THEN
                     ZUF = L%CF(I,J,1)*L%Z(IM1PN,JPNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JPNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JPNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JPNN,2)
				     ZUB = L%CB(I,J,1)*L%Z(IM1MN,JPNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JPNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JPNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JPNN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IM1PN,JMNN,2)			&
							+ L%CF(I,J,2)*L%Z(IPN,JMNN,2)			&
							+ L%CF(I,J,3)*L%Z(IP1PN,JMNN,2)			&
							+ L%CF(I,J,4)*L%Z(IP2PN,JMNN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IM1MN,JMNN,2)			&
							+ L%CB(I,J,2)*L%Z(IMN,JMNN,2)			&
							+ L%CB(I,J,3)*L%Z(IP1MN,JMNN,2)			&
							+ L%CB(I,J,4)*L%Z(IP2MN,JMNN,2)
                     ZUUF = L%CF(I,J,1)*L%Z(IM1PN,JPNN+NS,2)		&
							+ L%CF(I,J,2)*L%Z(IPN,JPNN+NS,2)		&
							+ L%CF(I,J,3)*L%Z(IP1PN,JPNN+NS,2)		&
							+ L%CF(I,J,4)*L%Z(IP2PN,JPNN+NS,2)
				     ZUUB = L%CB(I,J,1)*L%Z(IM1MN,JPNN+NS,2)		&
							+ L%CB(I,J,2)*L%Z(IMN,JPNN+NS,2)		&
							+ L%CB(I,J,3)*L%Z(IP1MN,JPNN+NS,2)		&
							+ L%CB(I,J,4)*L%Z(IP2MN,JPNN+NS,2)
                     ZLLF = L%CF(I,J,1)*L%Z(IM1PN,JMNN-NS,2)		&
							+ L%CF(I,J,2)*L%Z(IPN,JMNN-NS,2)		&
							+ L%CF(I,J,3)*L%Z(IP1PN,JMNN-NS,2)		&
							+ L%CF(I,J,4)*L%Z(IP2PN,JMNN-NS,2)
				     ZLLB = L%CB(I,J,1)*L%Z(IM1MN,JMNN-NS,2)		&
							+ L%CB(I,J,2)*L%Z(IMN,JMNN-NS,2)		&
							+ L%CB(I,J,3)*L%Z(IP1MN,JMNN-NS,2)		&
							+ L%CB(I,J,4)*L%Z(IP2MN,JMNN-NS,2)
				     DZUDX = (ZUF-ZUB)
				     DZLDX = (ZLF-ZLB)
					 DZUUDX = (ZUUF-ZUUB)
					 DZLLDX = (ZLLF-ZLLB)
				     ZADD = (-2.*DZUUDX+32.*DZUDX-60.0*DZMDX		&
										+32.*DZLDX-2.*DZLLDX)/24.
					 GAMMA = L%ALPHA(I,J)**2/(NS**2)
					 SCM = GAMMA*GRX/ALPHA*TWLVTH*HM*ZADD
 			         XM = XM - SCM
				  ENDIF
			   ENDIF
               !CALCULATE NONLINEAR TERMS AND NONLINEAR CORRECTIONS>>---
			   DUDX = 0.0
			   DWYDY = 0.0
			   VADD = 0.0
			   WADD = 0.0
			   FADD = 0.0
               IF (L%LAYGOV.EQ.3 .AND. L%MASK(I,J).EQ.1) THEN
				  IF (NL.EQ.0) THEN
				     DUDX = 0.0
					 DWYDY = 0.0
					 UADD = 0.0
					 WADD = 0.0 
				  ENDIF
			      IF (NL.EQ.1) THEN
			         UF = B1*U(IM1PN,J) + B2*U(IPN,J)
				     UB = B2*U(IM1MN,J) + B1*U(IMN,J)

				     WYF = B1*WY(I,JM1PN) + B2*WY(I,JPN)
				     WYB = B2*WY(I,JM1MN) + B1*WY(I,JMN)

                     UN = (B1*U(IPN,J)+B2*U(IP1PN,J))-(B2*U(IMN,J)	&
													+B1*U(IP1MN,J))
					 UN0 = (B1*U0(IM1PN,J)+B2*U0(IPN,J))			&
									-(B2*U0(IM1MN,J)+B1*U0(IMN,J))

                     WU = 0.5*((B1*W(IM1PN,J) + B2*W(IPN,J))		&
									+ (B1*W(I,JM1PN)+B2*W(I,JPN)))
				     WL = 0.5*((B1*W(IM1PN,JM1)+B2*W(IPN,JM1))		&
									+ (B2*W(I,JM1MN)+B1*W(I,JMN)))
				     W0U = 0.5*((B1*W0(IM1,JM1PN)+B2*W0(IM1,JPN))	&
								+ (B2*W0(IM1MN,J) + B1*W0(IMN,J)))
				     W0L = 0.5*((B2*W0(IM1,JM1MN)+B1*W0(IM1,JMN))	&
								+(B2*W0(IM1MN,JM1)+B1*W0(IMN,JM1)))

				     DUDX = (UF - UB)/DXS
				     DWYDY = (WYF - WYB)/DYS
				     UADD = (UN-UN0)/2.0/DXS
				     WADD = (WU-WL-W0U+W0L)/2.0/DYS
				  ENDIF
				  IF (NL.EQ.2) THEN
				     DUDX = (U(I,J)-U(IM1,J))/L%DX
					 DWYDY = (WY(I,J)-WY(I,JM1))/L%DY
					 UADD = (U(IP1,J)-U(I,J)-U0(I,J)+U0(IM1,J))		&
														/(2.0*L%DX)
					 WADD = (W(IP1,J)-W(IP1,JM1)-W0(I,J)+W0(I,JM1))	&
														/(2.0*L%DY)
				  ENDIF
			      IF (NL.EQ.3) THEN
				     IF (L%MODSCM .EQ. 1) A = 1.0
				     DUDX = (U(I,J)-U(IM1,J))/L%DX
					 DWYDY = (WY(I,J)-WY(I,JM1))/L%DY
					 UADD = (A*(U(IP1,J)-U(I,J))					&
								-(A-1.)*(U(I,J)-U(IM1,J))			&
								-(U0(I,J)-U0(IM1,J)))/(2.0*L%DX)
					 WADD = (A*(WY(IP1,J)-WY(IP1,JM1))				&
								-(A-1.)*(WY(I,J)-WY(I,JM1))			&
								-(WY0(I,J)-WY0(I,JM1)))/(2.0*L%DY)
				  ENDIF
                  XM = XM -L%DT*(DUDX+DWYDY) - L%DT*(UADD+WADD)
               ENDIF
			   !END OF NONLINEAR CALCULATION<<---
			   !CALCULATE 5TH-ORDER DERIVATIVES
			   IF (IT4.EQ.1) THEN
                  DZDX4 = ((L%Z(IP2,J,2)+L%Z(IM2,J,2))				&
								-4.*(L%Z(IP1,J,2)+L%Z(IM1,J,2))		&
								+6.*L%Z(I,J,2))/(L%DX**4)
			      DZDX3F = ((L%Z(IP2,J,2)-L%Z(IM2,J,2))				&
								-2.*(L%Z(IP1,J,2)-L%Z(IM1,J,2)))	&
								/(2.*L%DX**3)
			      DZDX3B = ((L%Z(IP2,J,1)-L%Z(IM2,J,1))				&
								-2.*(L%Z(IP1,J,1)-L%Z(IM1,J,1)))	&
								/(2.*L%DX**3)
				  DZDX3T = (DZDX3F - DZDX3B)/L%DT
				  DZDX2Y2 = (   (-2.*L%Z(IM2,JM1,2)					&
							+32.*L%Z(IM1,JM1,2)-60.*L%Z(I,JM1,2)	&
							+32.*L%Z(IP1,JM1,2)-2.*L%Z(IP2,JM1,2))	&
							-2.*(-2.*L%Z(IM2,J,2)+32.*L%Z(IM1,J,2)	&
							-60.*L%Z(I,J,2)+32.*L%Z(IP1,J,2)		&
							-2.*L%Z(IP2,J,2))+(-2.*L%Z(IM2,JP1,2)	&
							+32.*L%Z(IM1,JP1,2)-60.*L%Z(I,JP1,2)	&
							+32.*L%Z(IP1,JP1,2)-2.*L%Z(IP2,JP1,2))) &
							/(24.*L%DX**2*L%DY**2)

				  DZDXY2F = ((-2.*L%Z(IP1,JM2,2)+32.*L%Z(IP1,JM1,2)	&
							-60.*L%Z(IP1,J,2)+32.*L%Z(IP1,JP1,2)	&
							-2.*L%Z(IP1,JP2,2))-(-2.*L%Z(I,JM2,2)	&
							+32.*L%Z(I,JM1,2)-60.*L%Z(I,J,2)		&
							+32.*L%Z(I,JP1,2)-2.*L%Z(I,JP2,2)))		&
							/(24.*L%DX*L%DY**2)
				  DZDXY2B = ((-2.*L%Z(IP1,JM2,1)+32.*L%Z(IP1,JM1,1)	&
							-60.*L%Z(IP1,J,1)+32.*L%Z(IP1,JP1,1)	&
							-2.*L%Z(IP1,JP2,1))-(-2.*L%Z(I,JM2,1)	&
							+32.*L%Z(I,JM1,1)-60.*L%Z(I,J,1)		&
							+32.*L%Z(I,JP1,1)-2.*L%Z(I,JP2,1)))		&
							/(24.*L%DX*L%DY**2)
				  DZDXY2T = (DZDXY2F-DZDXY2B)/L%DT

				  TERM1 = (L%ALPHA(I,J)*L%DX)**3*DZDX4
				  TERM2 = (L%ALPHA(I,J)*L%DX)**2*L%DT*DZDX3T
				  TERM3 = -(L%ALPHA(I,J)*L%DX)*(L%DT**2)*(DZDX4+DZDX2Y2)
				  TERM4 = -(L%DT**3)*(DZDXY2T+DZDX3T)
				  FADD = -(HM*(TERM1+TERM2)+HM**2*(TERM3+TERM4))/48.

				  XM = XM - L%DT*FADD
			   ENDIF
			   !END OF CALCULATING 5TH-ORDER DERIVATIVES
			ELSE
			   XM = (1.-FF)*L%M(I,J,1)-GRX*HM*(L%Z(IP1,J,2)-L%Z(I,J,2))
			   IF (L%MODSCM .EQ. 0) THEN
			      SCM = GRX*TWLVTH*HM*((L%Z(IP1,JP1,2)			&
								-2*L%Z(IP1,J,2)+L%Z(IP1,JM1,2))		&
								-(L%Z(I,JP1,2)-2*L%Z(I,J,2)			&
								+L%Z(I,JM1,2)))
                  XM = XM - SCM
			   ENDIF
            ENDIF

            !USE BOTTOM FRICTION
			IF (L%FRIC_SWITCH.NE.1) XM = XM/(1.+FF)

            IF (ABS(XM) .LT. EPS) XM = ZERO
            L%M(I,J,2) = XM
          END IF
        END DO
      END DO
!
      DO J=JS,JYM1
        JP1 = J+1
		JP2 = J+2
		JM1 = J-1
		JM2 = J-2
		IF (JP1 .GE. L%NY) JP1 = L%NY
		IF (JP2 .GE. L%NY) JP2 = L%NY
		IF (JM1 .LE. 1) JM1 = 1
		IF (JM2 .LE. 1) JM2 = 1
        DO I=IS,L%NX
          IF ((L%H(I,J).GT.GX) .AND. (L%H(I,JP1).GT.GX)) THEN
            IM1 = I-1
			IM2 = I-2
			IP1 = I+1
			IP2 = I+2
			IF (IM1 .LT. 1) IM1 = 1
			IF (IM2 .LT. 1) IM2 = 1
            IF (IP1 .GE. L%NX) IP1 = L%NX
			IF (IP2 .GE. L%NX) IP2 = L%NX
			
			!WATER DEPTH AT DISCHARGE LOCATION P
			HN = 0.5*(L%H(I,J)+L%H(I,JP1)+L%Z(I,J,2)+L%Z(I,JP1,2)) 
            !USE BOTTOM FRICTION
			IF (L%FRIC_SWITCH.NE.1) THEN
			   IF (L%FRIC_SWITCH .EQ. 2) FM = L%FRIC_VCOEF(I,J)
			   XPP = 0.25*(L%M(I,J,1)+L%M(I,JP1,1)+L%M(IM1,J,1)		&
												+L%M(IM1,JP1,1))
!			   DF = 0.5*(L%H(I,J)+L%H(I,JP1))
			   IF (HN.GE.1.0E-5) THEN
			      FF = CONST*FM*FM*SQRT(L%N(I,J,1)**2+XPP**2)		&
									/HN**(2.3333333)
			   ELSE
			      FF = 0.0
			   ENDIF
            ENDIF !END OF APPLYING BOTTOM FRICTION

            ALPHA = L%ALPHA(I,J)
			SHIFT = 0.5*(ALPHA+1.0)
            NN = FLOOR(SHIFT)

			DXS = L%DX*L%ALPHA(I,J)
			DYS = L%DY*L%ALPHA(I,J)

            IF (NL.EQ.1) THEN
			   BETA = SHIFT - NN
			   B1 = (1.0-BETA)/DYS
			   B2 = BETA/DYS
			ENDIF

            JPN = J+NN
			JMN = J-NN
			JP1PN = JP1+NN
			JM1PN = JM1+NN
			JP1MN = JP1-NN
			JM1MN = JM1-NN
			JP2PN = JP2+NN
			JP2MN = JP2-NN

            NS = NINT(A)
			IF (NS.LT.1) NS = 1
		    IPNN = I+NS
			IMNN = I-NS


            GAMMA = L%ALPHA(I,J)**2
			SCM = 0.0

		    IF (L%MASK(I,J) .EQ. 1) THEN
               !CALCULATE DISPERSION CORRECTIONS
			   IF (INTP.EQ.0) THEN
				  ZMF = L%CF(I,J,1)*L%Z(I,JM1PN,2)					&
						+ L%CF(I,J,2)*L%Z(I,JPN,2)					&
						+ L%CF(I,J,3)*L%Z(I,JP1PN,2)				&
						+ L%CF(I,J,4)*L%Z(I,JP2PN,2)
				  ZMB = L%CB(I,J,1)*L%Z(I,JM1MN,2)					&
						+ L%CB(I,J,2)*L%Z(I,JMN,2)					&
						+ L%CB(I,J,3)*L%Z(I,JP1MN,2)				&
						+ L%CB(I,J,4)*L%Z(I,JP2MN,2)
                  DZMDY = (ZMF-ZMB)
                  XN = (1.-FF)*L%N(I,J,1)-GRY/ALPHA*HN*DZMDY
				  IF (L%MODSCM.EQ.0) THEN
				     ZADD = 1./L%DY*((L%Z(IP1,JP1,2)-2*L%Z(I,JP1,2)	&
									+L%Z(IM1,JP1,2))-(L%Z(IP1,J,2)	&
									-2*L%Z(I,J,2)+L%Z(IM1,J,2)))
					 GAMMA = L%ALPHA(I,J)**2
					 SCM = GAMMA*GRY*TWLVTH*HN*ZADD
 			         XN = XN - SCM
				  ENDIF
			   ENDIF
			   IF (INTP.EQ.1) THEN
		          A1Y = L%A1Y(I,J) !(ALPHA**2-1)/(24.*L%DY)
			      A2Y = L%A2Y(I,J) !3.*(9.-ALPHA**2)/(24.*L%DY)
			      DZDY = A1Y*(L%Z(I,JP2,2)-L%Z(I,J-1,2))			&
						+A2Y*(L%Z(I,JP1,2)-L%Z(I,J,2))
                  XN = (1.-FF)*L%N(I,J,1)-GRT*HN*DZDY
			      IF (L%MODSCM .EQ. 0) THEN
			         IF (ALPHA .LT. 1.732) THEN
			            GAMMA = ALPHA**2
			            DZUDY = A1Y*(L%Z(IP1,JP2,2)-L%Z(IP1,JM1,2))	&
								+A2Y*(L%Z(IP1,JP1,2)-L%Z(IP1,J,2))
				        DZLDY = A1Y*(L%Z(IM1,JM2,2)-L%Z(IM1,JM1,2))	&
								+A2Y*(L%Z(IM1,JP1,2)-L%Z(IM1,J,2))
			         ELSE
			            GAMMA = ALPHA**2/4.0
			            DZUDY = A1Y*(L%Z(IP2,JP2,2)-L%Z(IP2,JM1,2))	&
								+A2Y*(L%Z(IP2,JP1,2)-L%Z(IP2,J,2))
				        DZLDY = A1Y*(L%Z(IM2,JP2,2)-L%Z(IM2,JM1,2))	&
								+A2Y*(L%Z(IM2,JP1,2)-L%Z(IM2,J,2))				     
				     ENDIF
			         SCM = GAMMA*GRT*HN*TWLVTH					&
											*(DZUDY-2*DZDY+DZLDY)
					 XN = XN - SCM
                  ENDIF
			   ENDIF
			   IF (INTP.EQ.2) THEN
                  ZMF = L%CF(I,J,1)*L%Z(I,JM1PN,2)					&
						+ L%CF(I,J,2)*L%Z(I,JPN,2)					&
						+ L%CF(I,J,3)*L%Z(I,JP1PN,2)				&
						+ L%CF(I,J,4)*L%Z(I,JP2PN,2)
				  ZMB = L%CB(I,J,1)*L%Z(I,JM1MN,2)					&
						+ L%CB(I,J,2)*L%Z(I,JMN,2)					&
						+ L%CB(I,J,3)*L%Z(I,JP1MN,2)				&
						+ L%CB(I,J,4)*L%Z(I,JP2MN,2)
				  DZMDY = (ZMF-ZMB)
				  XN = (1.-FF)*L%N(I,J,1)-GRY/ALPHA*HN*DZMDY
				  IF (L%MODSCM.EQ.0) THEN
                     ZRF = L%CF(I,J,1)*L%Z(IPNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IPNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IPNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IPNN,JP2PN,2)
				     ZRB = L%CB(I,J,1)*L%Z(IPNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IPNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IPNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IPNN,JP2MN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IMNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IMNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IMNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IMNN,JP2PN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IMNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IMNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IMNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IMNN,JP2MN,2)
                     
				     DZRDY = (ZRF-ZRB)
				     DZLDY = (ZLF-ZLB)
				     ZADD = DZRDY-2.0*DZMDY+DZLDY
					 GAMMA = L%ALPHA(I,J)**2/NS**2
					 SCM = GAMMA*GRY*TWLVTH*HN*ZADD
 			         XN = XN - SCM
				  ENDIF
			   ENDIF
			   IF (INTP.EQ.3) THEN
                  ZMF = L%CF(I,J,1)*L%Z(I,JM1PN,2)					&
							+ L%CF(I,J,2)*L%Z(I,JPN,2)				&
							+ L%CF(I,J,3)*L%Z(I,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(I,JP2PN,2)
				  ZMB = L%CB(I,J,1)*L%Z(I,JM1MN,2)					&
							+ L%CB(I,J,2)*L%Z(I,JMN,2)				&
							+ L%CB(I,J,3)*L%Z(I,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(I,JP2MN,2)
				  DZMDY = (ZMF-ZMB)
				  XN = (1.-FF)*L%N(I,J,1)-GRY/ALPHA*HN*DZMDY
				  IF (L%MODSCM.EQ.0) THEN
                     ZRF = L%CF(I,J,1)*L%Z(IPNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IPNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IPNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IPNN,JP2PN,2)
				     ZRB = L%CB(I,J,1)*L%Z(IPNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IPNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IPNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IPNN,JP2MN,2)
                     ZRRF = L%CF(I,J,1)*L%Z(IPNN+NS,JM1PN,2)		&
							+ L%CF(I,J,2)*L%Z(IPNN+NS,JPN,2)		&
							+ L%CF(I,J,3)*L%Z(IPNN+NS,JP1PN,2)		&
							+ L%CF(I,J,4)*L%Z(IPNN+NS,JP2PN,2)
				     ZRRB = L%CB(I,J,1)*L%Z(IPNN+NS,JM1MN,2)		&
							+ L%CB(I,J,2)*L%Z(IPNN+NS,JMN,2)		&
							+ L%CB(I,J,3)*L%Z(IPNN+NS,JP1MN,2)		&
							+ L%CB(I,J,4)*L%Z(IPNN+NS,JP2MN,2)
                     ZLF = L%CF(I,J,1)*L%Z(IMNN,JM1PN,2)			&
							+ L%CF(I,J,2)*L%Z(IMNN,JPN,2)			&
							+ L%CF(I,J,3)*L%Z(IMNN,JP1PN,2)			&
							+ L%CF(I,J,4)*L%Z(IMNN,JP2PN,2)
				     ZLB = L%CB(I,J,1)*L%Z(IMNN,JM1MN,2)			&
							+ L%CB(I,J,2)*L%Z(IMNN,JMN,2)			&
							+ L%CB(I,J,3)*L%Z(IMNN,JP1MN,2)			&
							+ L%CB(I,J,4)*L%Z(IMNN,JP2MN,2)
                     ZLLF = L%CF(I,J,1)*L%Z(IMNN-NS,JM1PN,2)		&
							+ L%CF(I,J,2)*L%Z(IMNN-NS,JPN,2)		&
							+ L%CF(I,J,3)*L%Z(IMNN-NS,JP1PN,2)		&
							+ L%CF(I,J,4)*L%Z(IMNN-NS,JP2PN,2)
				     ZLLB = L%CB(I,J,1)*L%Z(IMNN-NS,JM1MN,2)		&
							+ L%CB(I,J,2)*L%Z(IMNN-NS,JMN,2)		&
							+ L%CB(I,J,3)*L%Z(IMNN-NS,JP1MN,2)		&
							+ L%CB(I,J,4)*L%Z(IMNN-NS,JP2MN,2)
                     
				     DZRDY = (ZRF-ZRB)
				     DZLDY = (ZLF-ZLB)
				     DZRRDY = (ZRRF-ZRRB)
				     DZLLDY = (ZLLF-ZLLB)
				     ZADD = (-2.*DZRRDY+32.*DZRDY-60.*DZMDY			&
										+32.*DZLDY-2.*DZLLDY)/24.
					 GAMMA = L%ALPHA(I,J)**2/NS**2
					 SCM = GAMMA*GRY*TWLVTH*HN*ZADD
 			         XN = XN - SCM
				  ENDIF
			   ENDIF
               
			   !START TO CALCULATING NONLINEAR TERMS
			   DVDY = 0.0
			   DWXDX = 0.0
			   VADD = 0.0
			   WADD = 0.0
			   FADD = 0.0
			   IF (L%LAYGOV.EQ.3 .AND. L%MASK(I,J).EQ.1) THEN
			      IF (NL.EQ.1) THEN
			         VF = B1*V(I,JM1PN) + B2*V(I,JPN)
				     VB = B2*V(I,JM1MN) + B1*V(I,JMN)

				     WXF = B1*WX(IM1PN,J) + B2*WX(IPN,J)
				     WXB = B2*WX(IM1MN,J) + B1*WX(IMN,J)

                     VN  = (B1*V(I,JPN)+B2*V(I,JP1PN))				&
									- (B2*V(I,JMN)+B1*V(I,JP1MN))
					 VN0 = (B1*V0(I,JM1PN)+B2*V0(I,JPN))			&
									- (B2*V0(I,JM1MN)+B1*V0(I,JMN))

                     WR  = 0.5*((B1*W(I,JM1PN) + B2*W(I,JPN))		&
									+ (B1*W(IM1PN,J)+B2*W(IPN,J)))
				     WL  = 0.5*((B1*W(IM1,JM1PN)+B2*W(IM1,JPN))		&
									+ (B2*W(IM1MN,J)+B1*W(IMN,J)))
				     W0R = 0.5*((B1*W0(IM1PN,JM1)+B2*W0(IPN,JM1))	&
								+ (B2*W0(IM1MN,J) + B1*W0(IMN,J)))
				     W0L = 0.5*((B2*W0(IM1MN,JM1)+B1*W0(IMN,JM1))	&
								+ (B2*W0(IM1,JM1MN)+B1*W0(IM1,JMN)))

				     DVDY = (VF - VB)/DYS
				     DWXDX = (WXF - WXB)/DXS
				     VADD = (VN-VN0)/2.0/DXS
				     WADD = (WR-WL-W0R+W0L)/2.0/DXS
				  ENDIF
				  IF (NL.EQ.2) THEN
				     DVDY = (V(I,J)-V(I,JM1))/L%DY
					 DWXDX = (WX(I,J)-WX(IM1,J))/L%DX
					 VADD = (V(I,JP1)-V(I,J)-V0(I,J)+V0(I,JM1))		&
														/(2.0*L%DY)
					 WADD = (W(I,JP1)-W(IM1,JP1)-W0(I,J)+W0(IM1,J))	&
														/(2.0*L%DX)
				  ENDIF
				  IF (NL.EQ.3) THEN
				     IF (L%MODSCM .EQ. 1) A = 1.0
				     DVDY = (V(I,J)-V(I,JM1))/L%DY
					 DWXDX = (WX(I,J)-WX(IM1,J))/L%DX
					 VADD = (A*(V(I,JP1)-V(I,J))					&
								-(A-1.)*(V(I,J)-V(I,JM1))			&
								-(V0(I,J)-V0(I,JM1)))/(2.0*L%DY)
					 WADD = (A*(WX(I,JP1)-WX(IM1,JP1))				&
								-(A-1.)*(WX(I,J)-WX(IM1,J))			&
								-(WX0(I,J)-WX0(IM1,J)))/(2.0*L%DX)
				  ENDIF
				  IF (NL.EQ.0) THEN
				     DVDX = 0.0
					 DWXDX = 0.0
					 VADD = 0.0
					 WADD = 0.0 
				  ENDIF
                  XN = XN - L%DT*(DVDY+DWXDX) -L%DT*(VADD+WADD)
               ENDIF
			   !END OF CALCULATING NONLINEAR TERMS
			   !CALCULATING 5TH-ORDER DERIVATIVES
			   IF (IT4.EQ.1) THEN
                     DZDY4 = ((L%Z(I,JP2,2)+L%Z(I,JM2,2))			&
								-4.*(L%Z(I,JP1,2)+L%Z(I,JM1,2))		&
								+6.*L%Z(I,J,2))/(L%DY**4)
					 DZDY3F = ((L%Z(I,JP2,2)-L%Z(I,JM2,2))			&
								-2.*(L%Z(I,JP1,2)-L%Z(I,JM1,2)))	&
								/(2.*L%DY**3)
					 DZDY3B = ((L%Z(I,JP2,1)-L%Z(I,JM2,1))			&
								-2.*(L%Z(I,JP1,1)-L%Z(I,JM1,1)))	&
								/(2.*L%DY**3)
					 DZDY3T = (DZDY3F - DZDY3B)/L%DT
					 DZDX2Y2 = (   (-2.*L%Z(IM1,JM2,2)				&
							+32.*L%Z(IM1,JM1,2)-60.*L%Z(IM1,J,2)	&
							+32.*L%Z(IM1,JP1,2)-2.*L%Z(IM1,JP2,2))	&
							-2.*(-2.*L%Z(I,JM2,2)+32.*L%Z(I,JM1,2)	&
							-60.*L%Z(I,J,2)+32.*L%Z(I,JP1,2)		&
							-2.*L%Z(I,JP2,2))+(-2.*L%Z(IP1,JM2,2)	&
							+32.*L%Z(IP1,JM1,2)-60.*L%Z(IP1,J,2)	&
							+32.*L%Z(IP1,JP1,2)-2.*L%Z(IP1,JP2,2))) &
							/(24.*L%DX**2*L%DY**2)

					 DZDX2YF = ((-2.*L%Z(IM2,JP1,2)					&
							+32.*L%Z(IM1,JP1,2)-60.*L%Z(I,JP1,2)	&
							+32.*L%Z(IP1,JP1,2)-2.*L%Z(IP2,JP1,2))	&
					        -(-2.*L%Z(IM2,J,2)+32.*L%Z(IM1,J,2)		&
							-60.*L%Z(I,J,2)+32.*L%Z(IP1,J,2)		&
							-2.*L%Z(IP2,J,2)))/(24.*L%DX**2*L%DY)
					 DZDX2YB = ((-2.*L%Z(IM2,JP1,1)					&
							+32.*L%Z(IM1,JP1,1)-60.*L%Z(I,JP1,1)	&
							+32.*L%Z(IP1,JP1,1)-2.*L%Z(IP2,JP1,1))	&
					        -(-2.*L%Z(IM2,J,1)+32.*L%Z(IM1,J,1)		&
							-60.*L%Z(I,J,1)+32.*L%Z(IP1,J,1)		&
							-2.*L%Z(IP2,J,1)))/(24.*L%DX**2*L%DY)
					 DZDX2YT = (DZDX2YF-DZDX2YB)/L%DT

					 TERM1 = (L%ALPHA(I,J)*L%DY)**3*DZDY4
					 TERM2 = (L%ALPHA(I,J)*L%DY)**2*L%DT*DZDY3T
					 TERM3 = -(L%ALPHA(I,J)*L%DY)*(L%DT**2)			&
													*(DZDY4+DZDX2Y2)
					 TERM4 = -(L%DT**3)*(DZDX2YT+DZDY3T)
					 FADD = -(HN*(TERM1+TERM2)+HN**2*(TERM3+TERM4))/48.

					 XN = XN - L%DT*FADD
			   ENDIF
			   !END OF CALCULATING 5TH-ORDER DERIVATIVES                    
			ELSE
			   XN = (1.-FF)*L%N(I,J,1)-GRY*HN*(L%Z(I,JP1,2)-L%Z(I,J,2))
			   IF (L%MODSCM .EQ. 0) THEN
			      SCM = GRY*TWLVTH*HN*((L%Z(IP1,JP1,2)			&
							-2*L%Z(I,JP1,2)+L%Z(IM1,JP1,2))			&
							-(L%Z(IP1,J,2)-2*L%Z(I,J,2)+L%Z(IM1,J,2)))
                  XN = XN - SCM
			   ENDIF
            ENDIF

!           !USE BOTTOM FRICTION            
            IF (L%FRIC_SWITCH.NE.1) XN = XN/(1.0+FF)
            IF (ABS(XN) .LT. EPS) XN = ZERO
            L%N(I,J,2) = XN
          END IF
        END DO
      END DO
	  ENDIF

	  IF (IDIM.EQ.1) THEN
	     J = FLOOR(L%NY/2.0)
	     IXM1 = L%NX-1
	     IF (L%ID .EQ. 1) THEN
	       IS = 1
	     ELSE
	       IS = 2
	     END IF
         DO I=IS,IXM1
            IP1 = I+1
		    IP2 = I+2
		    IM1 = I-1
		    IM2 = I-2
		    IF (IP1 .GE. L%NX) IP1 = L%NX
		    IF (IP2 .GE. L%NX) IP2 = L%NX
		    IF (IM1 .LE. 1) IM1 = 1
		    IF (IM2 .LE. 1) IM2 = 1
            IF ((L%H(I,J).GT.GX) .AND. (L%H(IP1,J).GT.GX)) THEN
               HM = L%HP(I,J)+0.5*(L%Z(I,J,2)+L%Z(IP1,J,2))
               A = L%ALPHA(I,J)
			   SHIFT = 0.5*(A+1.0)
               NN = FLOOR(SHIFT)

			   DXS = L%DX*L%ALPHA(I,J)
               IF (INTP.EQ.1) THEN
			      BETA = SHIFT - NN
			      B1 = (1.0-BETA)/DXS
			      B2 = BETA/DXS
			   ENDIF

               IPN = I+NN
			   IMN = I-NN

               NS = NINT(A)
			   IF (NS.LT.1) NS = 1

			   IP1PN = IP1+NN
			   IM1PN = IM1+NN
			   IP1MN = IP1-NN
			   IM1MN = IM1-NN
			   IP2PN = IP2+NN
			   IP2MN = IP2-NN

			   SCM = 0.0
			   ZX1 = 0.0
			   ZX2 = 0.0
			   ZX3 = 0.0
			   Z1 = 0.0
			   Z2 = 0.0
			   ZADD = 0.0
			   IF (L%MASK(I,J) .EQ. 1) THEN
			      IF (INTP.EQ.1) THEN  !YOON (2002) INTERPOLATION
		             A1X = L%A1X(I,J) !(ALPHA**2-1)/(24.*L%DX)
			         A2X = L%A2X(I,J) !3.*(9.-ALPHA**2)/(24.*L%DX)
			         DZMDX = A1X*(L%Z(I+2,J,2)-L%Z(I-1,J,2))		&
									+ A2X*(L%Z(I+1,J,2)-L%Z(I,J,2))
				     XM = (1.-FF)*L%M(I,J,1)-GRT*HM*DZMDX
				     IF (L%MODSCM .EQ. 0) THEN
						SCM = 0.0
					    XM = XM - SCM
				     ENDIF
                  ENDIF

			      IF (INTP.EQ.2) THEN !4-PTS LAGRANGE INTERPOLATION
                     ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)				&
								+ L%CF(I,J,2)*L%Z(IPN,J,2)			&
								+ L%CF(I,J,3)*L%Z(IP1PN,J,2)		&
								+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				     ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)				&
								+ L%CB(I,J,2)*L%Z(IMN,J,2)			&
								+ L%CB(I,J,3)*L%Z(IP1MN,J,2)		&
								+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
				     DZMDX = (ZMF-ZMB)
				     XM = (1.-FF)*L%M(I,J,1)-GRX/ALPHA*HM*DZMDX
				     IF (L%MODSCM .EQ. 0) THEN
					    SCM = 0.0
 			            XM = XM - SCM
				     ENDIF
			      ENDIF
				 ! COMPUTE NONLINEAR CONVECTION TERMS				     
				  IF (NL.EQ.2) THEN
					 DUDX = (U(I,J)-U(IM1,J))/L%DX
                     UADD = ((U(IP1,J)-U(I,J))						&
									-(U0(I,J)-U0(IM1,J)))/(2.0*L%DX)
					 XM = XM -L%DT*(DUDX) - L%DT*(UADD)
			      ENDIF
               ELSE
				  DZDX = (L%Z(IP1,J,2)-L%Z(I,J,2))/L%DX
			      XM = L%M(I,J,1)-GRT*HM*DZDX
				  IF (L%MODSCM .EQ. 0) THEN
					 SCM = 0.0
 			         XM = XM - SCM
				  ENDIF	
				  DUDX = (U(I,J)-U(IM1,J))/L%DX
				  UADD = ((U(IP1,J)-U(I,J))							&
									-(U0(I,J)-U0(IM1,J)))/(2.0*L%DX)			  
                  XM = XM -L%DT*(DUDX) - L%DT*(UADD)
			   ENDIF
            ELSE
			   XM = 0.0
			ENDIF
			IF (ABS(XM) .LT. EPS) XM = ZERO
            L%M(I,:,2) = XM
	     ENDDO
	  ENDIF

!
      RETURN
      END




!----------------------------------------------------------------------
      SUBROUTINE DPDX_DISP (L,DPDX,DQDY,I,J)		 
! ....SOLVE CONTINUITY EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!
!NOTES: 
!	  #. ADD SUPPORT FOR DX\=DY (FOR HIGH-LATITUDE, 05/04/2007)
!            RX = L%R
!            RY = L%DT/L%DY

!     SUBROUTINE DESIGNED TO INCCLUDE DISPERSION IN SWE
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL ALPHA, DPDX, DQDY, A1, A2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/

!*	  RX = L%RX
!*	  RY = L%RY
!*	  GRT = L%DT*GRAV
	  
!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
!*      IDIM = 2   
!.....INTERPOLATION METHOD. 
!		1 - YOON2002; 
!		2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2   

!
!*      IS = 2
!*      JS = 2
!*	  IF(L%ID .EQ. 1)THEN  !OUTTEST LAYER				 
!*	    IE = L%NX -1
!*	    JE = L%NY -1
!*	  ELSE				  ! INNER LAYER
!*	    IE = L%NX
!*	    JE = L%NY
!*	  ENDIF

!*      IF (IDIM.EQ.2) THEN
!*      DO J=JS,JE
        JM1 = J-1
		JM2 = J-2
		JP1 = J+1
		JP2 = J+2
		IF (JM1.LE.1) JM1 = 1
		IF (JM2.LE.1) JM2 = 1
		IF (JP1.GE.L%NY) JP1 = L%NY
		IF (JP2.GE.L%NY) JP2 = L%NY
!*        DO I=IS,IE
		  IM1 = I-1
		  IM2 = I-2
		  IP1 = I+1
		  IP2 = I+2
		  IF (IM1.LE.1) IM1 = 1
		  IF (IM2.LE.1) IM2 = 1
		  IF (IP1.GE.L%NX) IP1 = L%NX
		  IF (IP2.GE.L%NX) IP2 = L%NX

		    IF (L%MASK(I,J).EQ.1) THEN
			   IF (INTP.EQ.1) THEN
		          A1X = L%A1X(I,J)   !(ALPHA**2-1)/(24.*L%DX)
			      A1Y = L%A1Y(I,J)
			      A2X = L%A2X(I,J)   !3.*(9.-ALPHA**2)/(24.*L%DX)
			      A2Y = L%A2Y(I,J)
			      DPDX = A1X*(L%M(IP1,J,1)-L%M(IM2,J,1))			&
							+A2X*(L%M(I,J,1)-L%M(IM1,J,1))
                  DQDY = A1Y*(L%N(I,JP1,1)-L%N(I,JM2,1))			&
							+A2Y*(L%N(I,J,1)-L%N(I,JM1,1))
			   ENDIF
			   IF  (INTP.EQ.2) THEN
				  A = L%ALPHA(I,J)
			      SHIFT = 0.5*(A+1.0)
                  NN = FLOOR(SHIFT)

		          DXS = L%DX*L%ALPHA(I,J)
			      DYS = L%DY*L%ALPHA(I,J)

                  PF = L%CF(I,J,1)*L%M(IM1+NN,J,1)					&
						+ L%CF(I,J,2)*L%M(I+NN,J,1)					&
						+ L%CF(I,J,3)*L%M(IP1+NN,J,1)				&
						+ L%CF(I,J,4)*L%M(IP2+NN,J,1)
                  PB = L%CB(I,J,1)*L%M(IM1-NN,J,1)					&
						+ L%CB(I,J,2)*L%M(I-NN,J,1)					&
						+ L%CB(I,J,3)*L%M(IP1-NN,J,1)				&
						+ L%CB(I,J,4)*L%M(IP2-NN,J,1)
				  QF = L%CF(I,J,1)*L%N(I,JM1+NN,1)					&
						+ L%CF(I,J,2)*L%N(I,J+NN,1)					&
						+ L%CF(I,J,3)*L%N(I,JP1+NN,1)				&
						+ L%CF(I,J,4)*L%N(I,JP2+NN,1)
				  QB = L%CB(I,J,1)*L%N(I,JM1-NN,1)					&
						+ L%CB(I,J,2)*L%N(I,J-NN,1)					&
						+ L%CB(I,J,3)*L%N(I,JP1-NN,1)				&
						+ L%CB(I,J,4)*L%N(I,JP2-NN,1)
				  DP = (PF-PB)
			      DQ = (QF-QB)
				  DPDX = DP/DXS
				  DQDY = DQ/DYS
			   ENDIF
			ENDIF
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE DZDX_DISP (L,DZDX,I,J)		 
! ....SOLVE CONTINUITY EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!
!NOTES: 
!	  #. ADD SUPPORT FOR DX\=DY (FOR HIGH-LATITUDE, 05/04/2007)
!            RX = L%R
!            RY = L%DT/L%DY

!     SUBROUTINE DESIGNED TO INCCLUDE DISPERSION IN SWE
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL ALPHA, DPDX, DQDY, A1, A2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/

!*	  RX = L%RX
!*	  RY = L%RY
!*	  GRT = L%DT*GRAV
	  
!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
!*      IDIM = 2   
!.....INTERPOLATION METHOD. 
!		1 - YOON2002; 
!		2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2   

!
!*      IS = 2
!*      JS = 2
!*	  IF(L%ID .EQ. 1)THEN  !OUTTEST LAYER				 
!*	    IE = L%NX -1
!*	    JE = L%NY -1
!*	  ELSE				  ! INNER LAYER
!*	    IE = L%NX
!*	    JE = L%NY
!*	  ENDIF

!*      IF (IDIM.EQ.2) THEN
!*      DO J=JS,JE
        JM1 = J-1
		JM2 = J-2
		JP1 = J+1
		JP2 = J+2
		IF (JM1.LE.1) JM1 = 1
		IF (JM2.LE.1) JM2 = 1
		IF (JP1.GE.L%NY) JP1 = L%NY
		IF (JP2.GE.L%NY) JP2 = L%NY
!*        DO I=IS,IE
		  IM1 = I-1
		  IM2 = I-2
		  IP1 = I+1
		  IP2 = I+2
		  IF (IM1.LE.1) IM1 = 1
		  IF (IM2.LE.1) IM2 = 1
		  IF (IP1.GE.L%NX) IP1 = L%NX
		  IF (IP2.GE.L%NX) IP2 = L%NX

!*          IF (L%H(I,J) .GT. GX) THEN



		    IF (L%MASK(I,J).EQ.1) THEN
			   IF (INTP.EQ.1) THEN
		          A1X = L%A1X(I,J) !(ALPHA**2-1)/(24.*L%DX)
			      A2X = L%A2X(I,J) !3.*(9.-ALPHA**2)/(24.*L%DX)
			      DZDX = A1X*(L%Z(I+2,J,2)-L%Z(I-1,J,2))			&
									+ A2X*(L%Z(I+1,J,2)-L%Z(I,J,2))
			   ENDIF
			   IF  (INTP.EQ.2) THEN
                  A = L%ALPHA(I,J)
			      SHIFT = 0.5*(A+1.0)
                  NN = FLOOR(SHIFT)

                  IPN = I+NN
			      IMN = I-NN

                  NS = NINT(A)
			      IF (NS.LT.1) NS = 1
			      JPNN = J+NS
			      JMNN = J-NS

			      IP1PN = IP1+NN
			      IM1PN = IM1+NN
			      IP1MN = IP1-NN
			      IM1MN = IM1-NN
			      IP2PN = IP2+NN
			      IP2MN = IP2-NN

		          DXS = L%DX*L%ALPHA(I,J)
                  ZMF = L%CF(I,J,1)*L%Z(IM1PN,J,2)					&
						+ L%CF(I,J,2)*L%Z(IPN,J,2)					&
						+ L%CF(I,J,3)*L%Z(IP1PN,J,2)				&
						+ L%CF(I,J,4)*L%Z(IP2PN,J,2)
				  ZMB = L%CB(I,J,1)*L%Z(IM1MN,J,2)					&
						+ L%CB(I,J,2)*L%Z(IMN,J,2)					&
						+ L%CB(I,J,3)*L%Z(IP1MN,J,2)				&
						+ L%CB(I,J,4)*L%Z(IP2MN,J,2)
				  DZDX = (ZMF-ZMB)/DXS
			   ENDIF
			ENDIF

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE DZDY_DISP (L,DZDY,I,J)		 
! ....SOLVE CONTINUITY EQUATION (LINEAR) IN CARTESIAN COORD. 
!     LAYID = 1, OUTEST LAYER
!     OTHERWISE, INNER LAYER 
!
!NOTES: 
!	  #. ADD SUPPORT FOR DX\=DY (FOR HIGH-LATITUDE, 05/04/2007)
!            RX = L%R
!            RY = L%DT/L%DY

!     SUBROUTINE DESIGNED TO INCCLUDE DISPERSION IN SWE
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: L
	  REAL ALPHA, DPDX, DQDY, A1, A2
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA UB/99.0/

!*	  RX = L%RX
!*	  RY = L%RY
!*	  GRT = L%DT*GRAV
	  
!.....DIMENSIONALITY INCLUDED IN THE COMPUTATION
!*      IDIM = 2   
!.....INTERPOLATION METHOD. 
!		1 - YOON2002; 
!		2 - LAGRANGE 4 POINTS (CUBIC)
      INTP = 2  

!
!*      IS = 2
!*      JS = 2
!*	  IF(L%ID .EQ. 1)THEN  !OUTTEST LAYER				 
!*	    IE = L%NX -1
!*	    JE = L%NY -1
!*	  ELSE				  ! INNER LAYER
!*	    IE = L%NX
!*	    JE = L%NY
!*	  ENDIF

!*      IF (IDIM.EQ.2) THEN
!*      DO J=JS,JE
        JM1 = J-1
		JM2 = J-2
		JP1 = J+1
		JP2 = J+2
		IF (JM1.LE.1) JM1 = 1
		IF (JM2.LE.1) JM2 = 1
		IF (JP1.GE.L%NY) JP1 = L%NY
		IF (JP2.GE.L%NY) JP2 = L%NY
!*        DO I=IS,IE
		  IM1 = I-1
		  IM2 = I-2
		  IP1 = I+1
		  IP2 = I+2
		  IF (IM1.LE.1) IM1 = 1
		  IF (IM2.LE.1) IM2 = 1
		  IF (IP1.GE.L%NX) IP1 = L%NX
		  IF (IP2.GE.L%NX) IP2 = L%NX

!*          IF (L%H(I,J) .GT. GX) THEN

		    IF (L%MASK(I,J).EQ.1) THEN
			   IF (INTP.EQ.1) THEN
		          A1Y = L%A1Y(I,J) !(ALPHA**2-1)/(24.*L%DY)
			      A2Y = L%A2Y(I,J) !3.*(9.-ALPHA**2)/(24.*L%DY)
			      DZDY = A1Y*(L%Z(I,JP2,2)-L%Z(I,J-1,2))			&
						+A2Y*(L%Z(I,JP1,2)-L%Z(I,J,2))
			   ENDIF
			   IF (INTP.EQ.2) THEN
                  A = L%ALPHA(I,J)
			      SHIFT = 0.5*(A+1.0)
                  NN = FLOOR(SHIFT)

                  JPN = J+NN
			      JMN = J-NN
			      JP1PN = JP1+NN
			      JM1PN = JM1+NN
			      JP1MN = JP1-NN
			      JM1MN = JM1-NN
			      JP2PN = JP2+NN
			      JP2MN = JP2-NN

                  NS = NINT(A)
			      IF (NS.LT.1) NS = 1
		          IPNN = I+NS
			      IMNN = I-NS

			      DYS = L%DY*L%ALPHA(I,J)
                  ZMF = L%CF(I,J,1)*L%Z(I,JM1PN,2)					&
						+ L%CF(I,J,2)*L%Z(I,JPN,2)					&
						+ L%CF(I,J,3)*L%Z(I,JP1PN,2)				&
						+ L%CF(I,J,4)*L%Z(I,JP2PN,2)
				  ZMB = L%CB(I,J,1)*L%Z(I,JM1MN,2)					&
						+ L%CB(I,J,2)*L%Z(I,JMN,2)					&
						+ L%CB(I,J,3)*L%Z(I,JP1MN,2)				&
						+ L%CB(I,J,4)*L%Z(I,JP2MN,2)
				  DZDY = (ZMF-ZMB)/DYS
			   ENDIF
			ENDIF

      RETURN
      END
