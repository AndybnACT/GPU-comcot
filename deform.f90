!----------------------------------------------------------------------
      SUBROUTINE GET_FLOOR_DEFORM (LO,LA,FLT,TIME)
!......................................................................
!DESCRIPTION:
!     #. CALCULATE SEAFLOOR DEFORMATION (OKADA,1985) FOR MULTIPLE
!		 FAULT PLANE CONFIGURATION.
!     #. STEREOGRAPHIC PROJECTION IS IMPLEMENTED TO CREATE MORE
!		 ACCURATE MAPPING BETWEEN THE EARTH SURFACE AND THE PLANE USED
!		 IN OKADA (1985)
!INPUT:
!OUTPUT:
!NOTES:
!     #. CREATED ON DEC 21 2008 (XIAOMING WANG, GNS)
!     #. UPDATED ON JAN 05 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE FAULT_PARAMS
      TYPE (LAYER) 	:: LO
	  TYPE (LAYER),DIMENSION(NUM_GRID)  :: LA
      TYPE (FAULT), DIMENSION(NUM_FLT)  :: FLT
      REAL TIME, TEMP(LO%NX)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
	  TEMP = 0.0

	  NUM_FAULT = FLT(1)%NUM_FLT

	  DO K = 1, NUM_FAULT
	     IF (FLT(K)%T0.GE.TIME .AND. FLT(K)%T0.LT.TIME+LO%DT) THEN
             WRITE(*,*) "[WARNING] IF RAPTURE TIME != 0 MAY CAUSE DATA INCONSISTENCY BETWEEN HOST AND DEVICE"
             WRITE(*,*) "[WARNING] SHOULDN'T IGNORE THIS MESSAGE IF IT APPEARS IN THE SIMULATION LOOP"
			!CALCULATE SEAFLOOR DEFORMATION VIA FINITE FAULT THEORY
		    IF (FLT(K)%SWITCH.EQ.0 .OR. FLT(K)%SWITCH.EQ.9) THEN
			   IF (LO%INI_SWITCH .EQ. 0) THEN
	              CALL DEFORM_OKADA(LO,FLT(K))
			   ELSEIF (LO%INI_SWITCH .EQ. 9) THEN
				  CALL DEFORM_SMYLIE (LO,FLT(K))
			   ENDIF
			ENDIF
			!OBTAIN SEAFLOOR DEFORMATION FROM A DATA FILE
			IF (FLT(K)%SWITCH .EQ. 1) THEN
			   IF (FLT(K)%FS.EQ.0) CALL READ_COMCOT_DEFORM (LO,FLT(K))
			   IF (FLT(K)%FS.EQ.1) CALL READ_MOST_DEFORM (LO,FLT(K))
			   IF (FLT(K)%FS.EQ.2) CALL READ_XYZ_DEFORM (LO,FLT(K))
			ENDIF

			!WRITE SEAFLOOR DEFORMATION PROFILE INTO FILES
			IF (FLT(K)%FS .NE. 9) CALL WRITE_DEFORM (LO,K)

			!INTERPOLATING DEFORMATION FROM 1ST-LEVEL TO ALL SUB-LEVELS
			CALL ININTERP (LO,LA)

			!UPDATE BATHYMETRY (DUE TO SEAFLOOR DEFORMATION)
			CALL UPDATE_BATH (LO,LA)
			!APPLY SEAFLOOR DEFORMATION ONTO WATER SURFACE (INSTANTLY)
			IF (LO%LAYSWITCH .EQ. 0) THEN
               WRITE(*,*)"[WARNING] WRITING INFORMATION TO Z, SHOULDN'T IGNORE THIS MESSAGE IF IT APPERS IN THE SIMULATION LOOP"
			   LO%Z(:,:,1) = LO%Z(:,:,1) + LO%DEFORM(:,:)
			ENDIF
			DO I = 1,NUM_GRID
			   IF (LA(I)%LAYSWITCH .EQ. 0) THEN
		          LA(I)%Z(:,:,1) = LA(I)%Z(:,:,1) + LA(I)%DEFORM(:,:)

			   ENDIF
			ENDDO

		 ENDIF
	  ENDDO

	  RETURN
	  END

!----------------------------------------------------------------------
      SUBROUTINE DEFORM_OKADA (LO,FLT)
!......................................................................
!DESCRIPTION:
!     #. CALCULATE SEAFLOOR DEFORMATION VIA OKADA'S MODEL (1985);
!     #. STEREOGRAPHIC PROJECTION IS IMPLEMENTED TO CREATE MORE
!		 ACCURATE MAPPING BETWEEN THE EARTH SURFACE AND THE PLANE USED
!		 IN OKADA (1985)
!INPUT:
!	  #. FAULT PARAMETERS;
!OUTPUT:
!	  #. SEAFLOOR DEFORMATION;
!NOTES:
!     #. CREATED ON JUN ?? 2003 (XIAOMING WANG, CORNELL UNIVERSITY)
!     #. UPDATED ON DEC 18, 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON FEB 03 2009 (XIAOMING WANG, GNS)
!		 1. ADD DETECTION ON NAN/INF
!	  #. UPDATED ON FEB 16 2009 (XIAOMING WANG, GNS)
!		 1. ADD AN OPTION TO SELECT THE FOCUS LOCATION
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE FAULT_PARAMS
      TYPE (LAYER) 	:: LO
      TYPE (FAULT)  :: FLT
      REAL L, TEMP(LO%NX)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA OSIXTY/0.016666666667/, BIG/-999./
      DATA RAD/0.01745329252/

	  TEMP = 0.0

      NX = LO%NX
	  NY = LO%NY

!	  WRITE (*,*) FLT%D

      ANG_L = RAD_DEG*FLT%DL
      ANG_R = RAD_DEG*FLT%RD
      ANG_T = RAD_DEG*FLT%TH
      HALFL = 0.5*FLT%L
!.....CALCULATE FOCAL DEPTH USED FOR OKADA'S MODEL
	  HH = FLT%HH+0.5*FLT%W*SIN(ANG_L)
!.....DISPLACEMENT DUE TO DIFFERENT EPICENTER DEFINITION
!	  EPICENTER IS DEFINED AT THE CENTER OF FAULT PLANE
      DEL_X = 0.5*FLT%W*COS(ANG_L)*COS(ANG_T)
	  DEL_Y = 0.5*FLT%W*COS(ANG_L)*SIN(ANG_T)

      H1 = HH/SIN(ANG_L)
      H2 = HH/SIN(ANG_L)+FLT%W

!     EPICENTER IS DEFINED AT THE CENTER OF TOP EDGE OF FAULT PLANE
	  IF (FLT%SWITCH .EQ. 9) THEN
         HH = FLT%HH+FLT%W*SIN(ANG_L)
	     DEL_X = FLT%W*COS(ANG_L)*COS(ANG_T)
	     DEL_Y = FLT%W*COS(ANG_L)*SIN(ANG_T)
         H1 = HH/SIN(ANG_L)
         H2 = HH/SIN(ANG_L)+FLT%W
	  ENDIF

      DS = FLT%D*COS(ANG_R)
      DD = FLT%D*SIN(ANG_R)

      SN = SIN(ANG_L)
      CS = COS(ANG_L)

	  LO%DEFORM = 0.0
	  X_SHIFT = 0.0
	  Y_SHIFT = 0.0

      DO I=1,NX
         DO J=1,NY
		    IF (LO%LAYCORD .EQ. 0) THEN
			   CALL STEREO_PROJECTION (X_SHIFT,Y_SHIFT,LO%X(I),	&
										  LO%Y(J),FLT%X0,FLT%Y0)
!*			   CALL SPH_TO_UTM (X_SHIFT,Y_SHIFT,LO%X(I),LO%Y(J),&
!*								FLT%X0,FLT%Y0)
			   X_SHIFT = X_SHIFT - DEL_X
			   Y_SHIFT = Y_SHIFT + DEL_Y
		    ELSE
               YY = LO%DY*(FLOAT(J)-1.0)
               XX = LO%DX*(FLOAT(I)-1.0)
			   CALL STEREO_PROJECTION (XL,YL,FLT%XO,FLT%YO,		&
										  FLT%X0,FLT%Y0)
!*			   CALL SPH_TO_UTM (XL,YL,FLT%XO,FLT%YO,FLT%X0,FLT%Y0)
			   X_SHIFT = XX - XL - DEL_X
			   Y_SHIFT = YY - YL + DEL_Y
		    ENDIF
            X1 = X_SHIFT*SIN(ANG_T)+Y_SHIFT*COS(ANG_T)
            X2 = X_SHIFT*COS(ANG_T)-Y_SHIFT*SIN(ANG_T)
		    X2 = -X2
            X3 = ZERO
            P = X2*CS+HH*SN
            CALL STRIKE_SLIP (X1,X2,X3,X1+HALFL,P,ANG_L,HH,F1)
            CALL STRIKE_SLIP (X1,X2,X3,X1+HALFL,P-FLT%W,ANG_L,HH,F2)
            CALL STRIKE_SLIP (X1,X2,X3,X1-HALFL,P,ANG_L,HH,F3)
            CALL STRIKE_SLIP (X1,X2,X3,X1-HALFL,P-FLT%W,ANG_L,HH,F4)
            CALL DIP_SLIP (X1,X2,X3,X1+HALFL,P,ANG_L,HH,G1)
            CALL DIP_SLIP (X1,X2,X3,X1+HALFL,P-FLT%W,ANG_L,HH,G2)
            CALL DIP_SLIP (X1,X2,X3,X1-HALFL,P,ANG_L,HH,G3)
            CALL DIP_SLIP (X1,X2,X3,X1-HALFL,P-FLT%W,ANG_L,HH,G4)
            US = (F1-F2-F3+F4)*DS
            UD = (G1-G2-G3+G4)*DD

			IF (ABS(US) .GE. HUGE(US)) US = 0.0
			IF (ABS(UD) .GE. HUGE(UD)) UD = 0.0
			IF (ISNAN(US)) US = 0.0
			IF (ISNAN(UD)) UD = 0.0
			IF (ABS(US) .LE. GX) US = 0.0
			IF (ABS(UD) .LE. GX) UD = 0.0

            LO%DEFORM(I,J) = US + UD
         END DO
      END DO

!	  WRITE (*,*) 'SUBROUTINE OKADA HAS BEEN CALLED'

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE STRIKE_SLIP (X1,X2,X3,Y1,Y2,DP,DD,F)
!.....USED FOR OKADA'S MODEL (CREATED BY XIAOMING WANG IN JUN 2003)
!NOTE:
!	 #. UPDATED ON FEB04 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      SN = SIN(DP)
      CS = COS(DP)
      P = X2*CS + DD*SN
      Q = X2*SN - DD*CS
      D_BAR = Y2*SN - Q*CS
      R = SQRT(Y1**2 + Y2**2 + Q**2)
      XX = SQRT(Y1**2 + Q**2)
!*      A4 = 0.5*1/CS*(LOG(R+D_BAR) - SN*LOG(R+Y2))
	  TMP1 = R+D_BAR
	  TMP2 = R+Y2
	  IF (TMP1 .LE. EPS) TMP1 = EPS
	  IF (TMP2 .LE. EPS) TMP2 = EPS
      A4 = 0.5*1/CS*(LOG(TMP1) - SN*LOG(TMP2))
      F = -(D_BAR*Q/R/(R+Y2) + Q*SN/(R+Y2) + A4*SN)/2/3.14159

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE DIP_SLIP (X1,X2,X3,Y1,Y2,DP,DD,F)
!.....BASED ON OKADA'S PAPER (1985)
!.....CREATED BY XIAOMING WANG (JUN 2003)
!----------------------------------------------------------------------
      SN = SIN(DP)
      CS = COS(DP)

      P = X2*CS + DD*SN
      Q = X2*SN - DD*CS
      D_BAR = Y2*SN - Q*CS;
      R = SQRT(Y1**2 + Y2**2 + Q**2)
      XX = SQRT(Y1**2 + Q**2)
      A5 = 0.5*2/CS*ATAN((Y2*(XX+Q*CS)+XX*(R+XX)*SN)/Y1/(R+XX)/CS)
      F = -(D_BAR*Q/R/(R+Y1) + SN*ATAN(Y1*Y2/Q/R)					&
							 - A5*SN*CS)/2/3.14159

      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE DEFORM_SMYLIE (LAY,FLT)
!......................................................................
!DESCRIPTION:
!     #. CALCULATE SEAFLOOR DEFORMATION VIA THE ELASTIC FAULT PLANE
!		 MODEL BY MANSINHA AND SMYLIE (1971);
!     #. STEREOGRAPHIC PROJECTION IS IMPLEMENTED TO CREATE MORE
!		 ACCURATE MAPPING BETWEEN THE EARTH SURFACE AND THE PLANE USED
!		 BY MANSINHA AND SMYLIE (1971);
!INPUT:
!	  #. FAULT PARAMETERS;
!OUTPUT:
!	  #. SEAFLOOR DEFORMATION;
!NOTES:
!     #. UPDATED ON AUG 13 2003 (XIAOMING WANG, CORNELL UNIVERSITY)
!        1. BOTH CARTESIAN AND SPHERICAL COORD. ARE SUPPORTED
!        2. THE CALCULATION OF EPICENTER IS MODIFIED ;
!     #. UPDATED ON DEC 18, 2008 (XIAOMING WANG, GNS)
!	     1. OBLIQUE STEREOGRAPHIC PROJECTION IS IMPLEMENTED TO CREATE
!		    MORE ACCURATE MAPPING BETWEEN THE EARTH SURFACE AND THE
!		    PLANE USED IN MANSINHA AND SMYLIE (1971);
!		 2. EPICENTER IS CHOSEN AS THE TANGENTIAL POINT OF PROJECTION;
!	  #. UPDATED ON FEB 03 2009 (XIAOMING WANG, GNS)
!		 1. ADD DETECTION ON NAN OR INF
!	  #. UPDATED ON FEB 16 2009 (XIAOMING WANG, GNS)
!		 1. ADD AN OPTION TO SELECT THE FOCUS LOCATION
!----------------------------------------------------------------------
      USE LAYER_PARAMS
	  USE FAULT_PARAMS
      TYPE (LAYER) 	:: LAY
      TYPE (FAULT)  :: FLT
      REAL L, TEMP(LAY%NX)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA OSIXTY/0.016666666667/, BIG/-999./
      DATA RAD/0.01745329252/

	  TEMP = 0.0
!*	  ! CONVERT DEGREE TO RADIAN ! ON THE EAST SIDE BY TSO-REN
!*      !ANG_L = RAD*(180-FLT%DL)
!.....CONVERT DEGREE TO RADIAN
      ANG_L = RAD_DEG*FLT%DL
      ANG_R = RAD_DEG*FLT%RD
	  ANG_T = RAD_DEG*FLT%TH
      HALFL = 0.5*FLT%L

!.. . DEL_X,DEL_Y, SHIFTS DUE TO DIFFERENT EPICENTER DEFINITIONS
!     THE X, Y DISTANCE BETWEEN THE EPICENTER AND THE BREAK POINT

!     EPICENTER IS DEFINED AT THE CENTER OF FAULT PLANE
	  HH = FLT%HH - 0.5*FLT%W*SIN(ANG_L)
	  HTMP = FLT%HH
      DEL_X = HTMP*COS(ANG_L)/SIN(ANG_L)*COS(ANG_T)
      DEL_Y = HTMP*COS(ANG_L)/SIN(ANG_L)*SIN(ANG_T)
      H1 = HH/SIN(ANG_L)
      H2 = HH/SIN(ANG_L)+FLT%W

!	  EPICENTER IS DEFINED AT THE CENTER OF TOP EDGE OF FAULT PLANE
	  IF (FLT%SWITCH .EQ. 9) THEN
	     HTMP = FLT%HH
         DEL_X = HTMP*COS(ANG_L)/SIN(ANG_L)*COS(ANG_T)
         DEL_Y = HTMP*COS(ANG_L)/SIN(ANG_L)*SIN(ANG_T)
         H1 = FLT%HH/SIN(ANG_L)
         H2 = FLT%HH/SIN(ANG_L)+FLT%W
	  ENDIF

      DS = FLT%D*COS(ANG_R)
      DD = FLT%D*SIN(ANG_R)

      O12PI = 1./(12.*PI)

      NX = LAY%NX
      NY = LAY%NY
	  ! -----------------------------------------------
	  LAY%DEFORM = 0.0
      DO I=1,NX
!.......CALCULATE THE DISTANCE BETWEEN THE ORIGIN OF COMPUTATIONAL
!		DOMAIN AND THE ORIGIN OF FOCAL COORDINATE SYSTEM
         DO J=1,NY
		    IF (LAY%LAYCORD .EQ. 0) THEN
			   CALL STEREO_PROJECTION (X_SHIFT,Y_SHIFT,LAY%X(I),	&
										LAY%Y(J),FLT%X0,FLT%Y0)
			   X_SHIFT = X_SHIFT + DEL_X
			   Y_SHIFT = Y_SHIFT - DEL_Y
		    ELSE
               YY = LAY%DY*(FLOAT(J)-1.0)
               XX = LAY%DX*(FLOAT(I)-1.0)
			   CALL STEREO_PROJECTION (XL,YL,FLT%XO,FLT%YO,		&
										FLT%X0,FLT%Y0)
			   X_SHIFT = XX - XL + DEL_X
			   Y_SHIFT = YY - YL - DEL_Y
		    ENDIF
!...........CONVERTED TO FOCAL COORDINATES
            X1 = -(X_SHIFT*SIN(ANG_T)+Y_SHIFT*COS(ANG_T))
            X2 = X_SHIFT*COS(ANG_T)-Y_SHIFT*SIN(ANG_T)
            X3 = ZERO
            CALL USCAL (X1,X2,X3, HALFL,H2,ANG_L,F1)
            CALL USCAL (X1,X2,X3, HALFL,H1,ANG_L,F2)
            CALL USCAL (X1,X2,X3,-HALFL,H2,ANG_L,F3)
            CALL USCAL (X1,X2,X3,-HALFL,H1,ANG_L,F4)
            CALL UDCAL (X1,X2,X3, HALFL,H2,ANG_L,G1)
            CALL UDCAL (X1,X2,X3, HALFL,H1,ANG_L,G2)
            CALL UDCAL (X1,X2,X3,-HALFL,H2,ANG_L,G3)
            CALL UDCAL (X1,X2,X3,-HALFL,H1,ANG_L,G4)
            US = (F1-F2-F3+F4)*DS*O12PI
            UD = (G1-G2-G3+G4)*DD*O12PI

			IF (ABS(US).GE.HUGE(US)) US = 0.0
			IF (ABS(UD).GE.HUGE(UD)) UD = 0.0
			IF (ISNAN(US)) US = 0.0
			IF (ISNAN(UD)) UD = 0.0
			IF (ABS(US) .LE. GX) US = 0.0
			IF (ABS(UD) .LE. GX) UD = 0.0

            LAY%DEFORM(I,J) = US + UD
         END DO
      END DO
!
      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE USCAL (X1,X2,X3,C,CC,DP,F)
!.....CALCULATE STRIKE SLIP, USED FOR MANSINHA AND SMYLIE'S MODEL
!NOTE:
!	 #. UPDATED ON FEB04 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      SN = SIN(DP)
      CS = COS(DP)
      C1 = C
      C2 = CC*CS
      C3 = CC*SN
      R = SQRT((X1-C1)**2+(X2-C2)**2+(X3-C3)**2)
      Q = SQRT((X1-C1)**2+(X2-C2)**2+(X3+C3)**2)
      R2 = X2*SN-X3*CS
      R3 = X2*CS+X3*SN
      Q2 = X2*SN+X3*CS
      Q3 = -X2*CS+X3*SN
      H = SQRT(Q2**2+(Q3+CC)**2)
      K = SQRT((X1-C1)**2+Q2**2)
!*      A1 = LOG(R+R3-CC)
!*      A2 = LOG(Q+Q3+CC)
!*      A3 = LOG(Q+X3+C3)
	  TMP1 = R+R3-CC
	  TMP2 = Q+Q3+CC
	  TMP3 = Q+X3+C3
	  IF (TMP1 .LE. EPS) TMP1 = EPS
	  IF (TMP2 .LE. EPS) TMP2 = EPS
	  IF (TMP3 .LE. EPS) TMP3 = EPS
	  A1 = LOG(TMP1)
	  A2 = LOG(TMP2)
	  A3 = LOG(TMP3)
      B1 = 1+3.0*(TAN(DP))**2
      B2 = 3.0*TAN(DP)/CS
      B3 = 2.0*R2*SN
      B4 = Q2+X2*SN
      B5 = 2.0*R2**2*CS
      B6 = R*(R+R3-CC)
      B7 = 4.0*Q2*X3*SN**2
      B8 = 2.0*(Q2+X2*SN)*(X3+Q3*SN)
      B9 = Q*(Q+Q3+CC)
      B10 = 4.0*Q2*X3*SN
      B11 = (X3+C3)-Q3*SN
      B12 = 4.0*Q2**2*Q3*X3*CS*SN
      B13 = 2.0*Q+Q3+CC
      B14 = Q**3*(Q+Q3+CC)**2
      F = CS*(A1+B1*A2-B2*A3)+B3/R+2*SN*B4/Q-B5/B6+(B7-B8)/B9		&
			+ B10*B11/Q**3-B12*B13/B14
	  IF (ABS(R*Q*B9*B14) .LT. 1.0E-10) THEN
	     WRITE (*,*) 'DIVIDED BY ZERO'
	  END IF
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE UDCAL (X1,X2,X3,C,CC,DP,F)
!.....CALCULATE DIP-SLIP, USED FOR MANSINHA AND SMYLIE'S MODEL
!NOTE:
!	 #. UPDATED ON FEB04 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      SN = SIN(DP)
      CS = COS(DP)
      C1 = C
      C2 = CC*CS
      C3 = CC*SN
      R = SQRT((X1-C1)**2+(X2-C2)**2+(X3-C3)**2)
      Q = SQRT((X1-C1)**2+(X2-C2)**2+(X3+C3)**2)
      R2 = X2*SN-X3*CS
      R3 = X2*CS+X3*SN
      Q2 = X2*SN+X3*CS
      Q3 = -X2*CS+X3*SN
      H = SQRT(Q2**2+(Q3+CC)**2)
      K = SQRT((X1-C1)**2+Q2**2)
!*      A1 = LOG(R+X1-C1)
!*      A2 = LOG(Q+X1-C1)
	  TMP1 = R+X1-C1
	  TMP2 = Q+X1-C1
	  IF (TMP1 .LE. EPS) TMP1 = EPS
	  IF (TMP2 .LE. EPS) TMP2 = EPS
	  A1 = LOG(TMP1)
	  A2 = LOG(TMP2)
      B1 = Q*(Q+X1-C1)
      B2 = R*(R+X1-C1)
      B3 = Q*(Q+Q3+CC)
      D1 = X1-C1
      D2 = X2-C2
      D3 = X3-C3
      D4 = X3+C3
      D5 = R3-CC
      D6 = Q3+CC
      T1 = ATN(D1*D2,(H+D4)*(Q+H))
      T2 = ATN(D1*D5,R2*R)
      T3 = ATN(D1*D6,Q2*Q)
      F = SN*(D2*(2*D3/B2+4*D3/B1-4*C3*X3*D4*(2*Q+D1)/(B1**2*Q))	  &
			- 6*T1+3*T2-6*T3)+CS*(A1-A2-2*(D3**2)/B2				  &
			- 4*(D4**2-C3*X3)/B1-4*C3*X3*D4**2*(2*Q+X1-C1)/(B1**2*Q)) &
			+ 6*X3*(CS*SN*(2*D6/B1+D1/B3)-Q2*(SN**2-CS**2)/B1)

	  IF (ABS(B1**2*Q*B2*B3) .LT. 1.0E-10) THEN
	      WRITE(*,*) 'DIVIDED BY ZERO IN UDCAL'
	  ENDIF

      RETURN
      END

!----------------------------------------------------------------------
      REAL FUNCTION ATN (AX,AY)
!----------------------------------------------------------------------
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      AAX = ABS(AX)
      AAY = ABS(AY)
      P = AX*AY
      IF ((AAX.LE.EPS).AND.(AAY.LE.EPS)) THEN
        ATN = 0.2
        WRITE(6,1) AX,AY
 1      FORMAT('ATAN FIX --  AX =',E15.7,2X,'AY =',E15.7)
      ELSE
        SR = ATAN2(AAX,AAY)
        ATN = SIGN(SR,P)
      ENDIF
!
      RETURN
      END

!----------------------------------------------------------------------
      SUBROUTINE STEREO_PROJECTION (X,Y,LONIN,LATIN,LON0,LAT0)
!......................................................................
!DESCRIPTION:
!     # MAPPING A POINT ON THE ELLIPSOID SURFACE ONTO A PLANE;
!     # OBLIQUE STEREOGRAPHIC PROJECTION IS ADOPTED
!INPUT:
!     LATIN: LATITUDE IN DEGREES
!     LONIN: LONGITUDE IN DEGREES
!     LAT0: LATITUDE OF TANGENTIAL POINT IN DEGREES (E.G., EPICENTER)
!     LON0: LONGITUDE OF TANGENTIAL POINT IN DEGREES (E.G., EPICENTER)
!OUTPUT:
!     X: X COORDINATE/EASTING IN METERS RELATIVE TO ORIGIN (I.E., LON0)
!     Y: Y COORDINATE/NORTHING IN METERS RELATIVE TO ORIGIN (I.E., LAT0)
!REFERENCES:
!	  #. SNYDER, J.P. (1987). MAP PROJECTIONS - A WORKING MANUAL.
!                          USGS PROFESSIONAL PAPER 1395
!     #. ELLIPSOIDAL DATUM: WGS84
!WORKING NOTES:
!     CREATED ON DEC18 2008 (XIAOMING WANG, GNS)
!     UPDATED ON JAN02 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      REAL XOUT,YOUT,XS,YS,LAT0,LON0,LATIN,LONIN
	  REAL COS_X,COS_Y,SIN_X,SIN_Y
	  REAL LAT,LON,LT0,LN0,CS,SN,CS0,SN0,TMP,TMP0
	  REAL A,B,K0,E,ES,N,C,R,S1,S2,W1,W2,SA,SB,BETA,XI,LM,XI0,LM0
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  POLE = PI/2.0 - EPS	  !AVOID SINGULARITY AT POLES

	  ! CONVERT DEGREE TO RADIAN
	  LAT = LATIN*RAD_DEG
	  LON = LONIN*RAD_DEG
	  LT0 = LAT0*RAD_DEG
	  LN0 = LON0*RAD_DEG
	  IF (LAT .GT. POLE) LAT = POLE
	  IF (LAT .LT. -POLE) LAT = -POLE
	  IF (LT0 .GT. POLE) LT0 = POLE
	  IF (LT0 .LT. -POLE) LT0 = -POLE

	  CS = COS(LAT)
	  SN = SIN(LAT)
	  CS0 = COS(LT0)
	  SN0 = SIN(LT0)

	  ! PARAMETERS
	  XF = 0.0					! FALSE EASTING
	  YF = 0.0					! FALSE NORTHING

	  A = 6378137.0000          ! ELLIPSOIDAL SEMI-MAJOR AXIS
	  B = 6356752.3142			! ELLIPSOIDAL SEMI-MINOR AXIS
	  F = 0.00335281067183      ! FLATTENING, F = (A-B)/A
	  E = 0.08181919092891		! ECCENTRICITY, E = SQRT(2.0*F-F**2)
	  F2 = 0.00669438000426		! F2 = E**2
	  ES = 0.00673949675659		! 2ND ECCENTRICITY, ES = E**2/(1-E**2)

	  K0 = 0.9996				! SCALE FACTOR

	  TMP = SQRT(1.0-F2*SN**2)
	  TMP0 = SQRT(1.0-F2*SN0**2)
	  RHO0 = A*(1.0-F2)/TMP0**3
	  NU0 = A/TMP0
	  R = SQRT(RHO0*NU0)
	  N = SQRT(1.0+F2*CS0**4/(1.0-F2))

	  S1 = (1.0+SN0)/(1.0-SN0)
	  S2 = (1.0-E*SN0)/(1.0+E*SN0)
	  W1 = (S1*S2**E)**N
	  SN_XI0 = (W1-1.0)/(W1+1.0)
	  C = (N+SN0)*(1.0-SN_XI0)/(N-SN0)/(1.0+SN_XI0)

	  W2 = C*W1
	  SA = (1.0+SN)/(1.0-SN)
	  SB = (1.0-E*SN)/(1.0+E*SN)
	  W = C*(SA*SB**E)**N

	  XI0 = ASIN((W2-1.0)/(W2+1.0))
	  LM0 = LN0

	  LM = N*(LON-LM0)+LM0
	  XI = ASIN((W-1.0)/(W+1.0))

	  BETA = 1.0 + SIN(XI)*SIN(XI0) + COS(XI)*COS(XI0)*COS(LM-LM0)

	  Y = YF + 2.0*R*K0*(SIN(XI)*COS(XI0) &
						 - COS(XI)*SIN(XI0)*COS(LM-LM0))/BETA
	  X = XF + 2.0*R*K0*COS(XI)*SIN(LM-LM0)/BETA

!	  WRITE(*,*) LM,XI,X,Y

      RETURN
      END


!//////////////////////////////////////////////////////////////////////
! # SEAFLOOR DEFORMATION DATA, BATHYMETRIC AND TOPOGRAPHICAL DATA ARE
!   LOADED INTO COMCOT IN THE FOLLOWING SUBROUTINES
!//////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------------
      SUBROUTINE READ_COMCOT_DEFORM (LO,FAULT_INFO)
!.....READ DEFORMATION DATA FORMATTED FOR MOST MODEL
!     LAST REVISE: NOV 21 2008 (XIAOMING WANG)
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      USE FAULT_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (FAULT)	:: FAULT_INFO
	  REAL Z(LO%NX,LO%NY)
      INTEGER      :: ISTAT, IS, JS, IE, JE

	  IS = 1
	  JS = 1
	  IE = LO%NX
	  JE = LO%NY
	  Z = 0.0

      WRITE (*,*) '    READING COMCOT-FORMATTED DEFORMATION FROM FILE'
!	  WRITE (*,*) '       ',FAULT_INFO%DEFORM_NAME
      OPEN(99,FILE=FAULT_INFO%DEFORM_NAME,STATUS='OLD',IOSTAT=ISTAT)
      IF (ISTAT /=0) THEN
         PRINT *,"ERROR:: CAN'T OPEN DEFORMATION DATA FILE; EXITING."
         STOP
      END IF
      DO J=JS,JE
         READ (99,'(15F8.3)') (Z(I,J),I=IS,IE)
      END DO
      CLOSE (99)
	  LO%DEFORM(:,:) = Z(:,:)

!*	  DO I = IS,IE
!*	     DO J = JS,JE
!*		    IF ( ABS(LO%DEFORM(I,J)).LT.0.001 ) LO%DEFORM(I,J) = 0.0
!*		 END DO
!*	  END DO

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE READ_MOST_DEFORM (LO,FAULT_INFO)
!......................................................................
!DESCRIPTION:
!	  #. READ MOST-FORMATTED DEFORMATION DATA;
!     #. FIRST ROW OF DATA FILE DESCRIBE THE DIMENSION OF THE
!	     DEFORMATION GRID: NX, NY; LONGITUDE AND LATTIUDE OF THE CENTER
!		 OF DEFORMATION REGION AND ITS INDICES IN NUMERICAL DOMAIN;
!     #. GRID DATA IS WRITTEN ROW BY ROW (X FIRST)
!     #. NODATA TYPE, NAN, IS NOT ALLOWED
!     #. DEFORMATION DATA SHOULD HAVE THE SAME RESOLUTION AS THAT OF
!		 1ST-LEVEL GRIDS;
!NOTES:
!	  #. CREATED ON OCT 25 2008 (XIAOMING WANG, GNS)
!     #. UPDATED ON NOV 21 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      USE FAULT_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (FAULT)	:: FAULT_INFO
      INTEGER STAT, IS, JS, I, J, i0, j0, NX, NY
      INTEGER LENGTH, RC !, FLAG
	  REAL Z(LO%NX,LO%NY), X0, Y0
	  REAL,ALLOCATABLE :: TMP(:,:),TMP1(:,:)
	  REAL,ALLOCATABLE :: X(:),Y(:)
      CHARACTER(LEN=20) FNAME
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

      WRITE (*,*) '    READING MOST-FORMATTED DEFORMATION DATA FROM FILE'
!	  WRITE (*,*) '       ',FAULT_INFO%DEFORM_NAME

      OPEN (UNIT=23,FILE=FAULT_INFO%DEFORM_NAME,STATUS='OLD',		&
									IOSTAT=ISTAT,FORM='FORMATTED')
      IF (ISTAT /=0) THEN
         PRINT *,"ERROR:: CAN'T OPEN ",FAULT_INFO%DEFORM_NAME,"; EXITING."
         STOP
      END IF

      READ (23,*) NX,NY,X0,Y0,I0,J0
!	  WRITE(*,*) NX,NY,X0,Y0,I0,J0

!     (X0,Y0) DEFINES THE CENTER OF THE DEFORMATION AREA
!     (NX,NY) IS THE DIMENSION OF THE DEFORMATION GRIDS
      ALLOCATE(TMP(NX,NY))
	  ALLOCATE(TMP1(NX,NY))
	  ALLOCATE(X(NX))
	  ALLOCATE(Y(NY))
	  TMP = 0.0
	  TMP1 = 0.0
	  X = 0.0
	  Y = 0.0
	  Z = 0.0

      IF(X0.LT.LO%X(1) .OR. X0.GT.LO%X(LO%NX) .OR. 					&
	  		Y0.LT.LO%Y(1) .OR. Y0.GT.LO%Y(LO%NY)) THEN
         WRITE(*,*) 'PROBLEM: DEFORMATION AREA BEYOND THE DOMAIN !!!!'
      END IF

      DO J=1,NY
        READ (23,*) (TMP(I,J), I=1,NX)
		DO I = 1,NX
		   IF (ISNAN(TMP(I,J))) TMP(I,J) = 0.0
		   IF (ABS(TMP(I,J)).GE.HUGE(TMP(I,J))) TMP(I,J) = 0.0
		ENDDO
      END DO
      CLOSE (23)
	 !CREATE COORDINATES FOR THE DEFORMATION AREA
      DO I = 1,NX
		 X(I) = X0 - NX/2.0*LO%DX/60.0 + (I-1)*LO%DX/60.0
      ENDDO
	  DO J = 1,NY
	     Y(J) = Y0 - NY/2.0*LO%DY/60.0 + (J-1)*LO%DY/60.0
      ENDDO
!	  WRITE(*,*) X(1),X(NX)
!	  WRITE(*,*) Y(1),Y(NY)
!.....CONVERT THE FORMAT FROM MOST COORDINATES INTO COMCOT COORDINATES
!     NOTE: IN MOST DATA, Y POINTING TO THE SOUTH
!           IN COMCOT, Y POINTING TO THE NORTH
!!....DATA NEED TO FLIP
      ! FLIP DEFORMATION MATRIX
      DO I = 1,NX
	     DO J = 1,NY
		    K = NY - J + 1
			TMP1(I,K) = TMP(I,J)
		 END DO
	  END DO

!.....MAPPING THE DEFORM DATA ONTO THE NUMERICAL GRIDS
!*      CALL DEFORM_MAPPING (LO,TMP1,X,Y,NX,NY)
      CALL GRID_INTERP (Z,LO%X,LO%Y,LO%NX,LO%NY,TMP1,X,Y,NX,NY)

	  CALL SMOOTH_BC (Z,LO%X,LO%Y,LO%NX,LO%NY)

	  LO%DEFORM(:,:) = Z(:,:)

      DEALLOCATE(TMP,TMP1,X,Y,STAT = ISTAT)

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE READ_XYZ_DEFORM (LO,FAULT_INFO)
!......................................................................
!DESCRIPTION:
!	  #. READ XYZ-FORMAT DEFORMATION DATA;
!     #. GRIDDED DEFORMATION DATA CONTAINS 3 COLUMNS: X COORDINATES,
!		 Y COORDINATES, DEFORMATION (Z);
!	  #. COORDINATE SYSTEM IS DEFINED SO THAT X POINTING TO THE EAST
!		 (LONGITUDE) AND Y AXIS POINTING TO THE NORTH (LATITUDE);
!     #. GRID DATA IS WRITTEN ROW BY ROW (X FIRST) FROM WEST TO EAST,
!		 FROM SOUTH TO NORTH (OR FOR NORTH TO SOUTH);
!     #. NODATA TYPE, NAN, IS NOT ALLOWED (FOR COMPAQ COMPILER)
!NOTES:
!	  #. CREATED ON NOV 05 2008 (XIAOMING WANG, GNS)
!     #. UPDATED ON NOV 21 2008 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      USE FAULT_PARAMS
      TYPE (LAYER) :: LO
	  TYPE (FAULT)	:: FAULT_INFO
      INTEGER      :: ISTAT, IS, JS, IE, JE
	  REAL Z(LO%NX,LO%NY)
	  REAL,ALLOCATABLE :: TMP(:,:),DEFORM(:,:)
	  REAL,ALLOCATABLE :: XCOL(:),YCOL(:),ZCOL(:)
	  REAL,ALLOCATABLE :: X(:),Y(:),XTMP(:),YTMP(:)
      INTEGER	   LENGTH, RC, POS !, FLAG
	  INTEGER      COUNT
	  REAL         TEMP,TEMP1,TEMP2,TEMP3
      INTEGER :: RSTAT = 0
!      CHARACTER(LEN=20) FNAME
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  Z = 0.0
	  COUNT = -1

      WRITE (*,*) '    READING XYZ-FORMAT DEFORMATION DATA FILE...'
!	  WRITE (*,*) '       ',FAULT_INFO%DEFORM_NAME
      OPEN (UNIT=23,FILE=FAULT_INFO%DEFORM_NAME,STATUS='OLD',		&
									IOSTAT=ISTAT,FORM='FORMATTED')
      IF (ISTAT /=0) THEN
         PRINT *,"ERROR:: CAN'T OPEN DEFORMATION DATA FILE; EXITING."
         STOP
      END IF

!.....DETERMINE THE LENGTH OF BATHYMETRY DATA FILE
	  DO WHILE (RSTAT == 0)
		 COUNT = COUNT + 1
		 READ (23,*,IOSTAT=RSTAT) TEMP1,TEMP2,TEMP3
	  ENDDO
!*      CLOSE(23)
!*!	  WRITE (*,*) '       TOTAL NUMBER OF BATHYMETRY POINTS: ', COUNT
      NXY = COUNT
	  ALLOCATE(XCOL(NXY))
	  ALLOCATE(YCOL(NXY))
	  ALLOCATE(ZCOL(NXY))
	  XCOL = 0.0
	  YCOL = 0.0
	  ZCOL = 0.0

!*!.....READING DEFORMATION DATA
      REWIND(23)
	  DO I = 1,COUNT
         READ (23,*) XCOL(I), YCOL(I), ZCOL(I)
		 IF (ISNAN(ZCOL(I))) ZCOL(I) = 0.0
		 IF (ABS(ZCOL(I)).GE.HUGE(ZCOL(I))) ZCOL(I) = 0.0
	  END DO
	  CLOSE (23)

!.....DETERMINE GRID DIMENSION: NX,NY
      TMPX = XCOL(1)
	  TMPX1 = XCOL(2)
	  TMPY = YCOL(1)
	  TMPY1 = YCOL(2)
!	  WRITE(*,*) TMPX1,TMPX

!>>>
!.....CHECK IF THE DATA IS WRITTEN ROW BY ROW
	  IF (TMPX1.GT.TMPX .AND. ABS(TMPY1-TMPY).LT.EPS) THEN
!*	  IF (TMPX1.GT.TMPX .AND. TMPY1.EQ.TMPY) THEN
	     K = 1
	     DO WHILE (TMPX1.GT.TMPX)
	        K=K+1
	        TMPX1 = XCOL(K)
	     ENDDO
	     NX = K-1
!	     WRITE(*,*) NX
	     NY = NINT(DBLE(NXY/NX))

!	     WRITE (*,*) '       GRID DIMENSION OF DEFORMATION: ', NX,NY
         ALLOCATE(X(NX))
	     ALLOCATE(Y(NY))
	     ALLOCATE(YTMP(NY))
	     ALLOCATE(TMP(NX,NY))
	     ALLOCATE(DEFORM(NX,NY))
	     TMP = 0.0
	     X = 0.0
	     Y = 0.0
	     YTMP = 0.0
	     DEFORM = 0.0
!........OBTAINED X,Y COORDINATES
	     X(1:NX) = XCOL(1:NX)
	     DO J = 1,NY
	        K = (J-1)*NX + 1
	        YTMP(J) = YCOL(K)
	     END DO
	     !GENERATE GRID DATA
         DO J=1,NY
	        KS = (J-1)*NX + 1
		    KE = (J-1)*NX + NX
            TMP(1:NX,J) = ZCOL(KS:KE)
         END DO
	  ENDIF
!<<<

!>>>
!.....CHECK IF THE DATA IS WRITTEN COLUMN BY COLUMN
	  TMPX = XCOL(1)
	  TMPX1 = XCOL(2)
	  TMPY = YCOL(1)
	  TMPY1 = YCOL(2)
!	  write (*,*) TMPX,TMPX1,TMPY,TMPY1,NXY
	  IF (ABS(TMPX1-TMPX).LT.EPS .AND. ABS(TMPY1-TMPY).GT.EPS) THEN
!*	  IF (TMPX1.EQ.TMPX .AND. TMPY1.NE.TMPY) THEN
	     K = 1
	     DO WHILE (TMPX1.LE.TMPX)
	        K=K+1
	        TMPX1 = XCOL(K)
	     ENDDO
	     NY = K-1
!	     WRITE(*,*) NX
	     NX = NINT(DBLE(NXY/NY))

!*	     WRITE (*,*) '       GRID DIMENSION OF DEFORMATION: ', NX,NY
         ALLOCATE(X(NX))
	     ALLOCATE(Y(NY))
		 ALLOCATE(XTMP(NX))
	     ALLOCATE(YTMP(NY))
	     ALLOCATE(TMP(NX,NY))
	     ALLOCATE(DEFORM(NX,NY))
	     TMP = 0.0
	     X = 0.0
	     Y = 0.0
	     YTMP = 0.0
	     DEFORM = 0.0
!........OBTAINED X,Y COORDINATES
	     YTMP(1:NY) = YCOL(1:NY)
	     DO I = 1,NX
	        K = (I-1)*NY + 1
	        X(I) = XCOL(K)
	     END DO
	     !GENERATE GRID DATA
         DO I=1,NX
	        KS = (I-1)*NY + 1
		    KE = (I-1)*NY + NY
            TMP(I,1:NY) = ZCOL(KS:KE)
         END DO
	  ENDIF
!<<<

!........
!!....DETERMINE IF THE DATA NEED FLIP
!     E.G., Y COORDINATE IS FROM NORTH TO SOUTH OR FROM SOUTH TO NORTH
!     IFLIP = 0: FLIP; 1: NO FLIP OPERATION
      IFLIP = 0
	  IF (YTMP(NY).LT.YTMP(NY-1)) IFLIP = 1
!	  WRITE (*,*) IFLIP

	  IF (IFLIP .EQ. 1) THEN
         ! FLIP Y COORDINATES
         DO J = 1,NY
	        K = NY-J+1
		    Y(K) = YTMP(J)
	     END DO
         ! FLIP BATHYMETRY MATRIX
         DO I = 1,NX
	        DO J = 1,NY
		       K = NY - J + 1
			   DEFORM(I,K) = TMP(I,J)
		    END DO
	     END DO
	  ELSE
	     Y = YTMP
		 DEFORM = TMP
	  END IF
!      WRITE (*,*) NX,NY,X(1),X(NX),Y(1),Y(NY)
!      WRITE (*,*) LO%NX,LO%NY

!.....MAPPING THE DEFORM DATA ONTO THE NUMERICAL GRIDS
!*      CALL DEFORM_MAPPING (LO,DEFORM,X,Y,NX,NY)
      CALL GRID_INTERP (Z,LO%X,LO%Y,LO%NX,LO%NY,DEFORM,X,Y,NX,NY)

	  CALL SMOOTH_BC (Z,LO%X,LO%Y,LO%NX,LO%NY)

	  LO%DEFORM(:,:) = Z(:,:)

      DEALLOCATE(XCOL,YCOL,ZCOL,STAT = ISTAT)
      DEALLOCATE(TMP,XTMP,YTMP,STAT = ISTAT)
      DEALLOCATE(X,Y,DEFORM,STAT = ISTAT)


      RETURN
      END



!----------------------------------------------------------------------
      SUBROUTINE SMOOTH_BC (Z,X,Y,NX,NY)
!......................................................................
!DESCRIPTION:
!	  #. SMOOTH OUT DEFORMATION NEAR BOUNDARIES TO AVOID REFLECTION
!		 PROBLEMS FROM OPEN BOUNDARY;
!	  #.
!NOTES:
!     #. CREATED ON DEC 22 2008 (XIAOMING WANG, GNS)
!     #. UPDATED ON JAN 23 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
	  REAL Z(NX,NY),X(NX),Y(NY)
      INTEGER NN,NX,NY
	  REAL R_MAX,X_REL,M,COEF


	  N1 = 20
	  N2 = 20

	  IF (NX.LE.20) N1 = 5
	  IF (NY.LE.20) N2 = 5

	  M = 1.0
	  COEF = 1.0

	  DO J = 1,NY
		 R_MAX = N1
	  	 !LEFT BOUNDARY
		 DO I = 1,N1
		    X_REL = I-1.0
		    COEF = (EXP((X_REL/R_MAX)**M)-1.0)/(EXP(1.0)-1.0)
		    Z(I,J) = Z(I,J)*COEF
	     ENDDO
	     !RIGHT BOUNDARY
	     DO I = NX-N1,NX
		    X_REL = NX-I
		    COEF = (EXP((X_REL/R_MAX)**M)-1.0)/(EXP(1.0)-1.0)
		    Z(I,J) = Z(I,J)*COEF
	     ENDDO
	  ENDDO

	  DO I = 1,NX
	     R_MAX = N2
	     !TOP BOUNDARY
	     DO J = NY-N2,NY
		    X_REL = NY-J
		    COEF = (EXP((X_REL/R_MAX)**M)-1.0)/(EXP(1.0)-1.0)
		    Z(I,J) = Z(I,J)*COEF
	     ENDDO
	     !BOTTOM BOUNDARY
	     DO J = 1,N2
		    X_REL = J-1.0
		    COEF = (EXP((X_REL/R_MAX)**M)-1.0)/(EXP(1.0)-1.0)
		    Z(I,J) = Z(I,J)*COEF
	     ENDDO
	  ENDDO

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE WRITE_DEFORM (LO,K)
!DESCRIPTION:
!	  #. WRITE THE SEAFLOOR DEFORMATION INTO DATA FILE
!NOTES:
!	  #. CREATED ON DEC.28 2008 (XIAOMING WANG, GNS)
!	  #. UPDATED ON FEB 04 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: LO
      CHARACTER(LEN=40) FNAME
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN
      DATA OSIXTY/0.016666666667/, BIG/-999./

      WRITE (FNAME,1) K
 1    FORMAT('deform_seg',I3.3,'.dat')
      OPEN (25,FILE=FNAME,STATUS='UNKNOWN')
      DO J = 1,LO%NY
         WRITE (25,'(15F8.3)') (LO%DEFORM(I,J),I=1,LO%NX)
      ENDDO
      CLOSE (25)


      RETURN
      END


!*!----------------------------------------------------------------------
!*      SUBROUTINE DEFORM_MAPPING (LO,H,X,Y,NX,NY)
!*!......................................................................
!*!DESCRIPTION:
!*!	  #. MAPPING DEFORMATION DATA ONTO THE 1ST-LEVEL GRIDS
!*!     #. BILINEAR INTERPOLATION WILL BE IMPLEMENTED;
!*!NOTES:
!*!     #. CREATED ON OCT 25 2008 (XIAOMING WANG, GNS)
!*!     #. UPDATED ON NOV 21 2008 (XIAOMING WANG, GNS)
!*!----------------------------------------------------------------------
!*      USE LAYER_PARAMS
!*      TYPE (LAYER) :: LO
!*      INTEGER ISTAT, IS, JS, IE, JE, I0, J0, NX, NY
!*	  REAL H(NX,NY)
!*	  REAL X(NX), Y(NY)
!*	  REAL DEFORM(LO%NX,LO%NY)
!*
!*	  DEFORM = 0.0
!*
!*!      IF (X(1).GE.LO%X(LO%NX) .OR. X(NX).LE.LO%X(1)) THEN
!*!	     PRINT *,"ERROR:: DEFORMATION AREA OUTSIDE THE DOMAIN; EXITING."
!*!         STOP
!*!	  END IF
!*!	  IF (Y(1).GE.LO%Y(LO%NY) .OR. Y(NY).LE.LO%Y(LO%Y(1))) THEN
!*!	     PRINT *,"ERROR:: DEFORMATION AREA OUTSIDE THE DOMAIN; EXITING."
!*!         STOP
!*!	  END IF
!*
!*!	  !PERFORM BILINEAR INTERPOLATION
!*	  DO J = 1,LO%NY
!*	     DO I = 1,LO%NX
!*			KI = 0
!*			KJ = 0
!*		    DO KS = 1,NX-1
!*		       IF (LO%X(I).GE.X(KS) .AND. LO%X(I).LE.X(KS+1)) THEN
!*				  KI = KS
!*			   END IF
!*			END DO
!*		    DO KS = 1,NY-1
!*		       IF (LO%Y(J).GE.Y(KS) .AND. LO%Y(J).LE.Y(KS+1)) THEN
!*				  KJ = KS
!*			   END IF
!*			END DO
!*			IF (KI.GE.1 .AND. KI.LT.NX) THEN
!*			   IF (KJ.GE.1 .AND. KJ.LT.NY) THEN
!*			      DELTA_X = X(KI+1)-X(KI)
!*			      DELTA_Y = Y(KJ+1)-Y(KJ)
!*			      CX = (LO%X(I)-X(KI))/DELTA_X
!*			      CY = (LO%Y(J)-Y(KJ))/DELTA_Y
!*                  Z1 = H(KI,KJ)*(1.0-CX)*(1.0-CY)
!*			      Z2 = H(KI+1,KJ)*(CX)*(1-CY)
!*			      Z3 = H(KI,KJ+1)*(1-CX)*(CY)
!*			      Z4 = H(KI+1,KJ+1)*(CX)*(CY)
!*			      DEFORM(I,J) = Z1+Z2+Z3+Z4
!*			   ENDIF
!*			ENDIF
!*		 END DO
!*	  END DO
!*	  LO%DEFORM(:,:) = DEFORM(:,:)
!*
!*      RETURN
!*      END
