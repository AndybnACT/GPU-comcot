
!----------------------------------------------------------------------
      SUBROUTINE ALL_GRID (LO,LA)
!DESCRIPTION:
!	  #. SOLVE CONTINUITY AND MOMENTUM EQUATIONS FOR ALL SUB-LEVEL GRID
!		 LAYERS;
!INPUT:
!	  #. WATER SURFACE DISPLACEMENT AND VOLUME FLUXES OF LO;
!OUTPUT:
!	  #. WATER SURFACE DISPLACEMENT AND VOLUME FLUXES OF LA;
!NOTES:
!     #. NEWQ_S,NEWQ_C COMBINED INTO NEWQ (JUL 2003, XIAOMING WANG)
!     #. UPDATED ON SEP 17 2006 (XIAOMING WANG)
!     #. LEVEL2%DT = LEVEL1%DT/LEVEL2%REL_SIZE
!     #. UPDATED NOV.21 2008 (XIAOMING WANG, GNS)
!        1. UP TO 12 LEVELS OF SUB GRID LAYERS IMPLEMENTED
!     #. UPDATED DEC 22 2008 (XIAOMING WANG, GNS)
!		 1. TIME STEP SIZE RATIO IS NO LONGER FIXED AS 2
!			BUT DETERMINED FROM WATER DEPTH OF EACH GRID LAYER;
!        2. TIME STEP SIZE RATIO OF LO TO LA COULD BE ANY INTEGER;
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER) 	:: LO
	  TYPE (LAYER),DIMENSION(NUM_GRID)  :: LA
	  INTEGER K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12
	  INTEGER L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12
	  INTEGER NHALF
!	  LOGICAL UPZ(NUM_GRID)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!<<...START OF LEVEL 2
	  DO L2=1,NUM_GRID  
	     IF (LA(L2)%LAYSWITCH.EQ.0 .AND. LA(L2)%PARENT.EQ.LO%ID) THEN
	        DO K2 = 1,LA(L2)%REL_TIME
			   IF (K2 .EQ. 1) THEN
			      CALL JNQ (LO,LA(L2))
			   ELSE
			      CALL NEWQ (LO,LA(L2),K2)
			   ENDIF
               IF (LA(L2)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L2)%LAYCORD.NE.0) CALL SED_TRANSPORT (LA(L2))
			   CALL MASS (LA(L2))

!<<............START OF LEVEL 3
		       DO L3=1,NUM_GRID
	           IF (LA(L3)%LAYSWITCH.EQ.0 .AND.						&
									LA(L3)%PARENT.EQ.LA(L2)%ID) THEN
			   DO K3 = 1,LA(L3)%REL_TIME
			   IF (K3 .EQ. 1) THEN
				  CALL JNQ (LA(L2),LA(L3))
			   ELSE
				  CALL NEWQ (LA(L2),LA(L3),K3)
			   ENDIF
			   IF (LA(L3)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L3)%LAYCORD.NE.0) CALL SED_TRANSPORT (LA(L3))
			   CALL MASS (LA(L3))

!<<............START OF LEVEL 4 
               DO L4=1,NUM_GRID
               IF (LA(L4)%LAYSWITCH.EQ.0 .AND.						&
									LA(L4)%PARENT.EQ.LA(L3)%ID) THEN
			   DO K4 = 1,LA(L4)%REL_TIME
			   IF (K4 .EQ. 1) THEN
				  CALL JNQ (LA(L3),LA(L4))
			   ELSE
				  CALL NEWQ (LA(L3),LA(L4),K4)
			   ENDIF
			   IF (LA(L4)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L4)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L4))
			   CALL MASS (LA(L4))

!<<........... START OF LEVEL 5 
               DO L5=1,NUM_GRID
               IF (LA(L5)%LAYSWITCH.EQ.0 .AND.						&
									LA(L5)%PARENT.EQ.LA(L4)%ID) THEN
			   DO K5 = 1,LA(L5)%REL_TIME
			   IF (K5 .EQ. 1) THEN
				  CALL JNQ (LA(L4),LA(L5))
			   ELSE
			      CALL NEWQ (LA(L4),LA(L5),K5)
			   ENDIF
			   IF (LA(L5)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L5)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L5))
			   CALL MASS (LA(L5))

!<<............START OF LEVEL 6
               DO L6=1,NUM_GRID
               IF (LA(L6)%LAYSWITCH.EQ.0 .AND.						&
									LA(L6)%PARENT.EQ.LA(L5)%ID) THEN
			   DO K6 = 1,LA(L6)%REL_TIME
			   IF (K6 .EQ. 1) THEN
			      CALL JNQ (LA(L5),LA(L6))
			   ELSE
			      CALL NEWQ (LA(L5),LA(L6),K6)
			   ENDIF
			   IF (LA(L6)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L6)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L6))
			   CALL MASS (LA(L6))

!<<........... START OF LEVEL 7
               DO L7=1,NUM_GRID
               IF (LA(L7)%LAYSWITCH.EQ.0 .AND.						&
									LA(L7)%PARENT.EQ.LA(L6)%ID) THEN
			   DO K7 = 1,LA(L7)%REL_TIME
			   IF (K7 .EQ. 1) THEN
			      CALL JNQ (LA(L6),LA(L7))
			   ELSE
			      CALL NEWQ (LA(L6),LA(L7),K7)
			   ENDIF
			   IF (LA(L7)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L7)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L7))
			   CALL MASS (LA(L7))

!<<........... START OF LEVEL 8
               DO L8=1,NUM_GRID
               IF (LA(L8)%LAYSWITCH.EQ.0 .AND.						&
									LA(L8)%PARENT.EQ.LA(L7)%ID) THEN
			   DO K8 = 1,LA(L8)%REL_TIME
			   IF (K8 .EQ. 1) THEN
				  CALL JNQ (LA(L7),LA(L8))
			   ELSE
			      CALL NEWQ (LA(L7),LA(L8),K8)
			   ENDIF
			   IF (LA(L8)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L8)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L8))
			   CALL MASS (LA(L8))

!<<............START OF LEVEL 9
               DO L9=1,NUM_GRID
               IF (LA(L9)%LAYSWITCH.EQ.0 .AND.						&
									LA(L9)%PARENT.EQ.LA(L8)%ID) THEN
			   DO K9 = 1,LA(L9)%REL_TIME
			   IF (K9 .EQ. 1) THEN
				  CALL JNQ (LA(L8),LA(L9))
			   ELSE
				  CALL NEWQ (LA(L8),LA(L9),K9)
			   ENDIF
			   IF (LA(L9)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L9)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L9))
			   CALL MASS (LA(L9))

!<<............START OF LEVEL 10
               DO L10=1,NUM_GRID
               IF (LA(L10)%LAYSWITCH.EQ.0 .AND.						&
									LA(L10)%PARENT.EQ.LA(L9)%ID) THEN
			   DO K10 = 1,LA(L10)%REL_TIME
			   IF (K10 .EQ. 1) THEN
			      CALL JNQ (LA(L9),LA(L10))
			   ELSE
				  CALL NEWQ (LA(L9),LA(L10),K10)
			   ENDIF
			   IF (LA(L10)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L10)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L10))
			   CALL MASS (LA(L10))

!<<............START OF LEVEL 11
               DO L11=1,NUM_GRID
               IF (LA(L11)%LAYSWITCH.EQ.0 .AND.						&
								LA(L11)%PARENT.EQ.LA(L10)%ID) THEN
			   DO K11 = 1,LA(L11)%REL_TIME
			   IF (K11 .EQ. 1) THEN
				  CALL JNQ (LA(L10),LA(L11))
			   ELSE
				  CALL NEWQ (LA(L10),LA(L11),K11)
			   ENDIF
			   IF (LA(L11)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L11)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L11))
			   CALL MASS (LA(L11))

!<<........... START OF LEVEL 12
               DO L12=1,NUM_GRID
               IF (LA(L12)%LAYSWITCH.EQ.0 .AND.						&
								LA(L12)%PARENT.EQ.LA(L11)%ID) THEN
			   DO K12 = 1,LA(L12)%REL_TIME
			   IF (K12 .EQ. 1) THEN
			      CALL JNQ (LA(L11),LA(L12))
			   ELSE
				  CALL NEWQ (LA(L11),LA(L12),K12)
			   ENDIF
			   IF (LA(L12)%SEDI_SWITCH.EQ.0 .AND.					&
					LA(L12)%LAYCORD.NE.0 ) CALL SED_TRANSPORT (LA(L12))
			   CALL MASS (LA(L12))

!<<........... INSERT MORE LEVELS HERE

!>>........... END OF INSERTING MORE LEVLES

			   CALL MOMENT (LA(L12))
			   NHALF = FLOOR(LA(L12)%REL_TIME/2.0)+1
			   IF (K12.EQ.NHALF) THEN
				  CALL JNZ (LA(L11),LA(L12))
			   ENDIF
			   CALL UPDATE (LA(L12))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 12

			   CALL MOMENT (LA(L11))
			   NHALF = FLOOR(LA(L11)%REL_TIME/2.0)+1
			   IF (K11.EQ.NHALF) THEN
				  CALL JNZ (LA(L10),LA(L11))
			   ENDIF
			   CALL UPDATE (LA(L11))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 11

			   CALL MOMENT (LA(L10))
			   NHALF = FLOOR(LA(L10)%REL_TIME/2.0)+1
			   IF (K10 .EQ. NHALF) THEN
				  CALL JNZ (LA(L9),LA(L10))
			   ENDIF
			   CALL UPDATE (LA(L10))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 10

			   CALL MOMENT (LA(L9))
			   NHALF = FLOOR(LA(L9)%REL_TIME/2.0)+1
			   IF (K9 .EQ. NHALF) THEN
				  CALL JNZ (LA(L8),LA(L9))
			   ENDIF
			   CALL UPDATE (LA(L9))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 9

			   CALL MOMENT (LA(L8))
			   NHALF = FLOOR(LA(L8)%REL_TIME/2.0)+1
			   IF (K8 .EQ. NHALF) THEN
				  CALL JNZ (LA(L7),LA(L8))
			   ENDIF
			   CALL UPDATE (LA(L8))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 8

			   CALL MOMENT (LA(L7))
			   NHALF = FLOOR(LA(L7)%REL_TIME/2.0)+1
			   IF (K7 .EQ. NHALF) THEN
				  CALL JNZ (LA(L6),LA(L7))
			   ENDIF
			   CALL UPDATE (LA(L7))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 7

			   CALL MOMENT (LA(L6))
			   NHALF = FLOOR(LA(L6)%REL_TIME/2.0)+1
			   IF (K6 .EQ. NHALF) THEN
				  CALL JNZ (LA(L5),LA(L6))
			   ENDIF
			   CALL UPDATE (LA(L6))
			   ENDDO
			   ENDIF
			   ENDDO
!>>........... END OF LEVEL 6

			   CALL MOMENT (LA(L5))
			   NHALF = FLOOR(LA(L5)%REL_TIME/2.0)+1
			   IF (K5 .EQ. NHALF) THEN
				  CALL JNZ (LA(L4),LA(L5))
			   ENDIF
			   CALL UPDATE (LA(L5))
			   ENDDO
			   ENDIF
			   ENDDO
!<<........... END OF LEVEL 5

			   CALL MOMENT (LA(L4))
			   NHALF = FLOOR(LA(L4)%REL_TIME/2.0)+1
			   IF (K4 .EQ. NHALF) THEN
				  CALL JNZ (LA(L3),LA(L4))
			   ENDIF
			   CALL UPDATE (LA(L4))
			   ENDDO
			   ENDIF
			   ENDDO
!<<........... END OF LEVEL 4

               CALL MOMENT (LA(L3))
			   NHALF = FLOOR(LA(L3)%REL_TIME/2.0)+1
			   IF (K3 .EQ. NHALF) THEN
				  CALL JNZ (LA(L2),LA(L3))
			   ENDIF
               CALL UPDATE (LA(L3))
			   ENDDO
			   ENDIF
			   ENDDO
!>>............END OF LEVEL 3

               CALL MOMENT (LA(L2))
			   NHALF = FLOOR(LA(L2)%REL_TIME/2.0)+1
			   IF (K2 .EQ. NHALF) THEN
				  CALL JNZ (LO,LA(L2))
			   ENDIF
               CALL UPDATE (LA(L2))
			ENDDO
		 ENDIF
      ENDDO
!>>...END OF LEVEL 2

      RETURN
      END


!-----------------------------------------------------------------------
      SUBROUTINE UPDATE (LO)
!.....TRANSFER INFORMATION FROM LAST STEP TO NEXT STEP (INNER LAYER) 
!     UPDATED ON SEP 17 2006 BY XIAOMING WANG
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
	  TYPE (LAYER)  :: LO
 
      LO%Z(1,:,2) = LO%Z(2,:,2)
      LO%Z(:,1,2) = LO%Z(:,2,2)
    
      LO%Z(:,:,1) = LO%Z(:,:,2)
      LO%M(:,:,1) = LO%M(:,:,2)
      LO%N(:,:,1) = LO%N(:,:,2)
	  IF (LO%LAYGOV.GT.1) THEN
	     LO%M0(:,:)  = LO%M(:,:,1)
	  	 LO%N0(:,:)  = LO%N(:,:,1)
	  ENDIF
	  IF (LO%INI_SWITCH.EQ.3 .OR. LO%INI_SWITCH.EQ.4) THEN
	     LO%HT(:,:,1) = LO%HT(:,:,2)
	     LO%H(:,:) = LO%H(:,:) + LO%HT(:,:,2) - LO%HT(:,:,1)
	  ENDIF
      
      RETURN
      END
      
!----------------------------------------------------------------------
      SUBROUTINE JNZ (LO,LA)
!......................................................................
!DESCRIPTION:
!     #. THIS SUBROUTINE IS USED TO UPDATE FREE SURFACE ELEVATION 
!		 OF LO (LARGER GRIDS) WITH THAT OF LA (SMALLER GRIDS)
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTE:
!	  #. REVISED SIGNIFICANTLY ON NOV 2004 (XIAOMING WANG,CORNELL)
!     #. UPDATED ON JAN05 2009 (XIAOMING WANG, GNS) 
!	  #. UPDATED ON APR03 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA
	  INTEGER   IS,IE,JS,JE,IR
	  REAL HALF, SUM, COUNT
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  IF (LA%SC_OPTION .EQ. 0) THEN
      HALF = (LA%REL_SIZE*LA%REL_SIZE)/2.0
      IR=LA%REL_SIZE
      IS = LA%CORNERS(1)
      IE = LA%CORNERS(2)
      JS = LA%CORNERS(3)
      JE = LA%CORNERS(4)

	  I_SHIFT = -IS*IR+1
	  J_SHIFT = -JS*IR+1
      DO I = IS,IE
	     DO J = JS,JE
	        SUM = 0.0
		    COUNT = 0.0
			I0 = I*IR+I_SHIFT
			J0 = J*IR+J_SHIFT
			DO KI = 1,LA%REL_SIZE
			   DO KJ = 1,LA%REL_SIZE
!			      II=(I-IS)*IR+1+KI
!			      JJ=(J-JS)*IR+1+KJ
			      II = I0 + KI		
				  JJ = J0 + KJ		
!*			      IF (LA%H(II,JJ) .GT. GX) THEN
			      IF (LA%H(II,JJ)+LA%Z(II,JJ,2) .GT. GX) THEN
					 IF (MOD(LA%REL_TIME,2) .EQ. 0) THEN
						SUM = SUM+0.5*(LA%Z(II,JJ,1)+LA%Z(II,JJ,2))
					 ELSE
					    SUM = SUM+LA%Z(II,JJ,2)
					 ENDIF
!*			         SUM = SUM+LA%Z(IS,JS,2)
				     COUNT = COUNT + 1.0
			      ENDIF
			   ENDDO
		    ENDDO
            IF (COUNT .GT. HALF) THEN
               LO%Z(I,J,2) = SUM/COUNT
            ELSE
               LO%Z(I,J,2) = 0.0
            ENDIF
         ENDDO
      ENDDO
	  ELSE
	  CALL JNZ_SC (LO,LA)
	  ENDIF

      RETURN
      END

      
!-----------------------------------------------------------------------
      SUBROUTINE JNQ (LO,LA)
!DESCRIPTION:
!	  #. INTERPOLATE VOLUME FLUXES (HU,HV) FROM OUTER GRID REGION 
!		 (PARENT GRID LAYER, LO) INTO INNER/NESTED GRID LAYER
!		 (CHILD GRID LAYER, LA) ALONG CONNECTING BOUNDARIES AT 
!		  THE BEGINNING OF A TIME STEP;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTES:
!     #. REVISED SIGNIFICANTLY ON NOV 2004 (XIAOMING WANG, CORNELL)
!	  #. UPDATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR03 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA

	  IF (LA%SC_OPTION.EQ.0) THEN
	  IS = LA%CORNERS(1)
      IE = LA%CORNERS(2)
	  JS = LA%CORNERS(3)
      JE = LA%CORNERS(4)
      CALL EDGEINTERP_VERT (LO%M(IS-1,:,1),LO%NY,LA,0) !LEFT BOUNDARY
	  CALL EDGEINTERP_VERT (LO%M(IE,:,1),LO%NY,LA,1)   !RIGHT BOUNDARY
	  CALL EDGEINTERP_HORI (LO%N(:,JS-1,1),LO%NX,LA,0) !BOTTOM BOUNDARY
	  CALL EDGEINTERP_HORI (LO%N(:,JE,1),LO%NX,LA,1)   !TOP BOUNDARY
	  ELSE
	  CALL EDGE_INTERP_SC (LO,LA)
	  ENDIF

      RETURN
      END


!-----------------------------------------------------------------------
      SUBROUTINE EDGEINTERP_VERT (FLUX,NY,LA,SIDE)
!.....ADDED BY XIAOMING WANG (JAN 22, 2006)
!.....SIGNIFICANT CHANGE FROM OLD SUB. (BEFORE NOV 2005): 
!	  INTERPOLATE FLUX INSTEAD OF VELOCITY
!.....INTERPOLATING VOLUMN FLUX ALONG HORIZONTAL CONNECTING BOUNDARY
!     FROM COARSE GRIDS INTO FINER GRIDS
!         1 -- TOP CONNECTING BOUNDARY
!.....SIDE:
!         0 -- LEFT CONNECTING BOUNDARY
!         1 -- RIGHT CONNECTING BOUNDARY
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LA
	  INTEGER           :: NY,SIDE,NN,IR,XS,XE,YS,YE
!	  REAL, DIMENSION(2*LA%REL_SIZE+1)   :: EDGE
	  REAL, DIMENSION(NY)   :: FLUX
	  REAL C1,C2,FLUX_X,FLUX_E
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!	  DEGE = 0.0
      
	  NN = 2*LA%REL_SIZE
	  IR = LA%REL_SIZE
	  XS = LA%CORNERS(1)
      XE = LA%CORNERS(2)
	  YS = LA%CORNERS(3)
      YE = LA%CORNERS(4)

	  J_SHIFT = -YS*IR+1
	  IF (SIDE .EQ. 0) II = 1      !FOR LEFT BOUNDARY, SIDE=0
	  IF (SIDE .EQ. 1) II = LA%NX  !FOR RIGHT BOUNDARY, SIDE=1

      DO J=YS,YE
	     FLUX_S = 0.5*(FLUX(J-1)+FLUX(J))
		 FLUX_E = 0.5*(FLUX(J+1)+FLUX(J))
		 C1 = (FLUX_E-FLUX_S)/DBLE(NN)
		 C2 = -C1 + FLUX_S
!*		 DO K=1,NN+1
!*		    EDGE(K) = K*C1 + C2 !DBLE(K-1)*(FLUX_E-FLUX_S)/DBLE(NN)+FLUX_S
!*	     ENDDO
		 JS = J*IR + J_SHIFT !(J-YS)*IR+1
		 DO K=1,LA%REL_SIZE
			JJ = JS+K
		    LA%M(II,JJ,1) = 2.0*K*C1 + C2	!EDGE(2*K)
!*			IF (LA%H(II,JJ) .LE. GX) LA%M(II,JJ,1) = 0.0
			IF (LA%H(II,JJ)+LA%Z(II,JJ,1) .LE. GX) LA%M(II,JJ,1) = 0.0
	     ENDDO
      ENDDO

      RETURN
	  END

!-----------------------------------------------------------------------
      SUBROUTINE EDGEINTERP_HORI (FLUX,NX,LA,SIDE)
!.....ADDED BY XIAOMING WANG (JAN 22, 2006)
!.....SIGNIFICANT CHANGE FROM OLD SUB. (BEFORE NOV 2005): 
!	  INTERPOLATE FLUX INSTEAD OF VELOCITY
!.....INTERPOLATING VOLUMN FLUX ALONG HORIZONTAL CONNECTING BOUNDARY
!     FROM COARSE GRIDS INTO FINER GRIDS
!.....SIDE:
!         0 -- BOTTOM CONNECTING BOUNDARY
!         1 -- TOP CONNECTING BOUNDARY
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA
	  INTEGER           :: SIDE,XS,XE,YS,YE
!	  REAL, DIMENSION(2*LA%REL_SIZE+1)   :: EDGE
	  REAL, DIMENSION(NX)  :: FLUX
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN 
!	  EDGE = 0.0
 
	  NN = 2*LA%REL_SIZE
	  IR = LA%REL_SIZE
	  XS = LA%CORNERS(1)
      XE = LA%CORNERS(2)
	  YS = LA%CORNERS(3)
      YE = LA%CORNERS(4)

	  I_SHIFT = -XS*IR+1
	  IF (SIDE .EQ. 0) JJ = 1      !FOR BOTTOM BOUNDARY, SIDE=0
	  IF (SIDE .EQ. 1) JJ = LA%ny  !FOR TOP BOUNDARY, SIDE=1

      DO I=XS,XE
	     FLUX_S = 0.5*(FLUX(I-1)+FLUX(I))
		 FLUX_E = 0.5*(FLUX(I+1)+FLUX(I))
		 C1 = (FLUX_E-FLUX_S)/DBLE(NN)
		 C2 = -C1+FLUX_S
!*		 DO K=1,NN+1
!*		    EDGE(K) = K*C1+C2 !DBLE(K-1)*(FLUX_E-FLUX_S)/DBLE(NN)+FLUX_S
!*	     ENDDO
		 IS = I*IR+I_SHIFT	  !(I-XS)*IR+1
		 DO K=1,LA%REL_SIZE
			II = IS+K
		    LA%N(II,JJ,1) = 2.0*K*C1+C2	!EDGE(2*K)
!*			IF (LA%H(II,JJ) .LE. 0.0) LA%N(II,JJ,1) = 0.0
			IF (LA%H(II,JJ)+LA%Z(II,JJ,1) .LE. GX) LA%N(II,JJ,1) = 0.0
	     ENDDO
      ENDDO

      RETURN
	  END


!-----------------------------------------------------------------------
      SUBROUTINE NEWQ (LO, LA, T)
!.......................................................................
!DESCRIPTION:
!	  #. INTERPOLATING VOLUME FLUX FROM OUTER LAYER INTO INNER LAYER
!		 THROUGH FOUR CONNECTING BOUNDARIES (I.E., FOUR BOUNDARIES OF 
!		 GRID LAYER LA) WHICH SERVE AS THE FLUX  BOUNDARY CONDITION OF LA
!		 AT THE MOMENT BETWEEN N*DT AND (N+1)*DT WHEN TIME STEP SIZE 
!		 RATIO OF LO TO LA IS LARGER THAN 1;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!INPUT:
!	  #. LO: PARENT GRID INFO (OUTER GRID)
!	  #. LA: CHILD GRID INFO (INNER GRID)
!	  #. T: AN INTEGER AMONG 1 TO LA%REL_TIME
!NOTES:
!	  #. CREATED BY XIAOMING WANG (CORNELL, 2005)
!	  #. UPDATED ON JAN 25 2006 (XIAOMING WANG, CORNELL UNIVERSITY)
!	  #. UPDATED ON JAN09 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR03 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA
	  INTEGER T

	  IF (LA%SC_OPTION .EQ. 0) THEN
      CALL NEWQ_VERT(LO,LA,T,0) ! LEFT BOUNDARY
	  CALL NEWQ_VERT(LO,LA,T,1) ! RIGHT BOUNDARY
	  CALL NEWQ_HORI(LO,LA,T,0) ! BOTTOM BOUNDARY
	  CALL NEWQ_HORI(LO,LA,T,1) ! TOP BOUNDARY
	  ELSE
	  CALL NEWQ_SC (LO, LA, T)
	  ENDIF

	  RETURN
	  END


!-----------------------------------------------------------------------
      SUBROUTINE NEWQ_VERT (LO, LA, T, SIDE)
!.......................................................................
!DESCRIPTION:
!	  #. INTERPOLATING VOLUME FLUX FROM OUTER LAYER INTO INNER LAYER
!		 THROUGH LEFT AND RIGHT CONNECTING BOUNDARIES (I.E., LEFT AND 
!		 RIGHT BOUNDARIES OF GRID LAYER LA) WHICH SERVE AS THE FLUX
!		 BOUNDARY CONDITION OF LA AT THE MOMENT BETWEEN N*DT AND (N+1)*DT
!		 WHEN TIME STEP SIZE RATIO OF LO TO LA IS LARGER THAN 1;
!     #. SIDE:
!			0 -- LEFT CONNECTING BOUNDARY
!			1 -- RIGHT CONNECTING BOUNDARY
!NOTES:
!	  #. CREATED BY XIAOMING WANG (CORNELL, 2006)
!	  #. UPDATED ON JAN09 2009 (XIAOMING WANG, GNS) 
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA
	  INTEGER SIDE
	  INTEGER T
	  REAL HM, XM, GRX
	  REAL YFLUX(LO%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  GRX = LO%GRX
	  YFLUX = 0.0
	  XM = 0.0

      IS = LA%CORNERS(1)
      IE = LA%CORNERS(2)
      JS = LA%CORNERS(3)
      JE = LA%CORNERS(4)

	  IF (SIDE .EQ. 0) I = IS-1 !FOR LEFT BOUNDARY, SIDE=0
	  IF (SIDE .EQ. 1) I = IE   !FOR RIGHT BOUNDARY, SIDE=1

	  C1 = (T-1.0)/LA%REL_TIME
	  C2 = 1.0-C1

      IP1 = I+1
      DO J=JS-2,JE+2
         IF ((LO%H(I,J)+LO%Z(I,J,2).GT.GX)							&
					.AND. (LO%H(IP1,J)+LO%Z(IP1,J,2).GT.GX)) THEN
			IF (T.EQ.2) THEN !THIS ENSURES FORECAST CALCULATION ONLY BEING DONE ONCE
		       IF (LO%LAYCORD .EQ. 0) THEN
		          JM1 = J-1
                  TOT_N = LO%N(I,J,1)+LO%N(IP1,J,1)+LO%N(I,JM1,1)	&
													+LO%N(IP1,JM1,1)
                  XM = LO%M(I,J,1)-LO%R2(I,J)*(LO%Z(IP1,J,2)		&
									-LO%Z(I,J,2))+LO%R3(I,J)*TOT_N
			   ELSE
				  HM = LO%HP(I,J)+0.5*(LO%Z(I,J,2)+LO%Z(IP1,J,2))
                  XM = LO%M(I,J,1)-GRX*HM*(LO%Z(IP1,J,2)-LO%Z(I,J,2))
		       ENDIF
               IF (ABS(XM).LT.EPS) XM = ZERO
!*             LO%YFLUX(J) = 0.5*(XM+LO%M(I,J,1))
!              LO%YFLUX(J) = (T-1.0)/LA%REL_TIME*(XM-LO%M(I,J,1))+LO%M(I,J,1)
			   LO%YFLUX(J,SIDE+1) = XM
			ENDIF
		    YFLUX(J) = C1*LO%YFLUX(J,SIDE+1) + C2*LO%M(I,J,1)
         END IF
      END DO
      CALL EDGEINTERP_VERT(YFLUX,LO%NY,LA,SIDE)

      RETURN
	  END

!-----------------------------------------------------------------------
      SUBROUTINE NEWQ_HORI (LO, LA, T, SIDE)
!.......................................................................
!DESCRIPTION:
!	  #. INTERPOLATING VOLUME FLUX FROM OUTER LAYER INTO INNER LAYER
!		 THROUGH TOP AND BOTTOM CONNECTING BOUNDARIES (I.E., TOP AND 
!		 BOTTOM BOUNDARIES OF GRID LAYER LA) WHICH SERVE AS THE FLUX
!		 BOUNDARY CONDITION OF LA AT THE MOMENT BETWEEN N*DT AND (N+1)*DT
!		 WHEN TIME STEP SIZE RATIO OF LO TO LA IS LARGER THAN 1;
!     #. SIDE:
!			0 -- BOTTOM CONNECTING BOUNDARY
!			1 -- TOP CONNECTING BOUNDARY
!NOTES:
!	  #. CREATED BY XIAOMING WANG (CORNELL, 2006)
!	  #. UPDATED ON JAN09 2009 (XIAOMING WANG, GNS) 
!-----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)		:: LO, LA
	  INTEGER SIDE
	  INTEGER T
	  REAL HN, XN, GRY
	  REAL XFLUX(LO%NX)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  GRY = LO%GRY
	  XFLUX = 0.0
	  XN = 0.0

      IS = LA%CORNERS(1)
      IE = LA%CORNERS(2)
      JS = LA%CORNERS(3)
      JE = LA%CORNERS(4)

	  IF (SIDE .EQ. 0) J = JS-1 !FOR BOTTOM BOUNDARY, SIDE=0
	  IF (SIDE .EQ. 1) J = JE   !FOR TOP BOUNDARY, SIDE=1

	  C1 = (T-1.0)/LA%REL_TIME
	  C2 = 1-C1

      JP1 = J+1
      DO I=IS-2,IE+2
         IF ((LO%H(I,J)+LO%Z(I,J,2).GT.GX)							&
					.AND. (LO%H(I,JP1)+LO%Z(I,JP1,2).GT.GX)) THEN
			! THIS ENSURES FORECAST CALCULATION ONLY BEING DONE ONCE
	        IF (T.EQ.2) THEN 
		       IM1 = I-1
		       IF (LO%LAYCORD .EQ. 0) THEN
                  TOT_M = LO%M(IM1,J,1)+LO%M(IM1,JP1,1)+LO%M(I,J,1)	&
													+LO%M(I,JP1,1)
                  XN = LO%N(I,J,1)-LO%R4(I,J)*(LO%Z(I,JP1,2)		&
									-LO%Z(I,J,2))-LO%R5(I,J)*TOT_M	
			   ELSE
				  HN = LO%HQ(I,J)+0.5*(LO%Z(I,J,2)+LO%Z(I,JP1,2))
                  XN = LO%N(I,J,1)-GRY*HN*(LO%Z(I,JP1,2)-LO%Z(I,J,2))
		       ENDIF
               IF (ABS(XN).LT.EPS) XN = ZERO
			   LO%XFLUX(I,SIDE+1) = XN
!*             LO%XFLUX(I) = 0.5*(XN+LO%N(I,J,1))
!              LO%XFLUX(I) = (T-1.0)/LA%REL_TIME*(XN-LO%N(I,J,1))+LO%N(I,J,1)
	        ENDIF
		    XFLUX(I) = C1*LO%XFLUX(I,SIDE+1) + C2*LO%N(I,J,1)
         END IF
      END DO

      CALL EDGEINTERP_HORI(XFLUX,LO%NX,LA,SIDE)

      RETURN
	  END


!----------------------------------------------------------------------
	  SUBROUTINE ININTERP (LO,LA)
!DESCRIPTION:	
!	  #. INTERPOLATE DEFORMATION PROFILE FROM 1ST-LEVEL GRID ALYER 
!	     LO INTO ALL SUB-LEVEL GRID LAYERS, LA;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTES:
!	  #. UPDATED ON OCT 28 2004 (XIAOMING WANG)
!	  #. UPDATED ON SEP 17 2006 (XIAOMING WANG)
!	  #. LAST REVISE: JAN.05 2009 (XIAOMING WANG)
!	  #. UPDATED ON FEB03 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR03 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO
	  TYPE (LAYER),DIMENSION(NUM_GRID)  :: LA
	  INTEGER OPTION
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!.....INTERPOLATE FROM 1ST-LEVEL GRID REGION LO TO 2ND-LEVEL REGION LA
	  DO I = 1,NUM_GRID
	     IF (LA(I)%LAYSWITCH.EQ.0 .AND. LA(I)%PARENT.EQ.LO%ID) THEN
			IF (LA(I)%SC_OPTION .EQ. 0) THEN
		    CALL CELL_INTERP (LO,LA(I))
			ELSE
			CALL CELL_INTERP_SC (LO,LA(I))
			ENDIF
		 ENDIF
	  ENDDO
!.....INTERPOLATE FROM 2ND-LEVEL GRID REGION TO ALL SUB-LEVEL GRIDS      
	  DO L=2,NUM_GRID+1
		 DO J=1,NUM_GRID
		    IF (LA(J)%LAYSWITCH.EQ.0 .AND. LA(J)%LEVEL.EQ.L) THEN
               DO K=1,NUM_GRID
				  IF (LA(K)%LAYSWITCH.EQ.0							&
								.AND. LA(K)%PARENT.EQ.LA(J)%ID) THEN
					 IF (LA(K)%SC_OPTION .EQ. 0) THEN
					 CALL CELL_INTERP (LA(J),LA(K)) 
					 ELSE
					 CALL CELL_INTERP_SC (LA(J),LA(K))
					 ENDIF
				  ENDIF
			   ENDDO
		    ENDIF
		 ENDDO   
	  ENDDO	   

      RETURN
      END


!----------------------------------------------------------------------
      SUBROUTINE CELL_INTERP (LO,LA)
!......................................................................
!DESCRIPTION:
!	  #. INTERPOLATE INITIAL SURFACE DEFORMATION FROM OUTER REGION (LO)
!	     INTO INNER REGIONS (LA)
!NOTE: 
!	  #. CREATED ON OCT 28, 2004 (XIAOMING WANG, CORNELL UNIVERSITY)
!        ORIGINAL SUBROUTINE INI_TPL IS REMOVED TO 
!	     SIMPLIFY VARIABLE STRUCTURE
!	  #. UPDATED ON JAN06 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO, LA
      INTEGER       :: IS,JS,IE,JE,NX,IR
	  REAL TL,TR,BL,BR
	  REAL CELL(2*LA%REL_SIZE+1,2*LA%REL_SIZE+1)
	  REAL SMALL(LA%REL_SIZE,LA%REL_SIZE)
	  REAL DEFORM (LA%NX,LA%NY)

	  SMALL = 0.0
	  CELL = 0.0
	  DEFORM = 0.0

	  IS = LA%CORNERS(1)
	  IE = LA%CORNERS(2)
	  JS = LA%CORNERS(3)
	  JE = LA%CORNERS(4)

      NX = 2*LA%REL_SIZE
	  IR = LA%REL_SIZE
	  DO I = IS, IE
	     IP1 = I + 1
		 IM1 = I - 1
		 DO J = JS, JE
			JP1 = J + 1
			JM1 = J - 1
		    !GET DEFORMATION AT FOUR CORNERS OF A LARGE GRID CELL
            TL = (LO%DEFORM(I,J)+LO%DEFORM(I,JP1)+LO%DEFORM(IM1,J)	&
						+LO%DEFORM(IM1,JP1))/4.0 !TOP LEFT CORNER
            TR = (LO%DEFORM(I,J)+LO%DEFORM(I,JP1)+LO%DEFORM(IP1,J)	&
						+LO%DEFORM(IP1,JP1))/4.0 !TOP RIGHT CORNER
            BL = (LO%DEFORM(I,J)+LO%DEFORM(I,JM1)+LO%DEFORM(IM1,J)	&
						+LO%DEFORM(IM1,JM1))/4.0 !BOTTOM LEFT CORNER
            BR = (LO%DEFORM(I,J)+LO%DEFORM(I,JM1)+LO%DEFORM(IP1,J)	&
						+LO%DEFORM(IP1,JM1))/4.0 !BOTTOM RIGHT CORNER
!...........CELL STORES WATERDEPTH VALUES CENTERED AT LARGE GRID (I,J)
!			WHICH ARE INTERPOLATED FROM ITS ADJACENT LARGE GRIDS
!...........GET VALUES ALONG TWO VERTICAL BOUNDARY OF THE LARGE CELL
            DO K=1,NX+1
               CELL(1,K) = (K-1.0)*(TL-BL)/NX + BL
               CELL(NX+1,K) = (K-1.0)*(TR-BR)/NX + BR
            ENDDO
            !INTERPOLATED FROM TWO VERTICAL BOUNDARIES
            DO M=1, NX+1
               DO K=1, NX+1
                  CELL(M,K)=(M-1.)*(CELL(NX+1,K)-CELL(1,K))/NX+CELL(1,K)
               ENDDO
            ENDDO
!...........GET INTERPOLATED WATER DEPTH VALUES FOR NESTED GRIDS 
!			OVERLAPPED BY CELL (LAY1(I,J))
            DO K=1,LA%REL_SIZE
               DO M=1,LA%REL_SIZE
                  SMALL(K,M)=CELL(2*K,2*M)
               ENDDO
            ENDDO
            DEFORM((I-IS)*IR+2:(I-IS+1)*IR+1,						&
						(J-JS)*IR+2:(J-JS+1)*IR+1) = SMALL(:,:)
         ENDDO
      ENDDO 
	  LA%DEFORM(:,:) = DEFORM(:,:) 
!	  LA%Z(:,:,1) = LA%Z(:,:,1) + LA%DEFORM(:,:)

!*	  LA%Z((I-XS)*IR+2:(I-XS+1)*IR+1,(J-YS)*IR+2:(J-YS+1)*IR+1,1) = SMALL(:,:)
	   
	  RETURN 
	  END


!----------------------------------------------------------------------
      SUBROUTINE CELL_INTERP_SC (LO,LA)
!......................................................................
!DESCRIPTION:
!	  #. INTERPOLATE INITIAL SURFACE DEFORMATION FROM OUTER REGION (LO)
!	     INTO INNER REGIONS (LA)
!     #. DESIGNED FOR INTERPOLATION FROM SPHERICAL GRID LAYER, LO, TO 
!	     CARTESIAN GRID LAYERS, LA;
!	  #. REVERSE UNIVERSAL TRANSVERSE MERCATOR PROJECTION IS ADOPTED 
!		 TO PROJECT GRIDS IN CARTESIAN COORDINATE SYSTEM ONTO THE 
!		 EARTH SPHERE SURFACE (WGS84 ELLIPSOID);
!	  #. LA%UPZ = 
!			.TRUE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT 
!					THE SAME COORDINATE SYSTEM;
!			.FALSE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT  
!					DIFFERENT COORDINATE SYSTEM;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTE: 
!	  #. CREATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR06 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO, LA
      INTEGER KI,KJ,I,J
	  REAL Z1,Z2,Z3,Z4
	  REAL DEFORM(LA%NX,LA%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  DEFORM = 0.0

	  DO I = 1,LA%NX
		 DO J = 1,LA%NY
			KI = LA%POS(I,J,1)
			KJ = LA%POS(I,J,2)
			IF (KI.GE.1 .AND. KI.LT.LO%NX) THEN
			   IF (KJ.GE.1 .AND. KJ.LT.LO%NY) THEN
                  Z1 = LO%DEFORM(KI,KJ)*LA%CXY(I,J,1)
			      Z2 = LO%DEFORM(KI+1,KJ)*LA%CXY(I,J,2)
			      Z3 = LO%DEFORM(KI,KJ+1)*LA%CXY(I,J,3)
			      Z4 = LO%DEFORM(KI+1,KJ+1)*LA%CXY(I,J,4)
			      DEFORM(I,J) = Z1+Z2+Z3+Z4
			   ENDIF
			ENDIF
         ENDDO
      ENDDO 
!	  LA%Z(:,:,1) = LA%Z(:,:,1) + CELL(:,:)
	  LA%DEFORM(:,:) = DEFORM(:,:)
	  
	  RETURN  
	  END



!----------------------------------------------------------------------
      SUBROUTINE EDGE_INTERP_SC (LO,LA)
!......................................................................
!DESCRIPTION:
!	  #. VOLUME FLUXES ALONG FOUR BOUNDARIES OF CHILD GRID LAYER LA ARE 
!		 OBTAINED VIA INTERPOLATION FROM PARENT GRID LAYER LO AT 
!		 THE BEGINNING OF A TIME STEP;
!     #. DESIGNED FOR INTERPOLATION FROM SPHERICAL GRID LAYER, LO, TO 
!	     CARTESIAN GRID LAYERS, LA;
!	  #. REVERSE UNIVERSAL TRANSVERSE MERCATOR PROJECTION IS ADOPTED 
!		 TO PROJECT GRIDS IN CARTESIAN COORDINATE SYSTEM ONTO THE 
!		 EARTH SPHERE SURFACE (WGS84 ELLIPSOID);
!	  #. LA%UPZ = 
!				.TRUE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT 
!						  THE SAME COORDINATE SYSTEM;
!				.FALSE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT  
!						  DIFFERENT COORDINATE SYSTEM;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTE: 
!	  CREATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  UPDATED ON JAN06 2009 (XIAOMING WANG, GNS)
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO, LA
	  real Z1,Z2,Z3,Z4
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!.....GET FLUX THROUGH LEFT BOUNDARY OF LA
	  DO J = 1,LA%NY
		 KI = LA%POS(1,J,1)
		 KJ = LA%POS(1,J,2)
		 IF (LA%H(I,J)+LA%Z(I,J,2) .GT. GX) THEN
            Z1 = LO%M(KI,KJ,1)*LA%CXY(1,J,1)
		    Z2 = LO%M(KI+1,KJ,1)*LA%CXY(1,J,2)
		    Z3 = LO%M(KI,KJ+1,1)*LA%CXY(1,J,3)
		    Z4 = LO%M(KI+1,KJ+1,1)*LA%CXY(1,J,4)
		    LA%M(1,J,1) = Z1+Z2+Z3+Z4
		 ELSE
			LA%M(1,J,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH RIGHT BOUNDARY OF LA
	  DO J = 1,LA%NY
		 KI = LA%POS(LA%NX,J,1)
		 KJ = LA%POS(LA%NX,J,2)
		 IF (LA%H(I,J)+LA%Z(I,J,2) .GT. GX) THEN
            Z1 = LO%M(KI,KJ,1)*LA%CXY(LA%NX,J,1)
			Z2 = LO%M(KI+1,KJ,1)*LA%CXY(LA%NX,J,2)
			Z3 = LO%M(KI,KJ+1,1)*LA%CXY(LA%NX,J,3)
			Z4 = LO%M(KI+1,KJ+1,1)*LA%CXY(LA%NX,J,4)
			LA%M(LA%NX,J,1) = Z1+Z2+Z3+Z4
		 ELSE
			LA%M(LA%NX,J,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH BOTTOM BOUNDARY OF LA
	  DO I = 1,LA%NX
		 KI = LA%POS(I,1,1)
		 KJ = LA%POS(I,1,2)
		 IF (LA%H(I,J)+LA%Z(I,J,2) .GT. GX) THEN
            Z1 = LO%N(KI,KJ,1)*LA%CXY(I,1,1)
			Z2 = LO%N(KI+1,KJ,1)*LA%CXY(I,1,2)
			Z3 = LO%N(KI,KJ+1,1)*LA%CXY(I,1,3)
			Z4 = LO%N(KI+1,KJ+1,1)*LA%CXY(I,1,4)
			LA%N(I,1,1) = Z1+Z2+Z3+Z4
		 ELSE
			LA%N(I,1,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH TOP BOUNDARY OF LA
	  DO I = 1,LA%NX
		 KI = LA%POS(I,LA%NY,1)
		 KJ = LA%POS(I,LA%NY,2)
		 IF (LA%H(I,J)+LA%Z(I,J,2) .GT. GX) THEN
            Z1 = LO%N(KI,KJ,1)*LA%CXY(I,LA%NY,1)
			Z2 = LO%N(KI+1,KJ,1)*LA%CXY(I,LA%NY,2)
			Z3 = LO%N(KI,KJ+1,1)*LA%CXY(I,LA%NY,3)
			Z4 = LO%N(KI+1,KJ+1,1)*LA%CXY(I,LA%NY,4)
			LA%N(I,LA%NY,1) = Z1+Z2+Z3+Z4
		 ELSE
			LA%N(I,LA%NY,1) = 0.0
		 ENDIF
      ENDDO 

	  
	  RETURN  
	  END

!----------------------------------------------------------------------
      SUBROUTINE NEWQ_SC (LO,LA,T)
!......................................................................
!DESCRIPTION:
!	  #. VOLUME FLUXES ALONG FOUR BOUNDARIES OF CHILD GRID LAYER LA ARE 
!		 OBTAINED VIA INTERPOLATION FROM PARENT GRID LAYER LO AT TIME
!		 N*LO%DT AND (N+1)*LO%DT; VOLUME FLUXES OF LO AT THE BEGINNING
!		 OF A TIME STEP, N, IS KNOWN, BUT THEY ARE NOT KNOWN AT N+1.
!		 VOLUME FLUXES OF LO AT T = (N+1)*LO%DT ARE PREDICTED LOCALLY 
!		 ALONG THE FOUR CONNECTING BOUNDARIES OF BETWEEN LO AND LA. 
!		 THEN, INTERPOLATION IS CARRIED OUT TO GET THE FLUXES BETWEEN
!		 N*LO%DT AND (N+1)*LO%DT. THESE INTERPOLATED FLUXES ARE 
!		 INTERPOLATED AGAIN TO FORM THE FLUX BOUDNARY CONDITIONS OF LA.
!     #. DESIGNED FOR INTERPOLATION FROM SPHERICAL GRID LAYER, LO, TO 
!	     CARTESIAN GRID LAYERS, LA;
!	  #. REVERSE UNIVERSAL TRANSVERSE MERCATOR PROJECTION IS ADOPTED 
!		 TO PROJECT GRIDS IN CARTESIAN COORDINATE SYSTEM ONTO THE 
!		 EARTH SPHERE SURFACE (WGS84 ELLIPSOID);
!	  #. LA%UPZ = 
!				.TRUE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT 
!						  THE SAME COORDINATE SYSTEM;
!				.FALSE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT  
!						  DIFFERENT COORDINATE SYSTEM;
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTE: 
!	  #. CREATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN06 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR07 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE THE EFFICIENCY;
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO, LA
      INTEGER       :: KI,KJ,T
	  REAL XM1,XM2,YN1,YN2
	  REAL C1, C2
      DATA TWLVTH/0.08333333333333/
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

!.....FORECAST THE X FLUX THROUGH LEFT BOUNDARY AT T = N+1
	  IS = MIN0(LA%POS(1,1,1),LA%POS(1,LA%NY,1))-2
	  IE = MAX0(LA%POS(1,1,1),LA%POS(1,LA%NY,1))+2
	  JS = LA%CORNERS(3)-2
	  JE = LA%CORNERS(4)+2
	  SCM = 0.0
	  XM = 0.0
	  XN = 0.0
	  DO I=IS,IE
         IP1 = I+1
         DO J=JS,JE
			IF ((LO%H(I,J).GT.GX) .AND. (LO%H(IP1,J).GT.GX)) THEN
				JM1 = J-1
				JP1 = J+1
				IF (JM1.LE.1) JM1 = 1
				IF (JP1.GE.LO%NY) JP1 = LO%NY
				TOT_N = LO%N(I,J,1) + LO%N(IP1,J,1) + LO%N(I,JM1,1)	&
												+ LO%N(IP1,JM1,1)
				XM = LO%M(I,J,1) - LO%R2(I,J)*(LO%Z(IP1,J,2)		&
								- LO%Z(I,J,2)) + LO%R3(I,J)*TOT_N
				IF (LO%MODSCM .EQ. 0) THEN
					SCM =  LO%R2(I,J)*TWLVTH*((LO%Z(IP1,JP1,2)		&
								- 2*LO%Z(IP1,J,2)+LO%Z(IP1,JM1,2))	&
								- (LO%Z(I,JP1,2)-2*LO%Z(I,J,2)		&
								+ LO%Z(I,JM1,2)))
					XM = XM - SCM
				ENDIF
				IF (ABS(XM) .LT. EPS) XM = ZERO
				LO%M(I,J,2) = XM
			ELSE
			    LO%M(I,J,2) = 0.0
			END IF
         END DO
      END DO
!.....FORECAST THE X FLUX THROUGH RIGHT BOUNDARY AT T = N+1
	  IS = MIN0(LA%POS(LA%NX,1,1),LA%POS(LA%NX,LA%NY,1))-2
	  IE = MAX0(LA%POS(LA%NX,1,1),LA%POS(LA%NX,LA%NY,1))+2
	  JS = LA%CORNERS(3)-2
	  JE = LA%CORNERS(4)+2
	  SCM = 0.0
	  XM = 0.0
	  XN = 0.0
	  DO I=IS,IE
         IP1 = I+1
         DO J=JS,JE
			IF ((LO%H(I,J).GT.GX) .AND. (LO%H(IP1,J).GT.GX)) THEN
				JM1 = J-1
				JP1 = J+1
				IF (JM1.LE.1) JM1 = 1
				IF (JP1.GE.LO%NY) JP1 = LO%NY
				TOT_N = LO%N(I,J,1) + LO%N(IP1,J,1) + LO%N(I,JM1,1)	&
												+ LO%N(IP1,JM1,1)
				XM = LO%M(I,J,1) - LO%R2(I,J)*(LO%Z(IP1,J,2)		&
								- LO%Z(I,J,2)) + LO%R3(I,J)*TOT_N
				IF (LO%MODSCM .EQ. 0) THEN
					SCM =  LO%R2(I,J)*TWLVTH*((LO%Z(IP1,JP1,2)		&
								- 2*LO%Z(IP1,J,2)+LO%Z(IP1,JM1,2))	&
								- (LO%Z(I,JP1,2)-2*LO%Z(I,J,2)		&
								+ LO%Z(I,JM1,2)))
					XM = XM - SCM
				ENDIF
				IF (ABS(XM) .LT. EPS) XM = ZERO
				LO%M(I,J,2) = XM
			ELSE
			    LO%M(I,J,2) = 0.0
			END IF
         END DO
      END DO
!.....FORECAST THE Y FLUX THROUGH TOP BOUNDARY AT T = N+1
	  IS = LA%CORNERS(1)-2
	  IE = LA%CORNERS(2)+2
	  JS = MIN0(LA%POS(1,LA%NY,2),LA%POS(LA%NX,LA%NY,2))-2
	  JE = MAX0(LA%POS(1,LA%NY,2),LA%POS(LA%NX,LA%NY,2))+2
	  SCM = 0.0
	  XM = 0.0
	  XN = 0.0
      DO J=JS,JE
         JP1 = J+1
         DO I=IS,IE
			IF ((LO%H(I,J).GT.GX) .AND. (LO%H(I,JP1).GT.GX)) THEN
				IM1 = I-1
				IP1 = I+1
				IF (IM1.LE.1) IM1 = 1
				IF (IP1.GE.LO%NX) IP1 = LO%NX
				TOT_M = LO%M(IM1,J,1)+LO%M(IM1,JP1,1)+LO%M(I,J,1)	&
													+ LO%M(I,JP1,1)
				XN = LO%N(I,J,1) - LO%R4(I,J)*(LO%Z(I,JP1,2)		&
								- LO%Z(I,J,2))-LO%R5(I,J)*TOT_M
				IF (LO%MODSCM .EQ. 0) THEN
					SCM = LO%R4(I,J)*TWLVTH*((LO%Z(IP1,JP1,2)		&
								- 2*LO%Z(I,JP1,2)+LO%Z(IM1,JP1,2))	&
			                    - (LO%Z(IP1,J,2)-2*LO%Z(I,J,2)		&
								+ LO%Z(IM1,J,2)))
					XN = XN - SCM
				ENDIF
				IF (ABS(XN) .LT. EPS) XN = ZERO
				LO%N(I,J,2) = XN
			ELSE
			    LO%N(I,J,2) = 0.0
			END IF
         END DO
      END DO
!.....FORECAST THE Y FLUX THROUGH BOTTOM BOUNDARY AT T = N+1
	  IS = LA%CORNERS(1)-2
	  IE = LA%CORNERS(2)+2
	  JS = MIN0(LA%POS(1,1,2),LA%POS(LA%NX,1,2))-2
	  JE = MAX0(LA%POS(1,1,2),LA%POS(LA%NX,1,2))+2
	  SCM = 0.0
	  XM = 0.0
	  XN = 0.0
      DO J=JS,JE
         JP1 = J+1
         DO I=IS,IE
			IF ((LO%H(I,J).GT.GX) .AND. (LO%H(I,JP1).GT.GX)) THEN
				IM1 = I-1
				IP1 = I+1
				IF (IM1.LE.1) IM1 = 1
				IF (IP1.GE.LO%NX) IP1 = LO%NX
				TOT_M = LO%M(IM1,J,1)+LO%M(IM1,JP1,1)+LO%M(I,J,1)	&
													+ LO%M(I,JP1,1)
				XN = LO%N(I,J,1) - LO%R4(I,J)*(LO%Z(I,JP1,2)		&
								- LO%Z(I,J,2))-LO%R5(I,J)*TOT_M
				IF (LO%MODSCM .EQ. 0) THEN
					SCM = LO%R4(I,J)*TWLVTH*((LO%Z(IP1,JP1,2)		&
								- 2*LO%Z(I,JP1,2)+LO%Z(IM1,JP1,2))	&
			                    - (LO%Z(IP1,J,2)-2*LO%Z(I,J,2)		&
								+ LO%Z(IM1,J,2)))
					XN = XN - SCM
				ENDIF
				IF (ABS(XN) .LT. EPS) XN = ZERO
				LO%N(I,J,2) = XN
			ELSE
			    LO%N(I,J,2) = 0.0
			END IF
         END DO
      END DO

!.....INTERPOLATE FLUXES FROM LO TO LA THROUGH CONNECTING BOUNDARIES
	  C1 = (T-1.0)/LA%REL_TIME
	  C2 = 1.0 - C1
!.....GET FLUX THROUGH LEFT BOUNDARY OF LA
	  DO J = 1,LA%NY
		 KI = LA%POS(1,J,1)
		 KJ = LA%POS(1,J,2)
		 XM1 = 0.0
		 XM2 = 0.0
		 IF (LA%H(1,J)+LA%Z(1,J,2) .GT. GX) THEN
            Z1 = LO%M(KI,KJ,1)*LA%CXY(1,J,1)
			Z2 = LO%M(KI+1,KJ,1)*LA%CXY(1,J,2)
			Z3 = LO%M(KI,KJ+1,1)*LA%CXY(1,J,3)
			Z4 = LO%M(KI+1,KJ+1,1)*LA%CXY(1,J,4)
			XM1 = Z1+Z2+Z3+Z4
            Z1 = LO%M(KI,KJ,2)*LA%CXY(1,J,1)
			Z2 = LO%M(KI+1,KJ,2)*LA%CXY(1,J,2)
			Z3 = LO%M(KI,KJ+1,2)*LA%CXY(1,J,3)
			Z4 = LO%M(KI+1,KJ+1,2)*LA%CXY(1,J,4)
			XM2 = Z1+Z2+Z3+Z4
			LA%M(1,J,1) = C1*XM2 + C2*XM1
		 ELSE
			LA%M(1,J,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH RIGHT BOUNDARY OF LA
	  DO J = 1,LA%NY
		 KI = LA%POS(LA%NX,J,1)
		 KJ = LA%POS(LA%NX,J,2)
		 XM1 = 0.0
		 XM2 = 0.0
		 IF (LA%H(LA%NX,J)+LA%Z(LA%NX,J,2) .GT. GX) THEN
            Z1 = LO%M(KI,KJ,1)*LA%CXY(LA%NX,J,1)
			Z2 = LO%M(KI+1,KJ,1)*LA%CXY(LA%NX,J,2)
			Z3 = LO%M(KI,KJ+1,1)*LA%CXY(LA%NX,J,3)
			Z4 = LO%M(KI+1,KJ+1,1)*LA%CXY(LA%NX,J,4)
			XM1 = Z1+Z2+Z3+Z4
            Z1 = LO%M(KI,KJ,2)*LA%CXY(LA%NX,J,1)
			Z2 = LO%M(KI+1,KJ,2)*LA%CXY(LA%NX,J,2)
			Z3 = LO%M(KI,KJ+1,2)*LA%CXY(LA%NX,J,3)
			Z4 = LO%M(KI+1,KJ+1,2)*LA%CXY(LA%NX,J,4)
			XM2 = Z1+Z2+Z3+Z4
			LA%M(LA%NX,J,1) = C1*XM2 + C2*XM1
		 ELSE
			LA%M(LA%NX,J,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH BOTTOM BOUNDARY OF LA
	  DO I = 1,LA%NX
		 KI = LA%POS(I,1,1)
		 KJ = LA%POS(I,1,2)
		 YN1 = 0.0
		 YN2 = 0.0
		 IF (LA%H(I,1)+LA%Z(I,1,2) .GT. GX) THEN
            Z1 = LO%N(KI,KJ,1)*LA%CXY(I,1,1)
			Z2 = LO%N(KI+1,KJ,1)*LA%CXY(I,1,2)
			Z3 = LO%N(KI,KJ+1,1)*LA%CXY(I,1,3)
			Z4 = LO%N(KI+1,KJ+1,1)*LA%CXY(I,1,4)
			YN1 = Z1+Z2+Z3+Z4
            Z1 = LO%N(KI,KJ,2)*LA%CXY(I,1,1)
			Z2 = LO%N(KI+1,KJ,2)*LA%CXY(I,1,2)
			Z3 = LO%N(KI,KJ+1,2)*LA%CXY(I,1,3)
			Z4 = LO%N(KI+1,KJ+1,2)*LA%CXY(I,1,4)
			YN2 = Z1+Z2+Z3+Z4
			LA%N(I,1,1) = C1*YN2 + C2*YN1
		 ELSE
			LA%N(I,1,1) = 0.0
		 ENDIF
      ENDDO 

!.....GET FLUX THROUGH TOP BOUNDARY OF LA
	  DO I = 1,LA%NX
		 KI = LA%POS(I,LA%NY,1)
		 KJ = LA%POS(I,LA%NY,2)
		 YN1 = 0.0
		 YN2 = 0.0
		 IF (LA%H(I,LA%NY)+LA%Z(I,LA%NY,2) .GT. GX) THEN
            Z1 = LO%N(KI,KJ,1)*LA%CXY(I,LA%NY,1)
			Z2 = LO%N(KI+1,KJ,1)*LA%CXY(I,LA%NY,2)
			Z3 = LO%N(KI,KJ+1,1)*LA%CXY(I,LA%NY,3)
			Z4 = LO%N(KI+1,KJ+1,1)*LA%CXY(I,LA%NY,4)
			YN1 = Z1+Z2+Z3+Z4
            Z1 = LO%N(KI,KJ,2)*LA%CXY(I,LA%NY,1)
			Z2 = LO%N(KI+1,KJ,2)*LA%CXY(I,LA%NY,2)
			Z3 = LO%N(KI,KJ+1,2)*LA%CXY(I,LA%NY,3)
			Z4 = LO%N(KI+1,KJ+1,2)*LA%CXY(I,LA%NY,4)
			YN2 = Z1+Z2+Z3+Z4
!			LA%N(I,LA%NY,1) = ((T-1.0)/LA%REL_TIME)*(YN2-YN1)+YN1
			LA%N(I,LA%NY,1) = C1*YN2 + C1*YN1
		 ELSE
			LA%N(I,LA%NY,1) = 0.0
		 ENDIF
      ENDDO 
	  
	  RETURN  
	  END

!----------------------------------------------------------------------
      SUBROUTINE JNZ_SC (LO,LA)
!......................................................................
!DESCRIPTION:
!	  #. UPDATE FREE SURFACE DISPLACEMENT FROM INNER REGION (LA)  
!	     INTO OUTER REGION (LO)
!     #. DESIGNED FOR INTERPOLATION FROM SPHERICAL GRID LAYER, LO, TO 
!	     CARTESIAN GRID LAYERS, LA;
!	  #. REVERSE UNIVERSAL TRANSVERSE MERCATOR PROJECTION IS ADOPTED 
!		 TO PROJECT GRIDS IN CARTESIAN COORDINATE SYSTEM ONTO THE 
!		 EARTH SPHERE SURFACE (WGS84 ELLIPSOID);
!	  #. LA%UPZ = 
!				.TRUE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT 
!						  THE SAME COORDINATE SYSTEM;
!				.FALSE. - PARENT GRID, LO, AND CHILD GRID, LA, ADOPT  
!						  DIFFERENT COORDINATE SYSTEM;
!     #. RUNS MUCH SLOWER THAN JNZ, NOT EFFICIENT!!! (IMPROVED)
!	  #. SC_OPTION: COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!			 = 0: TRADITIONAL COUPLING SCHEME BETWEEEN SPH AND CART;
!			 = 1: IMPROVED COUPLING SCHEME BETWEEN SPH AND CART;
!NOTE: 
!	  #. CREATED ON JAN05 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON JAN06 2009 (XIAOMING WANG, GNS)
!	  #. UPDATED ON APR03 2009 (XIAOMING WANG, GNS)
!		 1. IMPROVE COUPLING SCHEME BETWEEN SPHERICAL AND CARTESIAN
!----------------------------------------------------------------------
      USE LAYER_PARAMS
      TYPE (LAYER)	:: LO, LA
      INTEGER COUNT,KI,KJ,I,J
	  REAL SUM, HALF
	  REAL Z(LO%NX,LO%NY),S(LO%NX,LO%NY)
	  COMMON /CONS/ ELMAX,GRAV,PI,R_EARTH,GX,EPS,ZERO,ONE,NUM_GRID,	&
					NUM_FLT,V_LIMIT,RAD_DEG,RAD_MIN

	  Z = 0.0
	  S = 0.0

      HALF = REAL(LA%REL_SIZE*LA%REL_SIZE)/2.0

	  DO I = 1,LA%NX
	     DO J = 1,LA%NY
			IF (LA%H(I,J)+LA%Z(I,J,2) .GT. GX) THEN
			   KI = LA%POS(I,J,1)
			   KJ = LA%POS(I,J,2)
			   Z(KI,KJ) = Z(KI,KJ) + LA%Z(I,J,2)
			   S(KI,KJ) = S(KI,KJ) + 1.0
			ENDIF
		 ENDDO
	  ENDDO

	  DO I = LA%CORNERS(1),LA%CORNERS(2)
		 DO J = LA%CORNERS(3),LA%CORNERS(4)
			IF (S(I,J) .GT. HALF) THEN
			   LO%Z(I,J,2) = Z(I,J)/S(I,J)
			ENDIF
		 ENDDO
	  ENDDO
	  
	  RETURN  
	  END

