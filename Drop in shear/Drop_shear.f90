
!****************************************************************
!
!	Improved conservative phase-field LBM solver
!	for two-phase flows in a periodic domain (2D)
!
!----------------------------------------------------------------
!	Based on the following paper:
!
!	A. Fakhari, T. Mitchell, C. Leonardi, and D. Bolster,
!	"Improved locality of the phase-field lattice-Boltzmann model
!	for immiscible fluids at high density ratios",
!	Physical Review E 96, 053301 (2017)
!----------------------------------------------------------------
!
!	written by Abbas Fakhari 10/30/2016
!	03/23/2017:	updated					 
!	10/18/2018:	minor updates
!
!****************************************************************

MODULE SHARE
    IMPLICIT NONE
    INTEGER :: t, X, Y

    INTEGER,PARAMETER :: L0 = 256  !512

    INTEGER,PARAMETER :: tf   = 50000
    INTEGER,PARAMETER :: step = tf/50

    INTEGER,PARAMETER :: Nx = 160
    INTEGER,PARAMETER :: Ny = 160

    INTEGER,PARAMETER :: X0 = Nx/2 + 1
    INTEGER,PARAMETER :: Y0 = Ny/2 + 1
 

    INTEGER,PARAMETER :: ex(0:8) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
    INTEGER,PARAMETER :: ey(0:8) = [0, 0, 1, 0,-1, 1, 1,-1,-1]
    REAL(8),PARAMETER :: Wa(0:8) = [16,4, 4, 4, 4, 1, 1, 1, 1] / 36.d0

    ! REAL(8),PARAMETER :: R = L0 ! Hit&trail
    REAL(8),PARAMETER :: R = 20 !L0/8.d0 ! original
!    REAL(8),PARAMETER :: a = R/64.d0 ! new added
!    REAL(8),PARAMETER :: b = R/32.d0 ! new added
    
    REAL(8),PARAMETER :: Rhol  = 0.1832d0
    REAL(8),PARAMETER :: Rhoh  = 0.1832d0
    REAL(8),PARAMETER :: dRho3 = (Rhoh - Rhol)/3
    REAL(8),PARAMETER :: dRho2 = (Rhoh - Rhol)/2 ! added for gravity or density force.

    REAL(8),PARAMETER :: Sigma = 0.01d0
    REAL(8),PARAMETER :: gY = 3.05e-5
    REAL(8),PARAMETER :: pai = 3.141592653589793d0 ! added value of pai

    REAL(8),PARAMETER :: tau_H = 1.d0 !0.8d0 !0.3d0 + 0.5d0 !(3.d0*U0*real(Nx))/Re mu was taken 0.005
    REAL(8),PARAMETER :: tau_L = 1.d0  !0.8d0  !0.3d0 + 0.5d0 !(3.d0*U0*real(Nx))/Re !10.d0*(3.d0*U0*real(Nx))/Re
    REAL(8),PARAMETER :: C_H = 1.d0 !0.3d0 + 0.5d0 !(3.d0*U0*real(Nx))/Re
    REAL(8),PARAMETER :: C_L = 0.d0  !0.3d0 + 0.5d0 !(3.d0*U0*real(Nx))/Re !10.d0*(3.d0*U0*real(Nx))/Re

!	REAL(8),PARAMETER :: tau = 0.3d0 + 0.5d0
!	REAL(8),PARAMETER :: s8 = 1.d0/tau

    REAL(8),PARAMETER :: W    = 3
    REAL(8),PARAMETER :: Beta = 12.d0 * Sigma/W
    REAL(8),PARAMETER :: k    = 1.5d0 * Sigma*W

    REAL(8),PARAMETER :: M   = 0.15d0
    REAL(8),PARAMETER :: w_c = 1.d0/(0.5d0 + 3.d0*M)

    REAL(8) :: h(0:8,0:Nx+1,0:Ny+1), g(0:8,0:Nx+1,0:Ny+1)
    REAL(8) :: Gamma(0:8), Ga_Wa(0:8), heq(0:8), geq(0:8), hlp(0:8), eF(0:8)
    REAL(8) :: C(0:Nx+1,0:Ny+1), P(Nx,Ny), mu(Nx,Ny), DcDx(Nx,Ny), DcDy(Nx,Ny)
    REAL(8) :: Rho(Nx,Ny), Ux(Nx, Ny), Uy(Nx,Ny), ni(Nx,Ny), nj(Nx,Ny)!, tau(Nx,Ny), s8(Nx,Ny)
    REAL(8) :: tmp, Ri, Fx, Fy

END MODULE SHARE

!**********************************************************************

PROGRAM Conservative_PhaseField_Periodic_2D
    USE SHARE
    IMPLICIT NONE

    CALL Initialize_distributions

!	OPEN (1, file = 'LBM.plt')
!	WRITE(1,*) 'Variables = X, Y, Ux, Uy, C, P, Rho'

    !=========================================================================================
!	PRINT '(/A/)', '   tf    Sigma     W      M      R      tau    s8    Rhol   Rhoh     L0'
!	PRINT '(I6,F8.4,F7.1,2F8.3,4F7.3,I7)', tf, Sigma, W, M, R, tau, s8, Rhol, Rhoh, L0

    PRINT '(/A6,5A12,A12/)', 't', 'C_min', 'C_max', 'Ux_max', 'Uy_max', '|U_max|', 'Mass_C'
    !=========================================================================================

    DO t = 0, tf

           IF( MOD(t,500)==0 )THEN
             CALL Track_Drop_Axes
           END IF

        IF( ISNAN(C(2,2)) )THEN
            PRINT '(/A,I5/)', '!!! THE PROGRAM DIVERGED AT t =', t
            STOP
        ELSEIF( MOD(t,step)==0 )THEN
            CALL RESULTS_Output
        !	CALL Interface_tracking
            PRINT '(I7,2F12.6,3E12.3,F11.2)', t, MINVAL(C), MAXVAL(C), MAXVAL(ABS(Ux)), MAXVAL(ABS(Uy)), DSQRT( MAXVAL(Ux**2+Uy**2) ), SUM(C(1:Nx,1:Ny))
        END IF

        CALL Improved_PhaseField_h_g

    END DO

    CALL CPU_TIME( tmp )
    PRINT '(/A)',' *****************************************'
    PRINT '(A,F12.1)', ' Time Elapsed:', tmp
    PRINT '( A)',' *****************************************'

END

!**********************************************************************

SUBROUTINE Initialize_distributions
    USE SHARE
    IMPLICIT NONE

    P  = 0
    Ux = 0
    Uy = 0

    DO Y = 0, Ny+1  !!1, Ny
    DO X = 0, Nx+1  !!1, Nx

     !	 Ri = 0.01d0*L0*dcos( real(X0) / L0 * 2.d0 * pai) ! was trying for RTI
        Ri = DSQRT( (X-(X0-0.5d0))**2.d0 + (Y-(Y0-0.5d0))**2.d0)    ! equation for circle.
     !    Ri = DSQRT( ((X-(X0))*(X-(X0)))/(a*a) + ((Y-(Y0))*(Y-(Y0)))/(b*b) )  ! equation for ellipse.
        C(X,Y) = 0.5d0 + 0.5d0 * TANH(2*(R-Ri)/W)   !drop
    !    C(X,Y) = 0.5d0 - 0.5d0 * TANH(2*(R-Ri)/W)    !bubble

    END DO
    END DO

    CALL Boundary_Conditions_C( C )

    CALL Chemical_Potential

    CALL Isotropic_Gradient( C, DcDx, DcDy )

    CALL normal_FD

    DO Y = 1, Ny
    DO X = 1, Nx

        Rho(X,Y) = Rhol + C(X,Y) * (Rhoh - Rhol)

        P(X,Y) = P(X,Y) + C(X,Y) * Sigma/R /(Rho(X,Y)/3)    !in 2D (drop)
    !    P(X,Y) = P(X,Y) - C(X,Y) * Sigma/R /(Rho(X,Y)/3)    !in 2D (bubble)

        CALL Equilibrium_new( Ux(X,Y), Uy(X,Y) )

        Gamma(:) = Ga_Wa(:) + Wa(:)

        !*******************	heq

        eF(:)  = ( 1.d0 - 4.d0*(C(X,Y)-0.5d0)**2.d0 )/W * ( ex(:)*ni(X,Y) + ey(:)*nj(X,Y) )

        hlp(:) = Wa(:) * eF(:)

        h(:, X,Y) = C(X,Y)*Gamma(:) - 0.5d0 * hlp(:)

        !*******************	geq

        g(:, X,Y) = P(X,Y) * Wa(:) + Ga_Wa(:)

    END DO
    END DO

END

!**********************************************************************

SUBROUTINE Equilibrium_new( U, V )
    USE SHARE, ONLY: ex, ey, Wa, Ga_Wa
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: U, V

    REAL(8) :: U2, eU(0:8)

    U2 = U*U + V*V

    eU(:) = ex(:) * U  + ey(:) * V

    Ga_Wa(:) = Wa(:) * ( eU(:)*(3.d0 + 4.5d0*eU(:)) - 1.5d0*U2 )

END

!**********************************************************************

SUBROUTINE Improved_PhaseField_h_g
    USE SHARE, ONLY: h, g
    IMPLICIT NONE

    CALL Collision_h_g

    CALL Streaming( h )
    CALL Streaming( g )

    CALL Boundary_Conditions_f( h )
    CALL Boundary_Conditions_f( g )

    CALL Macroscopic_Properties_h
    CALL Macroscopic_Properties_g

END

!**********************************************************************

SUBROUTINE Collision_h_g
    USE SHARE
    USE, INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE

    REAL(8) :: FpX, FpY, FmX, FmY, FbX, FbY
    REAL(8) :: gneq(9)
    REAL(8) :: tau, s8
!	REAL(8) :: C_local


!    CALL Boundary_Conditions_C( C )             !no need (already called)
!!	CALL Isotropic_Gradient( C, DcDx, DcDy )	!no need (already called)

    CALL normal_FD

    DO Y = 1, Ny
    DO X = 1, Nx

        tau = tau_L + (C(X,Y) - C_L) * (tau_H - tau_L) / (C_H - C_L)
    !    tau = tau_L + C(X,Y)  * (tau_H - tau_L) 
    !	tau = 1.0d0 / ( (1.0d0 - C(X,Y)-C_L)/tau_L + C(X,Y)/tau_H )
    !    s8(X,Y)  = 1.0d0 / tau(X,Y)

    !  IF (X == Nx/2 .AND. Y == Ny/2) PRINT *, "C, tau:", C(X,Y), tau

    !    C_local = MAX(0.0d0, MIN(1.0d0, C(X,Y)))
    !    tau = 1.0d0 / ( (1.0d0 - C_local)/tau_L + C_local/tau_H )
!
        s8  = 1.0d0/tau

        CALL Equilibrium_new( Ux(X,Y), Uy(X,Y))

        Gamma(:) = Ga_Wa(:) + Wa(:)

        !*******************	COLLISION (h)

        eF(:)  = ( 1.d0 - 4.d0*(C(X,Y)-0.5d0)**2.d0 )/W * ( ex(:)*ni(X,Y) + ey(:)*nj(X,Y) )

        hlp(:) = Wa(:) * eF(:)

        heq(:) = C(X,Y)*Gamma(:) - 0.5d0 * hlp(:)

        h(:, X,Y) = h(:, X,Y) * (1.d0-w_c) + heq(:) * w_c + hlp(:)

        !*******************	COLLISION (g)	******************
        !*******************	calculate the forcing terms

        FbX = 0.d0
    !	FbY = - ( Rho(X,Y) - dRho2 ) * gY
        FbY =  -(Rho(X,Y) * gY)

        FpX = - P(X,Y) * dRho3 * DcDx(X,Y)
        FpY = - P(X,Y) * dRho3 * DcDy(X,Y)

        geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:)

        IF (.NOT. IEEE_IS_FINITE(C(X,Y))) THEN
          PRINT *, "Bad C at X,Y=", X, Y, " C=", C(X,Y)
          STOP
        END IF

        gneq(:) = g(:,X,Y) - geq(:)
        CALL Compute_Viscous_Force_BGK(tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y), geq(:), FmX, FmY)
    !	CALL Calculate_Viscous_Force( tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y)-geq(:), FmX, FmY )
    !	CALL Calculate_Viscous_Force( tau_H, tau_L, DcDx(X,Y), DcDy(X,Y), gneq(:), C(X,Y), FmX, FmY )

        !Fx = mu(X,Y) * DcDx(X,Y) 
        Fx = mu(X,Y) * DcDx(X,Y) + FpX  + FmX ! 
        Fy = mu(X,Y) * DcDy(X,Y) + FpY  + FmY ! 
        !Fy = mu(X,Y) * DcDy(X,Y) 

        eF(:) = ex(:) * Fx + ey(:) * Fy

        hlp(:) = 3.d0 * Wa(:) * eF(:) / Rho(X,Y)

        geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:) - 0.5d0 * hlp(:)

        g(:, X,Y) = g(:, X,Y) * (1-s8) + geq(:) * s8 + hlp(:)
!		g(:, X,Y) = g(:, X,Y) * (1-s8(X,Y)) + geq(:) * s8(X,Y) + hlp(:)

    END DO
    END DO

END

!**********************************************************************

SUBROUTINE Boundary_Conditions_f( f )
    USE SHARE, ONLY: Nx, Ny, Rho, Wa, ex
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: f(0:8, 0:Nx+1,0:Ny+1)
    REAL(8), PARAMETER :: Uw = 0.5d0 !0.01652d0
    INTEGER, PARAMETER :: opp(0:8) = [0, 3, 4, 1, 2, 7, 8, 5, 6]
    INTEGER :: X


 !-- left and right boundaries
    f(:,  0  ,:) = f(:, Nx,:)   !periodic
    f(:, Nx+1,:) = f(:, 1 ,:)   !periodic

 !-- bottom and top boundaries
  !	f(:, :, 0  ) = f(:, :,Ny)	!periodic
!	f(:, :,Ny+1) = f(:, :,1 )	!periodic

! ******************* solid boundary condition bounce back *****************************
 DO X = 1, Nx
    ! Bottom wall (iY = 0)        
    f(2, X, 0) = f(4, X, 0)
    f(5, X, 0) = f(7, X, 0) + 6*Wa(5)*Rho(X,1)*(-Uw)*ex(5) 
    f(6, X, 0) = f(8, X, 0) + 6*Wa(6)*Rho(X,1)*(-Uw)*ex(6)

    ! Top wall (iY = Ny+1)
    f(4, X, Ny+1) = f(2, X, Ny+1)
    f(7, X, Ny+1) = f(5, X, Ny+1) + 6*Wa(7)*Rho(X,Ny)*Uw*ex(7)
    f(8, X, Ny+1) = f(6, X, Ny+1) + 6*Wa(8)*Rho(X,Ny)*Uw*ex(8)
  END DO
! ******************* solid boundary condition bounce back *****************************
END

!**********************************************************************

SUBROUTINE Boundary_Conditions_C( A )
    USE SHARE, ONLY: Nx, Ny
    IMPLICIT NONE
    REAL(8),INTENT(INOUT) :: A(0:Nx+1,0:Ny+1)

    CALL PeriodicXY_C( A )

END

!**********************************************************************

SUBROUTINE PeriodicXY_C( A )
    USE SHARE, ONLY: Nx, Ny
    IMPLICIT NONE
    REAL(8),INTENT(INOUT):: A(0:Nx+1,0:Ny+1)
    REAL(8), PARAMETER :: C0 = -0.5d0
    REAL(8), PARAMETER :: Uw = 0.005d0
    INTEGER :: X

 !-- bottom and top boundaries
  !  A(:,  0 ) = A(:,Ny  )   ! periodic
  !  A(:,Ny+1) = A(:, 1  )   ! periodic

 !-- left and right boundaries
     A(  0 ,:) = A(Nx  ,:)  ! periodic
     A(Nx+1,:) = A( 1  ,:)   ! periodic

! ******************* solid boundaru condition bounce back *****************************
    DO X = 1, Nx
!    ! Bottom wall (iY = 0)
        A(X, 0) = A(X, 1)
        A(X, 0) = A(X, 1)
        A(X, 0) = A(X, 1)

    !    A(X, 0)    = C0
    !    A(X, Ny+1) = C0


    !    Ux(X, 0)    = -Uw
    !    Ux(X, Ny+1) = +Uw
    !    Uy(X, 0)    = 0.d0
    !    Uy(X, Ny+1) = 0.d0

!    ! Top wall (iY = Ny+1)
        A(X, Ny+1) = A(X, Ny)
        A(X, Ny+1) = A(X, Ny)
        A(X, Ny+1) = A(X, Ny)
    END DO
! ******************* solid boundaru condition bounce back *****************************  

END

!**********************************************************************

SUBROUTINE Streaming( f )
    USE SHARE, ONLY: X, Y, Nx, Ny, ex, ey
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: f(0:8, 0:Nx+1,0:Ny+1)
    
    INTEGER :: I
    REAL(8) :: fnew(8,Nx,Ny)

    DO Y=1,Ny
    DO X=1,Nx
        DO I = 1, 8
            fnew(I,X,Y) = f(I,X-ex(I),Y-ey(I))
        END DO
    END DO
    END DO

    f(1:8,1:Nx,1:Ny) = fnew

END

!**********************************************************************

SUBROUTINE Macroscopic_Properties_h
    USE SHARE
    IMPLICIT NONE

    C(1:Nx,1:Ny) = SUM( h(:, 1:Nx,1:Ny), 1 )

    Rho = Rhol + C(1:Nx,1:Ny) * (Rhoh - Rhol)

END

!**********************************************************************

SUBROUTINE Macroscopic_Properties_g
    USE SHARE
    IMPLICIT NONE

    REAL(8) :: FpX, FpY, FmX, FmY, FbX, FbY
!    REAL(8) :: gneq(9)

    REAL(8) :: tau

!    CALL Macroscopic_Properties_h

    CALL Boundary_Conditions_C( C )

    CALL Chemical_Potential

    CALL Isotropic_Gradient( C, DcDx, DcDy )

    DO Y = 1, Ny
    DO X = 1, Nx

        tau = tau_L + (C(X,Y) - C_L) * (tau_H - tau_L) / (C_H - C_L)

    !    tau = tau_L + C(X,Y) * (tau_H - tau_L)
    !    s8  = 1.0d0 / tau

        P(X,Y) = SUM( g(:, X,Y) )

        FpX = - P(X,Y) * dRho3 * DcDx(X,Y)
        FpY = - P(X,Y) * dRho3 * DcDy(X,Y)

        FbX = 0.d0
    !	FbY = - ( Rho(X,Y) - dRho2 ) * gY
        FbY =  -(Rho(X,Y) * gY)

        CALL Equilibrium_new( Ux(X,Y), Uy(X,Y) )
        geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:)
    !	geq(:) = P(X,Y) * Wa(:) + Ga_Wa(:)
        
        CALL Compute_Viscous_Force_BGK(tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y), geq(:), FmX, FmY)

    !	CALL Calculate_Viscous_Force( tau, DcDx(X,Y), DcDy(X,Y), g(:,X,Y)-geq(:), FmX, FmY )
    !	CALL Calculate_Viscous_Force( tau_H, tau_L, DcDx(X,Y), DcDy(X,Y), gneq(:), C(X,Y), FmX, FmY )

        Fx = mu(X,Y) * DcDx(X,Y) + FpX + FmX
    !    Fx = mu(X,Y) * DcDx(X,Y) + FpX + FmX ! 
    !    Fy = mu(X,Y) * DcDy(X,Y) + FpY + FmY ! 
        Fy = mu(X,Y) * DcDy(X,Y) + FpY + FmY

        Ux(X,Y) = g(1,X,Y)-g(3,X,Y)+g(5,X,Y)-g(6,X,Y)-g(7,X,Y)+g(8,X,Y) + 0.5d0*Fx/Rho(X,Y)
        Uy(X,Y) = g(2,X,Y)-g(4,X,Y)+g(5,X,Y)+g(6,X,Y)-g(7,X,Y)-g(8,X,Y) + 0.5d0*Fy/Rho(X,Y)

    END DO
    END DO

END

!**********************************************************************

SUBROUTINE Chemical_Potential
    USE SHARE, ONLY: X, Y, Nx, Ny, Beta, k, C, mu
    IMPLICIT NONE

    REAL(8) :: D2C

    DO Y = 1, Ny
    DO X = 1, Nx

        D2C = ( C(X-1,Y-1)+C(X+1,Y-1)+C(X-1,Y+1)+C(X+1,Y+1) &
            +4*(C(X  ,Y-1)+C(X-1,Y  )+C(X+1,Y  )+C(X  ,Y+1)) - 20*C(X,Y) )/6

        mu(X,Y) = 4*Beta * C(X,Y) * (C(X,Y)-1.d0) * (C(X,Y)-0.5d0) - k*D2C

    END DO
    END DO

END

!**********************************************************************

SUBROUTINE normal_FD
    USE SHARE, ONLY: X, Y, Nx, Ny, DcDx, DcDy, tmp, ni, nj
    IMPLICIT NONE

    DO Y = 1, Ny
    DO X = 1, Nx

        tmp = DSQRT( DcDx(X,Y)**2 + DcDy(X,Y)**2 + 1.d-32 )

        ni(X,Y) = DcDx(X,Y) / tmp
        nj(X,Y) = DcDy(X,Y) / tmp

    END DO
    END DO

END

!**********************************************************************

SUBROUTINE Isotropic_Gradient( C, DcDx, DcDy )
    USE SHARE, ONLY: X, Y, Nx, Ny
    IMPLICIT NONE
    REAL(8),INTENT(IN) :: C(0:Nx+1,0:Ny+1)
    REAL(8),INTENT(OUT):: DcDx(Nx,Ny), DcDy(Nx,Ny)

    DO Y = 1, Ny
    DO X = 1, Nx

        DcDx(X,Y) = (C(X+1,Y  ) - C(X-1,Y  ))/3 + ( C(X+1,Y-1) + C(X+1,Y+1) - C(X-1,Y-1) - C(X-1,Y+1))/12
        DcDy(X,Y) = (C(X  ,Y+1) - C(X  ,Y-1))/3 + ( C(X-1,Y+1) + C(X+1,Y+1) - C(X-1,Y-1) - C(X+1,Y-1))/12

    END DO
    END DO

END

SUBROUTINE Compute_Viscous_Force_BGK(tau, DcDx, DcDy, g, geq, FmX, FmY)
  USE SHARE, ONLY: ex, ey, Rhoh, Rhol
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  REAL(8), INTENT(IN)  :: tau, DcDx, DcDy
  REAL(8), INTENT(IN)  :: g(0:8), geq(0:8)
  REAL(8), INTENT(OUT) :: FmX, FmY

  REAL(8) :: gneq(0:8)
  REAL(8) :: sxx, sxy, syy
!  INTEGER :: i

  gneq(:) = g(:) - geq(:)

  ! Compute stress tensor directly
 ! sxx = 0.0d0
 ! sxy = 0.0d0
 ! syy = 0.0d0

    sxx = SUM( gneq(1:) * ex(1:) * ex(1:) )
    sxy = SUM( gneq(1:) * ex(1:) * ey(1:) )
    syy = SUM( gneq(1:) * ey(1:) * ey(1:) )

  ! Final force calculation
  FmX = (0.5d0 - tau) / tau * (sxx * DcDx + sxy * DcDy) * (Rhoh - Rhol)
  FmY = (0.5d0 - tau) / tau * (sxy * DcDx + syy * DcDy) * (Rhoh - Rhol)

  ! Safety check
  IF (.NOT. IEEE_IS_FINITE(FmX) .OR. .NOT. IEEE_IS_FINITE(FmY)) THEN
    PRINT *, "NaN in FmX or FmY! tau=", tau
    STOP
  END IF
END SUBROUTINE


SUBROUTINE Track_Drop_Axes
  USE SHARE, ONLY: X, Y, Nx, Ny, T, L0, Rho, Rhoh, Rhol
  IMPLICIT NONE

  INTEGER :: Xc, Yc
  REAL(8) :: rho_d, rho_u, pos_d, pos_u
  REAL(8) :: x_left, x_right, y_bot, y_top
  REAL(8) :: a_axis, b_axis
  REAL(8) :: rho_mid

  ! Center of the drop (assumed initialized at domain center)
  Xc = Nx / 2
  Yc = Ny / 2
  rho_mid = 0.5d0 * (Rhoh + Rhol)

  ! === Vertical (Y-direction) interface detection ===
  DO Y = 1, Ny
    IF (Rho(Xc, Y) < rho_mid) THEN
      rho_d = Rho(Xc, Y)
      pos_d = Y
      EXIT
    END IF
  END DO
  DO Y = Ny, 1, -1
    IF (Rho(Xc, Y) > rho_mid) THEN
      rho_u = Rho(Xc, Y)
      pos_u = Y
      EXIT
    END IF
  END DO
  y_bot = pos_d !+ (pos_u - pos_d) / (rho_u - rho_d) * (rho_mid - rho_d)
  y_top = pos_u !+ (pos_d - pos_u) / (rho_d - rho_u) * (rho_mid - rho_u)

  ! === Horizontal (X-direction) interface detection ===
  DO X = 1, Nx
    IF (Rho(X, Yc) < rho_mid) THEN
      rho_d = Rho(X, Yc)
      pos_d = X
      EXIT
    END IF
  END DO
  DO X = Nx, 1, -1
    IF (Rho(X, Yc) > rho_mid) THEN
      rho_u = Rho(X, Yc)
      pos_u = X
      EXIT
    END IF
  END DO
  x_left  = pos_d !+ (pos_u - pos_d) / (rho_u - rho_d) * (rho_mid - rho_d)
  x_right = pos_u !+ (pos_d - pos_u) / (rho_d - rho_u) * (rho_mid - rho_u)

  ! === Compute full axis lengths (in lattice units) ===
!  a_axis = (x_right - x_left) / L0
!  b_axis = (y_top - y_bot) / L0

  a_axis = (x_right - x_left) 
  b_axis = (y_top - y_bot) 
  ! === Write to file ===
  OPEN(99, FILE="drop_axes.dat", ACCESS="APPEND")
  WRITE(99,'(I10,2F15.6)') T, a_axis, b_axis
  CLOSE(99)

END SUBROUTINE Track_Drop_Axes





!**********************************************************************
!SUBROUTINE Calculate_Viscous_Force( tau, DcDx, DcDy, gneq, FmX, FmY )
!!SUBROUTINE Calculate_Viscous_Force( tau_H, tau_L, DcDx, DcDy, gneq, C, FmX, FmY )
!	IMPLICIT NONE
!!	REAL(8),INTENT(IN)  :: tau_H, tau_L, DcDx, DcDy, gneq(9), C
!    REAL(8),INTENT(IN)  :: DcDx, DcDy, gneq(0:8)
!	REAL(8),INTENT(OUT) :: FmX, FmY
!	REAL(8) :: tau
!	
!
!	 CALL Calculate_Viscous_Force_BGK( tau, DcDx, DcDy, gneq, FmX, FmY )
!	 PRINT *, "NaN or Inf in tau: ", tau
!	!CALL Calculate_Viscous_Force_BGK( tau_H, tau_L, DcDx, DcDy, gneq, C, FmX, FmY )
!
!END

!**********************************************************************

!SUBROUTINE Calculate_Viscous_Force_BGK( tau_H, tau_L, DcDx, DcDy, gneq, FmX, FmY, C )
!	USE SHARE, ONLY: X, Y, Nx, Ny, Rhoh, Rhol
!!	USE SHARE
!	IMPLICIT NONE
!	REAL(8),INTENT(IN)  :: tau_H, tau_L, DcDx, DcDy, gneq(0:8), C(0:Nx+1,0:Ny+1)
!	REAL(8),INTENT(OUT) :: FmX, FmY
!
!	REAL(8) :: tau, s8, sxx, sxy, syy
!
!	DO Y = 1, Ny
!	DO X = 1, Nx
!
!	  tau = tau_L + C(X,Y) * (tau_H - tau_L)
!      s8  = 1.0d0 / tau
!	
!	CALL Calculate_Stress_Tensor_BGK( gneq(1:), sxx, sxy, syy )
!
!	FmX = (0.5d0-tau)/tau * (sxx*DcDx+sxy*DcDy) * (Rhoh-Rhol)
!	FmY = (0.5d0-tau)/tau * (sxy*DcDx+syy*DcDy) * (Rhoh-Rhol)
!
!	END DO
!    END DO
!END


!SUBROUTINE Calculate_Viscous_Force_BGK( tau, DcDx, DcDy, gneq, FmX, FmY )
!!SUBROUTINE Calculate_Viscous_Force_BGK(tau_H, tau_L, DcDx, DcDy, gneq, C, FmX, FmY)
!  USE SHARE, ONLY: Rhoh, Rhol, C_L, C_H
!  USE, INTRINSIC :: IEEE_ARITHMETIC
!
!  IMPLICIT NONE
!  REAL(8),INTENT(IN)  :: DcDx, DcDy, gneq(0:8)
!!  REAL(8), INTENT(IN)  :: tau_H, tau_L, DcDx, DcDy, gneq(9), C
!  REAL(8), INTENT(OUT) :: FmX, FmY
!  REAL(8) :: sxx, sxy, syy, tau
!!  REAL(8) ::  tau, C_clamped
!
!!  C_local = MAX(0.0d0, MIN(1.0d0, C_local))
!!   C_clamped = MAX(0.0d0, MIN(1.0d0, C))
!!    IF (.NOT. IEEE_IS_FINITE(C)) THEN
!!    PRINT *, "NaN or Inf in C: ", C
!!    STOP
!!  END IF
!!  !tau = tau_L + (C(X,Y) - C_L) * (tau_H - tau_L) / (C_H - C_L)
!!  tau = tau_L + (C_clamped - C_L) * (tau_H - tau_L) / (C_H - C_L)
!!!  tau = 1.0d0 / ( (1.0d0 - C_local)/tau_L + C_local/tau_H )
!!
!
!  
!CALL Calculate_Stress_Tensor_BGK(gneq, sxx, sxy, syy)
!
!  FmX = (0.5d0 - tau) / tau * (sxx * DcDx + sxy * DcDy) * (Rhoh - Rhol)
!  FmY = (0.5d0 - tau) / tau * (sxy * DcDx + syy * DcDy) * (Rhoh - Rhol)
!
!!   IF (.NOT. IEEE_IS_FINITE(FmX) .OR. .NOT. IEEE_IS_FINITE(FmY)) THEN
!!    PRINT *, "NaN in viscous force! tau=", tau, "C=", C
!!    STOP
!!  END IF
!END SUBROUTINE


!**********************************************************************

!SUBROUTINE Calculate_Stress_Tensor_BGK( gneq, sxx, sxy, syy )
!	USE SHARE, ONLY: ex, ey
!	IMPLICIT NONE
!	REAL(8),INTENT(IN) :: gneq(1:8)
!	REAL(8),INTENT(OUT):: sxx, sxy, syy
!
!	sxx = SUM( gneq(1:) * ex(1:) * ex(1:) )
!	sxy = SUM( gneq(1:) * ex(1:) * ey(1:) )
!	syy = SUM( gneq(1:) * ey(1:) * ey(1:) )
!
!END

!**********************************************************************



SUBROUTINE RESULTS_Output
    USE SHARE, ONLY: X, Y, Nx, Ny, Ux, Uy, C, step, T, L0, P, Rho
    IMPLICIT NONE
    CHARACTER(LEN=16) :: filename
    CHARACTER(LEN=6) :: B

    
 ! Format timestep into 6-digit padded string
    WRITE(B, '(I6.6)') t
    filename = 'out/2D' // B // '.plt'

    OPEN(1, FILE=filename, STATUS='REPLACE')
    WRITE(1,*) 'Variables = X, Y, Ux, Uy, C, P, Rho'
    WRITE(1,*) 'Zone T = "', T, '" F=Point, I=', Nx, ', J=', Ny

!	WRITE(1,*) 'Zone T = "', T, '" F=Point, I=', Nx, ', J=', Ny

    DO Y = 1, Ny
    DO X = 1, Nx
        WRITE(1,'(2F8.3,5E14.6)') (X-0.5d0)/L0, (Y-0.5d0)/L0, Ux(X,Y), Uy(X,Y), C(X,Y)-0.5d0, P(X,Y)*Rho(X,Y)/3, Rho(X,Y)
    END DO
    END DO
    CLOSE(1)

END

! SUBROUTINE Interface_tracking
!	USE SHARE, ONLY: X, Y, Nx, Ny, T, L0, Rho, Rhoh, Rhol
!	IMPLICIT NONE
!	real*8 ::rho_d, pos_d, rho_u, pos_u, pos, pos2
!
!	WRITE(1,*) 'Zone T = "', T, '" F=Point, I=', Nx, ', J=', Ny
!	      open(45, file='posi.dat', access='append')
!      X = Nx / 2
!
!      do Y = 1, Ny
!        if (Rho(X, Y) .lt. (Rhoh + Rhol) / 2.d0) then
!            rho_d = Rho(X, Y)
!            pos_d = Y
!        endif
!      enddo
!
!      do y = Ny, 1, -1
!        if (Rho(X, Y) .gt. (Rhoh + Rhol) / 2.d0) then
!            rho_u = Rho(X, Y)
!            pos_u = Y
!        endif
!      enddo
!
!      pos = pos_d + (pos_u - pos_d) / (rho_u - rho_d) * ((Rhoh + Rhol)/ 2. - rho_d)
!      X = 1
!
!      do y = 1, Ny
!        if (Rho(X, Y) .lt. (Rhoh + Rhol) / 2.d0) then
!            rho_d = Rho(X, Y)
!            pos_d = y
!        endif
!      enddo
!
!      do y = Ny, 1, -1
!        if (Rho(X, Y)  .gt. (Rhoh + Rhol) / 2.d0) then
!            rho_u = Rho(X, Y)
!            pos_u = Y
!        endif
!      enddo
!
!             pos2 = pos_d + (pos_u - pos_d) / (rho_u - rho_d) * ((Rhoh + Rhol)/ 2. - rho_d )
!      write(45, '(2f15.4)')  pos, pos2
!
!
!END