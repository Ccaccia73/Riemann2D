MODULE mod_exact_riemann

    USE mod_euler_flux_jacobian
    USE mod_thermodynamics

CONTAINS

    SUBROUTINE exact_riemann (wl, wr,  wCl, wCr, xdir)

        ! risolve il sistema di equazioni non lineari:

        !               P(vCl, wl)  =    P(vCr, wr)
        !             u_1(vCl, wl)  =  u_3(vCr, wr)
        !

        ! in realtà il sistema che risolviamo è:

        !             P(vCl, wl)  -    P(vCr, wr)  =  0
        !           u_1(vCl, wl)  -  u_3(vCr, wr)  =  0

        IMPLICIT NONE

        REAL(dp_kind), DIMENSION(:), INTENT(IN)  :: wl, wr ! stati destro e sinistro
        REAL(dp_kind), DIMENSION(:), INTENT(OUT) :: wCl, wCr ! 2 stati intermedi

        LOGICAL, INTENT(IN) :: xdir

        REAL(dp_kind) :: nu,  nu_vacuum,  &
            rhol, vl, ml, ul, el, cl, &
            rhor, vr, mr, ur, er, cr, &
            U2r, U2l, &
            trasv_ul, trasv_ur, &
            u1, v1, P1, Dv1, &
            u2, v2, P2, Dv2, &
            a11, a21, a12, a22, rhs1, rhs2, det

        REAL(dp_kind), PARAMETER :: rel_tol = 1.0d-10
        INTEGER, PARAMETER :: nmax = 1000

        INTEGER :: n

        ! risolviamo un pb di riemann in 1D usando u o v a seconda della direzione

        rhol = wl(1);   vl = 1/rhol
        rhor = wr(1);   vr = 1/rhor

        IF ( xdir) THEN
            ml = wl(2)
            mr = wr(2)
            trasv_ul = wl(3)*vl
            trasv_ur = wr(3)*vr
        ELSE
            ml = wl(3)
            mr = wr(3)
            trasv_ul = wl(2)*vl
            trasv_ur = wr(2)*vr
        END IF

        ul = ml * vl
        ur = mr * vr

        ! parametri limite
        ! due fan
        !    nu_2r =
        ! due shock
        !    nu_2s =

        ! modulo della velocità al quadrato

        U2l = ul**2 + trasv_ul**2
        U2r = ur**2 + trasv_ur**2

        el = wl(4) / rhol  -  U2l**2 / 2;   cl = sound_speed(el, rhol)
        er = wr(4) / rhor  -  U2r**2 / 2;   cr = sound_speed(er, rhor)

        nu_vacuum = 2 * (cl + cr) / (gamma -1)

        nu = ur - ul ! il confronto fra le velocità limite potrebbe
                     ! essere usato per scegliere il guess iniziale
                     ! lo faremo dopo.
                     ! per ora usiamo la media

        IF (nu >= nu_vacuum) THEN

            WRITE(*,*) 'formazione del vuoto'
            WRITE(*,*) 'STOP'
            STOP

        ENDIF

        ! guess iniziale media per entrambe le incognite del sistema 2x2
        v1 = (vl + vr) / 2
        v2 = (vl + vr) / 2

        DO n = 1, nmax

            ! rete di protezione per evitare di andare al di sotto
            ! del limite dell'adiabatica di RH del gas ideal politropico
            !v_asintoto = ((gamma-1)/(gamma+1)) * v_
            !maledetto gamma
            !questi controlli sono limitati al modello del PIG

            IF (v1 < ((gamma-1)/(gamma+1)) * vl) THEN

                v1 = (gamma/(gamma+1)) * vl

                WRITE(*,*) 'metodo di Newton nella subroutine'
                WRITE(*,*) 'exact_riemann genera un volume specifico v1'
                WRITE(*,*) 'inferiore al limite asintotico della adiabatica di RH'
                WRITE(*,*) 'correzione ad hoc forzata in modo ordinario'
                WRITE(*,*) 'v1 = (gamma/(gamma+1)) * vl'

            ENDIF

            CALL loci_rar_RH (1, wl, v1,  P1, u1, a11, a21, xdir)

            IF (v2 < ((gamma-1)/(gamma+1)) * vr) THEN

                v2 = (gamma/(gamma+1)) * vr

                WRITE(*,*) 'metodo di Newton nella subroutine'
                WRITE(*,*) 'exact_riemann genera un volume specifico v2'
                WRITE(*,*) 'inferiore al limite asintotico della adiabatica di RH'
                WRITE(*,*) 'correzione ad hoc forzata in modo ordinario'
                WRITE(*,*) 'v2 = (gamma/(gamma+1)) * vr'

            ENDIF

            CALL loci_rar_RH (3, wr, v2,  P2, u2, a12, a22, xdir)

            ! devo cambiare il segno di a11 e a22
            a12 = -a12;   a22 = -a22

            ! metodo di Newton in forma incrementale per sistema
            rhs1 = - (P1 - P2)
            rhs2 = - (u1 - u2)
            ! dobbiamo risolvere un sistema lineare 2x2

            det = a11 * a22  -  a12 * a21

            Dv1 = (a22*rhs1 - a12*rhs2) / det ! regola di Cramer
            Dv2 = (a11*rhs2 - a21*rhs1) / det

            v1 = v1 + Dv1
            v2 = v2 + Dv2

            IF (SQRT(Dv1**2 + Dv2**2) <=  rel_tol * SQRT(v1**2 + v2**2)) THEN

                wCl(1) = 1/v1;     wCr(1) = 1/v2

                CALL loci_rar_RH(1, wl, v1,  P1, u1, a11, a21, xdir)
                    ! non servono a11 e a21

                IF ( xdir ) THEN
                    wCl(2) = u1/v1;    wCr(2) = u1/v2  ! u2 == u1
                    wCl(3) = trasv_ul/v1;    wCR(3) =  trasv_ur/v2
                ELSE
                    wCl(2) = trasv_ul/v1;    wCR(2) =  trasv_ur/v2
                    wCl(3) = u1/v1;    wCr(3) = u1/v2  ! u2 == u1
                END IF

                wCl(4) = (specific_energy(P1, 1/v1)  +  (u1**2 + trasv_ul**2)/2) / v1 ! controllare se e(P, v) o e(P, rho)
                wCr(4) = (specific_energy(P1, 1/v2)  +  (u1**2 + trasv_ur**2)/2) / v2 ! P1 == P2


                WRITE(*,*) 'metodo di Newton converge in', n, 'iterazioni'
                WRITE(*,*)

                RETURN

            ENDIF

        END DO

        WRITE(*,*)
        WRITE(*,*) 'metodo di Newton non converge in', nmax, 'iterazioni'
        WRITE(*,*) 'STOP nella SUBROUTINE exact_riemann'
        STOP

    END SUBROUTINE exact_riemann


    SUBROUTINE loci_rar_RH (i, w_, v,  P, u, dP_dv, du_dv, xdir)

        IMPLICIT NONE

        INTEGER,                     INTENT(IN) :: i
        REAL(dp_kind), DIMENSION(:), INTENT(IN) :: w_  ! stato pivot
        REAL(dp_kind),               INTENT(IN) :: v ! var indipendente
        REAL(dp_kind),               INTENT(OUT) :: P, u, dP_dv, du_dv
        LOGICAL,                     INTENT(IN) :: xdir

        REAL(dp_kind) :: s, v_, u_, P_, DP, Dv, S_DD

        SELECT CASE (i)

            CASE(1);  s = -1

            CASE(3);  s = 1

            CASE DEFAULT

                WRITE(*,*) 'i deve essere 1 o 3 in loci_rar_RH'
                WRITE(*,*) 'STOP'
                STOP

        END SELECT

        ! definire energia di stato post discontinuità in funzione
        ! di volume specifico e pressione

        v_ = 1/w_(1)

        IF ( xdir ) THEN
            u_ = w_(2) * v_
        ELSE
            u_ = w_(3) * v_
        END IF

        P_ = Pgreco(w_)

        IF (v >= v_) THEN ! rarefazione

            P = pressure_RAR(P_, v_, v, dP_dv)

            u = u_  -  s * integral_c_v_isentropic(P_, v_, v)

            du_dv = -s * sound_speed_isentropic(P_, v_, v) / v

        ELSE ! urto

            P = pressure_RH(P_, v_, v, dP_dv)

            DP = P - P_;   Dv = v - v_

            S_DD = SQRT(-DP * Dv)

            u = u_  +  s * S_DD

            IF (S_DD /= 0) THEN
                du_dv = -s * (Dv * dP_dv  +  DP)/(2 * S_DD)
            ELSE
                du_dv = 0
            ENDIF

        END IF

    END SUBROUTINE loci_rar_RH

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !================================================================

    FUNCTION ws_sonic_state(i, w, xdir) RESULT(ws)

        ! Sonic values of the rarefaction wave
        ! similarity solution at xi = 0

        IMPLICIT NONE

        INTEGER,                     INTENT(IN)  :: i
        REAL(dp_kind), DIMENSION(:), INTENT(IN)  :: w
        LOGICAL,                     INTENT(IN)  :: xdir

        REAL(dp_kind), DIMENSION(SIZE(w)) :: ws

        REAL(KIND=8) :: s,  rho, u, trasv_u, e, P, c, lambda,  &
            rhos, us, Ets, Ps


        SELECT CASE(i)

            CASE(1);  s = -1
            CASE(4);  s = 1

            CASE DEFAULT

                WRITE (*,*) 'In FUNCTION ws_sonic_state, i'
                WRITE (*,*) 'must be either 1 or 4.  STOP.'
                STOP

        END SELECT


        rho = w(1)

        IF ( xdir ) THEN
           u = w(2)/rho
           trasv_u = w(3)/rho
        ELSE
           u = w(3)/rho
           trasv_u = w(2)/rho
        END IF



        e = w(4)/w(1) - (u**2 + trasv_u**2)/2

        P = Pgreco(w)

        c = sound_speed(e, rho)

        lambda = u  +  s * c

        rhos = rho * (1  -  s * ((gamma-1)/(gamma+1)) * lambda/c)**(2/(gamma-1))

        us = u  -  (2/(gamma+1)) * lambda

        Ps = P * (rhos/rho)**gamma

        Ets = rhos * (specific_energy(Ps, rhos)  +  (us**2 + trasv_u**2)/2)


        IF ( xdir ) THEN
            ws = [rhos, rhos*us, rhos*trasv_u, Ets]
        ELSE
            ws = [rhos, rhos*trasv_u, rhos*us, Ets]
        END IF


    END FUNCTION ws_sonic_state

!================================================================

END MODULE mod_exact_riemann
