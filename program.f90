PROGRAM tubePrandtlNavierStokes
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureDensityCoupling
    INTEGER(2), PARAMETER :: io = 12
    INTEGER(4) :: ni, nj
    INTEGER(4) :: i, j, s_max
    REAL(8) :: l, h, dx, dy, u_0, rho_0, mu, rho, eps, dt, c, gamma, cfl, a
    REAL(8), ALLOCATABLE :: x_node(:,:), y_node(:,:), x_cell(:,:), y_cell(:,:)
    REAL(8), ALLOCATABLE :: u_n(:,:), v_n(:,:), p_n(:,:), rho_n(:,:)
    REAL(8), ALLOCATABLE :: u_c(:,:), v_c(:,:), p_c(:,:), rho_c(:,:)


    CALL DataInput(io, l, h, ni, nj, u_0, rho_0, mu, s_max, eps, cfl, c, gamma)

    ALLOCATE(x_node(ni,nj))
    ALLOCATE(y_node(ni,nj))
    ALLOCATE(x_cell(0:ni, 0:nj))
    ALLOCATE(y_cell(0:ni, 0:nj))

    ALLOCATE(u_n(ni,nj))   
    ALLOCATE(v_n(ni,nj))   
    ALLOCATE(p_n(ni,nj))
    ALLOCATE(rho_n(ni,nj))

    ALLOCATE(u_c(0:ni,0:nj))   
    ALLOCATE(v_c(0:ni,0:nj))   
    ALLOCATE(p_c(0:ni,0:nj))
    ALLOCATE(rho_c(0:ni,0:nj))

    CALL MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell, u_0, rho_0, a, cfl, dt, mu)

    CALL InitialConditionsPrandtl(ni, nj, u_0, rho_0, u_n, v_n, p_n, rho_n, c, gamma)

    CALL InitialConditionsNavierStokes(ni, nj, u_0, rho_0, u_c, v_c, p_c, rho_c, c, gamma)

    CALL SolverPrandtl(ni, nj, s_max, dx, dy, mu, eps, u_n, v_n, rho_n, p_n, u_0, rho_0, h, gamma, c)

    CALL SolverNavierStokes(ni, nj, s_max, dx, dy, mu, eps, u_0, rho_0, u_c, v_c, p_c, rho_c, dt, io, c, gamma, a)

    CALL OutputFieldsNode(io, ni, nj, x_node, y_node, u_n, v_n, p_n - &
        PressureDensityCoupling(c, rho_n(1,1), gamma), rho_n - rho_0, c, gamma)

    CALL OutputFieldsCell(io, ni, nj, x_cell, y_cell, u_c, v_c, p_c - &
        PressureDensityCoupling(c, rho_n(1,1), gamma), rho_c - rho_0, c, gamma)

    DEALLOCATE(x_node)
    DEALLOCATE(y_node)
    DEALLOCATE(x_cell)
    DEALLOCATE(y_cell)
    
    DEALLOCATE(u_n)   
    DEALLOCATE(v_n)   
    DEALLOCATE(p_n)
    DEALLOCATE(rho_n)

    DEALLOCATE(u_c)   
    DEALLOCATE(v_c)   
    DEALLOCATE(p_c)
    DEALLOCATE(rho_c)


END PROGRAM


SUBROUTINE DataInput(io, l, h, ni, nj, u_0, rho_0, mu, s_max, eps, cfl, c, gamma)
    ! Takes input data from file input.txt
    IMPLICIT NONE
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj, s_max
    REAL(8) :: l, h, u_0, rho_0, mu, eps, cfl, c, gamma
    INTENT(IN) io
    INTENT(OUT) l, h, ni, nj, u_0, rho_0, mu, s_max, eps, cfl, c, gamma

    WRITE(*,*) 'READING INPUT FILE'
    OPEN(io,FILE='INPUT.TXT')
    READ(io,*) l
    READ(io,*) h
    READ(io,*) ni
    READ(io,*) nj
    READ(io,*) u_0
    READ(io,*) rho_0
    READ(io,*) mu
    READ(io,*) s_max
    READ(io,*) eps
    READ(io,*) cfl
    READ(io,*) c
    READ(io,*) gamma

    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE MeshMaking(ni, nj, l, h, dx, dy, x_node, y_node, x_cell, y_cell, u_0, rho_0, a, cfl, dt, mu)
    ! Makes mesh for numerical solution Prandtl (node) and 
    ! Navier-Stokes (cell) systems of equations
    IMPLICIT NONE
    INTEGER(4) :: ni, nj, i, j
    REAL(8) :: l, h, dx, dy, u_0, rho_0, a, cfl, dt, mu
    REAL(8), DIMENSION(ni,nj) :: x_node, y_node
    REAL(8), DIMENSION(0:ni,0:nj) :: x_cell, y_cell
    INTENT(IN) l, h, ni, nj
    INTENT(OUT) dx, dy, x_node, y_node, x_cell, y_cell

    WRITE(*,*) 'MESH MAKING'

    dx = l / (ni - 1)
    dy = h / (nj - 1)
    dt = cfl * MIN(0.5D0 * dx * dx * rho_0 / mu, 0.5D0 * dy * dy * rho_0 / mu, dx / u_0)
    a = 1D0 / (u_0 * u_0)

    DO i = 1, ni
        DO j = 1, nj
            x_node(i,j) = (i - 1) * dx
            y_node(i,j) = (j - 1) * dy
        END DO
    END DO

    DO i = 0, ni
        DO j = 0, nj
            x_cell(i,j) = (i - 5D-1) * dx
            y_cell(i,j) = (j - 5D-1) * dy
        END DO
    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE InitialConditionsPrandtl(ni, nj, u_0, rho_0, u, v, p, rho, c, gamma)
    ! Initial uniform velocity, density and pressure conditions 
    ! in the inlet
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureDensityCoupling
    INTEGER(4) :: i, j, ni, nj
    REAL(8) :: u_0, rho_0, c, gamma
    REAL(8), DIMENSION(ni,nj) :: u, p, rho, v
    INTENT(IN) ni, nj, u_0, rho_0, c, gamma
    INTENT(OUT) u, p, rho

    WRITE(*,*) 'INITIAL CONDITIONS APPLYING (PRANDTL)'

    rho(1,:) = rho_0

    p(1,:) = PressureDensityCoupling(c, rho_0, gamma)

    u(1,:) = u_0

    v(1,:) = 0D0
    
    WRITE(*,*) 'SUCCESS'
    
    END SUBROUTINE


SUBROUTINE BoundaryConditionsNavierStokes(ni, nj, u_0, rho_0, u, v, p, rho, c, gamma)
    ! Boundary conditions for velocities and pressure
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureDensityCoupling
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0, rho_0, c, gamma
    REAL(8), DIMENSION(0:ni,0:nj) :: u, v, p, rho
    INTENT(IN) ni, nj, u_0, rho_0
    INTENT(OUT) u, v, p, rho
    
    u(0, 1:nj) = u_0
    v(0, 1:nj) = 0D0
    p(0, 1:nj) =  p(1, 1:nj)
    rho(0, 1:nj) = rho(1, 1:nj)

    u(ni, 1:nj) = u(ni - 1, 1:nj)
    v(ni, 1:nj) = v(ni - 1, 1:nj)
    rho(ni, 1:nj) = rho_0 
    p(ni, 1:nj) = PressureDensityCoupling(c, rho_0, gamma)

    u(1:ni, 0) = - u(1:ni, 1)
    v(1:ni, 0) = - v(1:ni, 1)
    rho(1:ni, 0) = rho(1:ni, 1)
    p(1:ni, 0) = p(1:ni, 1)

    u(1:ni, nj) = u(1:ni, nj - 1)
    v(1:ni, nj) = - v(1:ni, nj - 1)
    rho(1:ni, nj) = rho(1:ni, nj - 1)
    p(1:ni, nj) = p(1:ni, nj - 1)
   
    END SUBROUTINE


SUBROUTINE InitialConditionsNavierStokes(ni, nj, u_0, rho_0, u, v, p, rho, c, gamma)
    ! Initial velocities and pressure conditions in area
    IMPLICIT NONE
    REAL(8), EXTERNAL :: PressureDensityCoupling
    INTEGER(4) :: ni, nj
    REAL(8) :: u_0, rho_0, c, gamma
    REAL(8), DIMENSION(0:ni,0:nj) :: u, v, p, rho
    INTENT(IN) ni, nj, u_0
    INTENT(OUT) u, v, p

    WRITE(*,*) 'INITIAL CONDITIONS APPLYING (NAVIER-STOKES)'
    
    u = u_0

    v = 0D0

    rho = rho_0       

    p = PressureDensityCoupling(c, rho_0, gamma)   

    WRITE(*,*) 'SUCCESS'
    
    END SUBROUTINE


SUBROUTINE ThomasAlgorithm(k_max, a, b, c, d, res)
    ! Solution of tridiagonal system 
    ! a_{k} * x_{k - 1} + b_{k} * x_{k} + c_{k} * x_{k + 1} = d_{k}
    ! with k_max unknowns (a_{1} = 0, c_{k_max} = 0)
    IMPLICIT NONE
    INTEGER(4) :: k, k_max
    REAL(8), DIMENSION(k_max) :: a, b, c, d, alpha, beta, res
    INTENT(IN) k_max, a, b, c, d
    INTENT(OUT) res
   
    alpha(2) = - c(1) / b(1)
    beta(2) = d(1) / b(1)

    DO k = 2, k_max - 1
        alpha(k + 1) = - c(k) / (b(k) + a(k) * alpha(k))
        beta(k + 1) = (d(k) - a(k) * beta(k)) / (b(k) + a(k) * alpha(k))
    END DO

    res(k_max) = (d(k) - a(k) * beta(k)) / (b(k) + a(k) * alpha(k))

    DO k = k_max - 1, 1, -1
        res(k) = alpha(k + 1) * res(k + 1) + beta(k + 1)
    END DO

    END SUBROUTINE


SUBROUTINE SolverPrandtl(ni, nj, s_max, dx, dy, mu, eps, u, v, rho, p, u_0, rho_0, h, gamma, const)
    ! Solver for Prandtl (Simplified Navier-Stokes) equations system
    IMPLICIT NONE
    LOGICAL(1), EXTERNAL :: ConvergenceCheckPrandtl
    REAL(8), EXTERNAL :: ResidualPrandtl
    REAL(8), EXTERNAL :: DensityPressureCoupling
    INTEGER(4) :: i, j, s, ni, nj, s_max, k
    REAL(8) :: dx, dy, mu, eps, g_0, g, h, u_0, gamma, const, rho_0
    REAL(8), DIMENSION(ni,nj) :: u, v, p, rho, u_temp, v_temp, p_temp, rho_temp
    REAL(8), DIMENSION(nj) :: a, b, c, d
    INTENT(IN) :: ni, nj, s_max, dx, dy, mu, eps
    INTENT(OUT) :: u, v

    WRITE(*,*) 'SOLVING EQUATIONS (PRANDTL)'

    g_0 = rho_0 * u_0 * h

    u_temp = 0D0
    v_temp = 0D0

    do i = 2, ni

        u(i,:) = u(i - 1,:)
        v(i,:) = v(i - 1,:)
        p(i,:) = p(i - 1,:)
        rho(i,:) = rho(i - 1,:)

        DO k = 1, 1000

            u(i,:) = u(i - 1,:)
            v(i,:) = v(i - 1,:)

            DO s = 1, s_max

                a(1) = 0D0
                b(1) = 1D0
                c(1) = 0D0
                d(1) = 0D0

                DO j = 2, nj - 1
                    a(j) = - rho(i, j - 1) * v(i, j - 1) / (2D0 * dy) - mu / dy**2
                    b(j) = rho(i, j) * u(i, j) / dx + 2D0 * mu / dy**2
                    c(j) = rho(i, j + 1) * v(i, j + 1) / (2D0 * dy) - mu / dy**2
                    d(j) = rho(i - 1, j) * u(i - 1, j)**2 / dx - (p(i,j) - p(i - 1, j)) / dx
                END DO

                a(nj) = -1D0
                b(nj) = 1D0
                c(nj) = 0D0
                d(nj) = 0D0

                CALL ThomasAlgorithm(nj, a, b, c, d, u(i, :))

                DO j = 2, nj
                    v(i,j) = v(i, j - 1) - dy / (2D0 * dx) * (u(i, j) + u(i, j - 1) - &
                        rho(i - 1, j) / rho(i,j) * (u(i - 1,j) + u(i - 1, j - 1)))
                END DO

                IF ((ConvergenceCheckPrandtl(u(i, :), u_temp(i, :), nj, eps)) .AND. &
                        (ConvergenceCheckPrandtl(v(i, :), v_temp(i, :), nj, eps))) THEN
                    ! WRITE(*,*) 'SOLUTION CONVERGED BY RESIDUALS, NODE №', I, ', s = ', s
                    EXIT 
                END IF

                ! IF (s == s_max) THEN
                    ! WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY, NODE №', I
                ! END IF

                v_temp = v
                u_temp = u

            END DO

            g = 0D0
            DO j = 2, nj
                g = g + u(i,j) * rho(i,j)
            END DO
            g = g * dy

            write(*,*) i, k, g, p(i, 1)
            
            IF (ABS(g_0 - g) / g_0 < 1D-6) THEN
                EXIT
            END IF

            p(i,:) = p(i,:) + 1d-1 * g_0 / h**2 * (g - g_0)

        END DO

    END DO

    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE


SUBROUTINE SolverNavierStokes(ni, nj, s_max, dx, dy, mu, eps, u_0, rho_0, u, v, p, rho, dt, io, c, gamma, a)
    !Solver for Navier-Stokes system of equations
    IMPLICIT NONE
    REAL(8), EXTERNAL :: HalfIndexValue, ResidualNavierStokes, DensityPressureCoupling
    INTEGER(2) :: io
    INTEGER(4) :: i, j, ni, nj, s, s_max
    REAL(8) :: u_hat_left, u_hat_right
    REAL(8) :: v_hat_top, v_hat_bot
    REAL(8) :: u_left, u_right, u_top, u_bot
    REAL(8) :: v_left, v_right, v_top, v_bot
    REAL(8) :: rho_left, rho_right, rho_top, rho_bot
    REAL(8) :: p_left, p_right, p_top, p_bot
    REAL(8) :: dx, dy, mu, eps, u_0, a, dt, res_u, res_v, res_p, c, gamma, rho_0
    REAL(8), DIMENSION(0:ni,0:nj) :: u_old, v_old, p_old, rho_old, u, v, p, rho
    INTENT(IN) ni, nj, s_max, dx, dy, mu, eps, u_0, io
    INTENT(OUT) u, v, p

    WRITE(*,*) 'SOLVING EQUATIONS (NAVIER-STOKES)'

    OPEN(io, FILE='RESIDUAL_LOGS.plt')
    WRITE(io,*) 'VARIABLES = "S", "RES_U", "RES_V", "RES_P"'

    DO s = 1, s_max

        u_old = u
        v_old = v
        p_old = p
        rho_old = rho

        DO i = 1, ni - 1
            DO j = 1, nj - 1

                u_hat_left = 5D-1 * (u(i - 1, j) + u(i,j))
                u_hat_right = 5D-1 * (u(i,j) + u(i + 1,j))
                
                v_hat_top = 5D-1 * (v(i,j) + v(i,j + 1))
                v_hat_bot = 5D-1 * (v(i,j - 1) + v(i,j))

                u_left = HalfIndexValue(u_hat_left, u(i, j), u(i - 1,j))
                v_left = HalfIndexValue(u_hat_left, v(i, j), v(i - 1,j))
                rho_left = HalfIndexValue(u_hat_left, rho(i - 1, j), rho(i,j))
                p_left = HalfIndexValue(u_hat_left, p(i - 1, j), p(i,j))
                
                u_right = HalfIndexValue(u_hat_right, u(i + 1, j), u(i,j))
                v_right = HalfIndexValue(u_hat_right, v(i + 1, j), v(i,j))
                rho_right = HalfIndexValue(u_hat_right, rho(i, j), rho(i + 1,j))
                p_right = HalfIndexValue(u_hat_right, p(i, j), p(i + 1,j))
                
                u_top = HalfIndexValue(v_hat_top, u(i, j + 1), u(i,j))
                v_top = HalfIndexValue(v_hat_top, v(i, j + 1), v(i,j))
                rho_top = HalfIndexValue(v_hat_top, rho(i, j), rho(i,j + 1))
                p_top = HalfIndexValue(v_hat_top, p(i, j), p(i,j + 1))
                
                u_bot = HalfIndexValue(v_hat_bot, u(i, j), u(i,j - 1))
                v_bot = HalfIndexValue(v_hat_bot, v(i, j), v(i,j - 1))
                rho_bot = HalfIndexValue(v_hat_bot, rho(i, j - 1), rho(i,j))
                p_bot = HalfIndexValue(v_hat_bot, p(i, j - 1), p(i,j))

                IF (j == nj - 1) THEN
                    p(i,j) = p(i,j) - dt / a * ((rho_right * u_right - rho_left * u_left) / dx &
                        - rho_bot * v_bot / dy)
                ELSE
                    p(i,j) = p(i,j) - dt / a * ((rho_right * u_right - rho_left * u_left) / dx &
                        + (rho_top * v_top - rho_bot * v_bot) / dy)
                END IF

                ! WRITE(*,*) 'rho'
                ! WRITE(*,*) rho(i,j)                

                ! WRITE(*,*) 'p'
                ! WRITE(*,*) p(i,j)

                u(i,j) = u(i,j) - dt / rho(i,j) * ((u_hat_right * u_right * rho_right - u_hat_left * u_left * rho_left) / dx &
                    + (v_hat_top * u_top * rho_top - v_hat_bot * u_bot * rho_bot) / dy &
                    + (p_right - p_left) / dx &
                    - 4D0 / 3D0 * mu * (u(i + 1, j) - 2D0 * u(i,j) + u(i - 1, j)) / dx ** 2 &
                    - mu * (u(i, j + 1) - 2D0 * u(i,j) + u(i, j - 1)) / dy ** 2) 

                ! WRITE(*,*) 'u'
                ! WRITE(*,*) u(i,j)

                v(i,j) = v(i,j) - dt / rho(i,j) * ((u_hat_right * v_right * rho_right - u_hat_left * v_left * rho_left) / dx &
                    + (v_hat_top * v_top * rho_top - v_hat_bot * v_bot * rho_bot) / dy &
                    + (p_top - p_bot) / dy &
                    - mu * (v(i, j + 1) - 2D0 * v(i,j) + v(i, j - 1)) / dy ** 2) 
                    - 4D0 / 3D0 * mu * (v(i + 1, j) - 2D0 * v(i,j) + v(i - 1, j)) / dx ** 2)

                rho(i,j) = DensityPressureCoupling(c, p(i,j), gamma)
                ! WRITE(*,*) 'i, j, rho, p, u, v'
                ! WRITE(*,*) i, j, rho(i,j), p(i,j), u(i,j), v(i,j)

            END DO
        END DO

        CALL BoundaryConditionsNavierStokes(ni, nj, u_0, rho_0, u, v, p, rho, c, gamma)

        res_u = ResidualNavierStokes(u(1:ni - 1, 1:nj - 1), u_old(1:ni - 1, 1:nj - 1), ni - 1, nj - 1, dt)
        res_v = ResidualNavierStokes(v(1:ni - 1, 1:nj - 1), v_old(1:ni - 1, 1:nj - 1), ni - 1, nj - 1, dt)
        res_p = ResidualNavierStokes(p(1:ni - 1, 1:nj - 1), p_old(1:ni - 1, 1:nj - 1), ni - 1, nj - 1, dt)
        
        IF (s > 2) THEN
            WRITE(io,*) s, res_u, res_v, res_p
        END IF

        IF ((res_u < eps) .AND. (res_v < eps) .AND. (res_p < eps) .AND. (s > 1)) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY RESIDUALS'
            EXIT
        END IF

        IF (MOD(s,100) == 0) THEN
            WRITE(*,*) s, 'ITERATIONS MADE'
            WRITE(*,*) res_u, res_v, res_p
        END IF

        IF (s == s_max) THEN
            WRITE(*,*) 'SOLUTION CONVERGED BY ITERATIONS BOUNDARY'
        END IF

    END DO

    CLOSE(io)

    WRITE(*,*) 'SUCCESS'



    END SUBROUTINE


SUBROUTINE OutputFieldsNode(io, ni, nj, x, y, u, v, p, rho, c, gamma)
    ! Nodes-based results output
    IMPLICIT NONE
    REAL(8), EXTERNAL :: DensityPressureCoupling
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8) :: c, gamma
    REAL(8), DIMENSION(ni,nj) :: x, y
    REAL(8), DIMENSION(ni,nj) :: u, v, p, rho
    INTENT(IN) io, ni, nj, x, y, u, v, p, rho
    
    WRITE(*,*) 'RESULTS OUTPUT (PRANDTL) ' 
    OPEN(io,FILE='RES_PR.PLT')
    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P", "RHO"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') u(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') v(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') p(1:ni, 1:nj)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') rho(1:ni, 1:nj)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE 


SUBROUTINE OutputFieldsCell(io, ni, nj, x, y, u, v, p, rho, c, gamma)
    ! Cells-based results output
    IMPLICIT NONE
    REAL(8), EXTERNAL :: DensityPressureCoupling
    INTEGER(2) :: io
    INTEGER(4) :: ni, nj
    REAL(8) :: c, gamma
    REAL(8), DIMENSION(0:ni, 0:nj) :: x, y
    REAL(8), DIMENSION(0:ni, 0:nj) :: u, v, p, rho
    INTENT(IN) io, ni, nj, x, y, u, v, p, rho
        
    WRITE(*,*) 'RESULTS OUTPUT (NAVIER-STOKES)' 
    OPEN(io, FILE='RES_NS.PLT')
    WRITE(io,*) 'VARIABLES = "X", "Y", "U", "V", "P", "RHO"' 
    WRITE(io,*) 'ZONE I=', ni, ', J=', nj, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
    WRITE(io,'(100E25.16)') x(1:ni, 1:nj) 
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') y(1:ni, 1:nj)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') u(1:ni - 1, 1:nj - 1)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') v(1:ni - 1, 1:nj - 1)
    WRITE(io,*) 
    WRITE(io,'(100E25.16)') p(1:ni - 1, 1:nj - 1)
    WRITE(io,*)
    WRITE(io,'(100E25.16)') rho(1:ni - 1, 1:nj - 1)
    CLOSE(io)
    WRITE(*,*) 'SUCCESS'

    END SUBROUTINE 


REAL(8) FUNCTION HalfIndexValue(arg, minus_res, plus_res)
    !Upwind scheme direction definition 
    IMPLICIT NONE
    REAL(8) :: arg, minus_res, plus_res

    IF (arg < 0D0) THEN
        HalfIndexValue = minus_res
    ELSE
        HalfIndexValue = plus_res
    END IF

    END FUNCTION


LOGICAL(1) FUNCTION ConvergenceCheckPrandtl(a, b, n, eps)
    !Convergence check for iteration process of solving
    !system of Prandtl equations
    IMPLICIT NONE
    REAL(8), EXTERNAL :: ResidualPrandtl
    REAL(8), DIMENSION(n) :: a, b
    REAL(8) :: eps
    INTEGER(4) :: n

    ConvergenceCheckPrandtl = (ResidualPrandtl(a, b, n) < eps)

    END FUNCTION


REAL(8) FUNCTION ResidualPrandtl(a, b, n)
    !Calculation residuals for Prandtl system of equations
    IMPLICIT NONE
    REAL(8), DIMENSION(n) :: a, b, dif
    INTEGER(4) :: i, n

    DO i = 1, n
        dif(i) = abs(a(i) - b(i)) / abs(a(i))
    END DO

    ResidualPrandtl = MAXVAL(dif)

    END FUNCTION


REAL(8) FUNCTION PressureDensityCoupling(c, x, gamma)
    ! Calculation pressure from density using barotropic law
    IMPLICIT NONE
    REAL(8) :: c, x, gamma

    PressureDensityCoupling = c * x**gamma

    END FUNCTION


REAL(8) FUNCTION DensityPressureCoupling(c, x, gamma)
    ! Calculation density from pressure using barotropic law
    IMPLICIT NONE
    REAL(8) :: c, x, gamma

    DensityPressureCoupling = (x / c)**(1D0 / gamma)

    END FUNCTION


REAL(8) FUNCTION ResidualNavierStokes(a, b, ni, nj, dt)
    !Calculation residuals for Navier-Stokes system of equations
    IMPLICIT NONE
    REAL(8), DIMENSION(ni, nj) :: a, b, dif
    REAL(8) :: dt
    INTEGER(4) :: i, j, ni, nj

    DO i = 1, ni
        DO j = 1, nj
            dif(i,j) = abs(a(i,j) - b(i,j))
        END DO
    END DO

    ResidualNavierStokes = MAXVAL(dif) / dt 

    END FUNCTION
