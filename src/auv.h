
interface
subroutine lagrangian(u_t, u_x, u_y, alpha, L)
implicit none
REAL*8, intent(in) :: u_t
REAL*8, intent(in) :: u_x
REAL*8, intent(in) :: u_y
REAL*8, intent(in) :: alpha
REAL*8, intent(out) :: L
end subroutine
end interface
interface
subroutine hamiltonian(x, y, z, v_x, v_y, v_z, q_r, q_x, q_y, q_z, &
lambda_x, lambda_y, lambda_z, lambda_v_x, lambda_v_y, &
lambda_v_z, lambda_q_r, lambda_q_x, lambda_q_y, lambda_q_z, &
u_t, u_x, u_y, alpha, v_xf, v_yf, v_zf, m, rho, cd, A, T, &
omega, H)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: v_x
REAL*8, intent(in) :: v_y
REAL*8, intent(in) :: v_z
REAL*8, intent(in) :: q_r
REAL*8, intent(in) :: q_x
REAL*8, intent(in) :: q_y
REAL*8, intent(in) :: q_z
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_z
REAL*8, intent(in) :: lambda_v_x
REAL*8, intent(in) :: lambda_v_y
REAL*8, intent(in) :: lambda_v_z
REAL*8, intent(in) :: lambda_q_r
REAL*8, intent(in) :: lambda_q_x
REAL*8, intent(in) :: lambda_q_y
REAL*8, intent(in) :: lambda_q_z
REAL*8, intent(in) :: u_t
REAL*8, intent(in) :: u_x
REAL*8, intent(in) :: u_y
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: v_xf
REAL*8, intent(in) :: v_yf
REAL*8, intent(in) :: v_zf
REAL*8, intent(in) :: m
REAL*8, intent(in) :: rho
REAL*8, intent(in) :: cd
REAL*8, intent(in) :: A
REAL*8, intent(in) :: T
REAL*8, intent(in) :: omega
REAL*8, intent(out) :: H
end subroutine
end interface
interface
subroutine pontryagin(x, y, z, v_x, v_y, v_z, q_r, q_x, q_y, q_z, &
lambda_x, lambda_y, lambda_z, lambda_v_x, lambda_v_y, &
lambda_v_z, lambda_q_r, lambda_q_x, lambda_q_y, lambda_q_z, &
alpha, v_xf, v_yf, v_zf, m, rho, cd, A, T, omega, u)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: v_x
REAL*8, intent(in) :: v_y
REAL*8, intent(in) :: v_z
REAL*8, intent(in) :: q_r
REAL*8, intent(in) :: q_x
REAL*8, intent(in) :: q_y
REAL*8, intent(in) :: q_z
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_z
REAL*8, intent(in) :: lambda_v_x
REAL*8, intent(in) :: lambda_v_y
REAL*8, intent(in) :: lambda_v_z
REAL*8, intent(in) :: lambda_q_r
REAL*8, intent(in) :: lambda_q_x
REAL*8, intent(in) :: lambda_q_y
REAL*8, intent(in) :: lambda_q_z
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: v_xf
REAL*8, intent(in) :: v_yf
REAL*8, intent(in) :: v_zf
REAL*8, intent(in) :: m
REAL*8, intent(in) :: rho
REAL*8, intent(in) :: cd
REAL*8, intent(in) :: A
REAL*8, intent(in) :: T
REAL*8, intent(in) :: omega
REAL*8, intent(out), dimension(1:3, 1:1) :: u
end subroutine
end interface
interface
subroutine eom_fullstate(x, y, z, v_x, v_y, v_z, q_r, q_x, q_y, q_z, &
lambda_x, lambda_y, lambda_z, lambda_v_x, lambda_v_y, &
lambda_v_z, lambda_q_r, lambda_q_x, lambda_q_y, lambda_q_z, &
u_t, u_x, u_y, v_xf, v_yf, v_zf, m, rho, cd, A, T, omega, &
dfs)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: v_x
REAL*8, intent(in) :: v_y
REAL*8, intent(in) :: v_z
REAL*8, intent(in) :: q_r
REAL*8, intent(in) :: q_x
REAL*8, intent(in) :: q_y
REAL*8, intent(in) :: q_z
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_z
REAL*8, intent(in) :: lambda_v_x
REAL*8, intent(in) :: lambda_v_y
REAL*8, intent(in) :: lambda_v_z
REAL*8, intent(in) :: lambda_q_r
REAL*8, intent(in) :: lambda_q_x
REAL*8, intent(in) :: lambda_q_y
REAL*8, intent(in) :: lambda_q_z
REAL*8, intent(in) :: u_t
REAL*8, intent(in) :: u_x
REAL*8, intent(in) :: u_y
REAL*8, intent(in) :: v_xf
REAL*8, intent(in) :: v_yf
REAL*8, intent(in) :: v_zf
REAL*8, intent(in) :: m
REAL*8, intent(in) :: rho
REAL*8, intent(in) :: cd
REAL*8, intent(in) :: A
REAL*8, intent(in) :: T
REAL*8, intent(in) :: omega
REAL*8, intent(out), dimension(1:20, 1:1) :: dfs
end subroutine
end interface
interface
subroutine eom_fullstate_jac(x, y, z, v_x, v_y, v_z, q_r, q_x, q_y, &
q_z, lambda_x, lambda_y, lambda_z, lambda_v_x, lambda_v_y, &
lambda_v_z, lambda_q_r, lambda_q_x, lambda_q_y, lambda_q_z, &
u_t, u_x, u_y, v_xf, v_yf, v_zf, m, rho, cd, A, T, omega, &
ddfs)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: v_x
REAL*8, intent(in) :: v_y
REAL*8, intent(in) :: v_z
REAL*8, intent(in) :: q_r
REAL*8, intent(in) :: q_x
REAL*8, intent(in) :: q_y
REAL*8, intent(in) :: q_z
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_z
REAL*8, intent(in) :: lambda_v_x
REAL*8, intent(in) :: lambda_v_y
REAL*8, intent(in) :: lambda_v_z
REAL*8, intent(in) :: lambda_q_r
REAL*8, intent(in) :: lambda_q_x
REAL*8, intent(in) :: lambda_q_y
REAL*8, intent(in) :: lambda_q_z
REAL*8, intent(in) :: u_t
REAL*8, intent(in) :: u_x
REAL*8, intent(in) :: u_y
REAL*8, intent(in) :: v_xf
REAL*8, intent(in) :: v_yf
REAL*8, intent(in) :: v_zf
REAL*8, intent(in) :: m
REAL*8, intent(in) :: rho
REAL*8, intent(in) :: cd
REAL*8, intent(in) :: A
REAL*8, intent(in) :: T
REAL*8, intent(in) :: omega
REAL*8, intent(out), dimension(1:20, 1:20) :: ddfs
end subroutine
end interface

