
interface
subroutine lagrangian(v_omega_1, v_omega_2, alpha, tpr, b, r, L)
implicit none
REAL*8, intent(in) :: v_omega_1
REAL*8, intent(in) :: v_omega_2
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out) :: L
end subroutine
end interface
interface
subroutine hamiltonian(x, y, theta, lambda_x, lambda_y, lambda_theta, &
v_omega_1, v_omega_2, alpha, tpr, b, r, H)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_theta
REAL*8, intent(in) :: v_omega_1
REAL*8, intent(in) :: v_omega_2
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out) :: H
end subroutine
end interface
interface
subroutine pontryagin(x, y, theta, lambda_x, lambda_y, lambda_theta, &
alpha, tpr, b, r, u)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_theta
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out), dimension(1:2, 1:1) :: u
end subroutine
end interface
interface
subroutine eom_fullstate(x, y, theta, lambda_x, lambda_y, lambda_theta, &
v_omega_1, v_omega_2, tpr, b, r, dfs)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_theta
REAL*8, intent(in) :: v_omega_1
REAL*8, intent(in) :: v_omega_2
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out), dimension(1:6, 1:1) :: dfs
end subroutine
end interface
interface
subroutine eom_fullstate_jac(x, y, theta, lambda_x, lambda_y, &
lambda_theta, v_omega_1, v_omega_2, tpr, b, r, ddfs)
implicit none
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: lambda_x
REAL*8, intent(in) :: lambda_y
REAL*8, intent(in) :: lambda_theta
REAL*8, intent(in) :: v_omega_1
REAL*8, intent(in) :: v_omega_2
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out), dimension(1:6, 1:6) :: ddfs
end subroutine
end interface

