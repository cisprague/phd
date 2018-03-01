
subroutine lagrangian(v_omega_1, v_omega_2, alpha, tpr, b, r, L)
implicit none
REAL*8, intent(in) :: v_omega_1
REAL*8, intent(in) :: v_omega_2
REAL*8, intent(in) :: alpha
REAL*8, intent(in) :: tpr
REAL*8, intent(in) :: b
REAL*8, intent(in) :: r
REAL*8, intent(out) :: L

L = alpha*(v_omega_1 + v_omega_2 + 1) + (-alpha + 1)*(v_omega_1**2 + &
      v_omega_2**2 + 1)

end subroutine

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

H = alpha*(v_omega_1 + v_omega_2 + 1) + lambda_x*((1.0d0/2.0d0)* &
      v_omega_1 + (1.0d0/2.0d0)*v_omega_2)*cos(theta) + (-alpha + 1)*( &
      v_omega_1**2 + v_omega_2**2 + 1) + (1.0d0/2.0d0)*lambda_theta*( &
      -v_omega_1 + v_omega_2)/b + (1.0d0/2.0d0)*lambda_y*(-v_omega_1 + &
      v_omega_2)*sin(theta)/b

end subroutine

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

u(1, 1) = (1.0d0/4.0d0)*(2*alpha*b + b*lambda_x*cos(theta) + &
      lambda_theta + lambda_y*sin(theta))/(b*(alpha - 1))
u(2, 1) = (1.0d0/4.0d0)*(2*alpha*b + b*lambda_x*cos(theta) - &
      lambda_theta - lambda_y*sin(theta))/(b*(alpha - 1))

end subroutine

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

dfs(1, 1) = ((1.0d0/2.0d0)*v_omega_1 + (1.0d0/2.0d0)*v_omega_2)*cos( &
      theta)
dfs(2, 1) = (1.0d0/2.0d0)*(-v_omega_1 + v_omega_2)*sin(theta)/b
dfs(3, 1) = (1.0d0/2.0d0)*(-v_omega_1 + v_omega_2)/b
dfs(4, 1) = 0
dfs(5, 1) = 0
dfs(6, 1) = lambda_x*((1.0d0/2.0d0)*v_omega_1 + (1.0d0/2.0d0)*v_omega_2) &
      *sin(theta) - 1.0d0/2.0d0*lambda_y*(-v_omega_1 + v_omega_2)*cos( &
      theta)/b

end subroutine

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

ddfs(1, 1) = 0
ddfs(2, 1) = 0
ddfs(3, 1) = 0
ddfs(4, 1) = 0
ddfs(5, 1) = 0
ddfs(6, 1) = 0
ddfs(1, 2) = 0
ddfs(2, 2) = 0
ddfs(3, 2) = 0
ddfs(4, 2) = 0
ddfs(5, 2) = 0
ddfs(6, 2) = 0
ddfs(1, 3) = -((1.0d0/2.0d0)*v_omega_1 + (1.0d0/2.0d0)*v_omega_2)*sin( &
      theta)
ddfs(2, 3) = (1.0d0/2.0d0)*(-v_omega_1 + v_omega_2)*cos(theta)/b
ddfs(3, 3) = 0
ddfs(4, 3) = 0
ddfs(5, 3) = 0
ddfs(6, 3) = lambda_x*((1.0d0/2.0d0)*v_omega_1 + (1.0d0/2.0d0)*v_omega_2 &
      )*cos(theta) + (1.0d0/2.0d0)*lambda_y*(-v_omega_1 + v_omega_2)* &
      sin(theta)/b
ddfs(1, 4) = 0
ddfs(2, 4) = 0
ddfs(3, 4) = 0
ddfs(4, 4) = 0
ddfs(5, 4) = 0
ddfs(6, 4) = ((1.0d0/2.0d0)*v_omega_1 + (1.0d0/2.0d0)*v_omega_2)*sin( &
      theta)
ddfs(1, 5) = 0
ddfs(2, 5) = 0
ddfs(3, 5) = 0
ddfs(4, 5) = 0
ddfs(5, 5) = 0
ddfs(6, 5) = -1.0d0/2.0d0*(-v_omega_1 + v_omega_2)*cos(theta)/b
ddfs(1, 6) = 0
ddfs(2, 6) = 0
ddfs(3, 6) = 0
ddfs(4, 6) = 0
ddfs(5, 6) = 0
ddfs(6, 6) = 0

end subroutine
