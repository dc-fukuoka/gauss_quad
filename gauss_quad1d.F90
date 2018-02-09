! 1D
module params
  implicit none
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: npts = 2
  real(dp),dimension(npts),parameter :: weights = (/1.0d0, 1.0d0/)
  real(dp),dimension(npts),parameter :: gauss_nodes = (/-1.0d0/sqrt(3.0d0), 1.0d0/sqrt(3.0d0)/)
end module params

module subs
  use params
  implicit none
  abstract interface
     function func(alpha, beta, gamma, x) result(res)
       implicit none
       real(8),intent(in) :: alpha, beta, gamma, x
       real(8) :: res
     end function func
  end interface
  
contains
  function f(alpha, beta, gamma, x) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, gamma, x
    real(dp) :: res
    
    res = alpha*x**2+beta*x+gamma
  end function f
  function ff(alpha, beta, gamma, x) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, gamma, x
    real(dp) :: res
    
    res = alpha*x**3/3.0d0+beta*x**2/2.0d0+gamma*x
  end function ff

  function analytic(alpha, beta, gamma, a, b) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, gamma
    real(dp),intent(in) :: a, b
    procedure(func),pointer :: pf
    real(dp) :: res

    pf => ff
    res = pf(alpha, beta, gamma, b) - pf(alpha, beta, gamma, a)
    
  end function analytic

  function gauss_quad_1d(alpha, beta,gamma, a, b, pf) result(res)
    implicit none
    real(dp),intent(in) :: alpha, beta, gamma
    real(dp),intent(in) :: a, b
    procedure(func),pointer,intent(in) :: pf
    real(dp) :: res
    real(dp) :: x
    integer :: i

    res = 0.0d0
    do i = 1, npts
       x = (b-a)/2.0d0*gauss_nodes(i) + (b+a)/2.0d0
       res = res + (b-a)/2.0d0*weights(i)*(pf(1.0d0, 1.0d0, 1.0d0, x))
    end do
  end function gauss_quad_1d
end module subs

program main
  use subs
  implicit none
  procedure(func),pointer :: pf => null()
  real(dp) :: alpha, beta, gamma
  real(dp) :: a, b
  real(dp),dimension(2) :: sol

  pf => f
  
  alpha = 1.0d0
  beta  = 1.0d0
  gamma = 1.0d0

  a = -4.0d0
  b =  3.0d0
  write(6, *) "integration alpha*x^2 + beta*x + gamma from a to b"
  write(6, '(a, 3(1pe14.5))') "alpha, beta, gamma:", alpha, beta, gamma
  write(6, '(a, 2(1pe14.5))') "a, b:", a, b
  sol(1) = analytic(alpha, beta, gamma, a, b)
  sol(2) = gauss_quad_1d(alpha, beta, gamma, a, b, pf)
  write(6, '(a, 1pe14.5)') "analytic solution  :", sol(1)
  write(6, '(a, 1pe14.5)') "gaussian quadrature:", sol(2)
  write(6, '(a, 1pe14.5)') "diff:", sol(2) - sol(1)
  
  stop
end program main
  
