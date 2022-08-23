module run_module
	use my_prec ! mp
	implicit none
	contains
	subroutine run_solver(M, D, x)     
			!solve Ax=D
			!input :: M - nx5 matrix
			!      :: D - nx1 vector
			!output:: x - nx1 vector
		real(mp), intent(in) :: M(:, :), D(:)
		real(mp), intent(out) :: x(:)
		real(mp), dimension(1:size(M, dim=1)) :: b, a, c 
		real(mp), dimension(1:size(M, dim=1)) :: beta, alpha, p, q, r 
		integer :: n, i
		n = size(x) - 1

! 		print*, M, shape(M)
! 		do i=1, n+1
! 			print*, M(i,:)
! 		enddo
		! extract a, b and c from M
		a(1) = M(1, 1); a(n) = M(n, 4)
		a(2) = M(2, 2); a(n+1) = M(n+1, 5)
		b(1) = M(1, 2); b(n) = M(n, 5)
		b(2) = M(2, 3); b(n+1) = 0
		c(1) = M(1, 3); c(n) = 0
		c(2) = M(2, 4); c(n+1) = 0
		do i=3, n-1
			a(i) = M(i, 3)
			b(i) = M(i, 4)
			c(i) = M(i, 5)
		enddo

! 		do i=1, n+1
! 			print*, a(i), b(i), c(i)
! 		enddo	

		p = 0; q = 0; r = 0
		beta(1) = 0; beta(2) = b(1)
		alpha(1) = a(1) 
		p(1) = b(1)/alpha(1)
		alpha(2) = a(2) - p(1)*beta(1)
		q(1) = c(1) / alpha(1); q(2) = c(2) / alpha(2)
		p(2) = (b(2) - q(1)*beta(2)) / alpha(2)
		r(1) = d(1) / alpha(1); r(2) = (d(2) - r(1)*beta(2)) / alpha(2)
		do i=3, n+1
			beta(i) = b(i-1) - p(i-2)*c(i-2)
			alpha(i) = a(i) - p(i-1)*beta(i) - q(i-2)*c(i-2)
			p(i) = (b(i) - q(i-1)*beta(i))/alpha(i)
			q(i) = c(i) / alpha(i)
			r(i) = (d(i) - r(i-1)*beta(i) - r(i-2)*c(i-2)) / alpha(i)
		enddo	

! 		print*, r
		x(n+1) = r(n+1)
		x(n) = r(n) - p(n)*x(n+1)
		do i=n-1,1,-1
			x(i) = r(i) - p(i)*x(i+1) - q(i)*x(i+2)
		enddo
! 		do i=1, n+1
! 			print*, x(i)
! 		enddo	

		return
	end subroutine run_solver
end module run_module