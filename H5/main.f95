program main
	use my_prec ! mp
	use spline_module ! f, f_derivative
	implicit none
	real(mp), allocatable, dimension(:, :) :: input
	real(mp), allocatable, dimension(:) :: X, Y, P
	real(mp), allocatable, dimension(:) :: x_res, y_res
	integer :: n, i

	open(unit=1, file='data.dat', action="read")
		read(1,'(2X, I8)') n 
		allocate(input(3, n+1))
		read(1, *) input;
		allocate(X(0:n), Y(0:n), P(0:n))
		allocate(x_res(0:100*n-1), y_res(0:100*n-1))
	close(1)

	X = input(1, :)
	Y = input(2, :)
	P = input(3, :)

	call f(X, Y, P, x_res, y_res)

	open(unit=1, file='result.dat', action="write")
		do i=0, 100*n -1
			write(1, *) x_res(i), y_res(i)
		enddo
	close(1)

end program main