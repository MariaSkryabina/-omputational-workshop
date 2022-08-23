program main
	use funs
	implicit none
	real(8), allocatable :: A(:,:), B(:), x(:), C(:, :)
	integer :: i, n
	character(len=2) :: value !Ключ, который определяет способ решения СЛАУ
	call getarg(1, value) 

	open(1, file='data.dat')             
	    read(1,'(2X, I8)') n
	    allocate(A(n,n), B(n), x(n))
	    read(1,*) A; A = transpose(A)
	    read(1,*) B 
	close(1)

	C = reshape([A, B], shape=(/n,n+1/))
	if (value == '-G') then
		call gauss(C, x)
	elseif (value == '-J') then
		call jordan(C, x)
	elseif (value == '-L') then
		call choice(C, x)
	else
		write(*, *) 'Choose key : -G  Gauss way, -J  Jordan way, -L Lead way'
	endif
	open(1, file='result.dat')
		write(1, *) '# ', n
		do i=1,n
			write(1,*) x(i)
		enddo
	close(1)

	deallocate(A, B, x)

end program main
