program newtonslau
        use newtonmethod
        use func

        implicit none
        real, allocatable :: X0(:), X(:) 
        integer :: n, i, nummax=1000

        open(unit=1, file='data.dat', action="read")
                read(1,'(2x,i8)') n
                allocate(X0(n),X(n))
                read(1,*) X0
        close(1)

        ! Получаем решение системы, используя начальное приближение X
        call newton(X0,nummax,X,Ffun) 
	 write(*,*) '|f(X)|=', sqrt(dot_product(Ffun(X),Ffun(X))) 
        open(unit=1, file='result.dat', action="write")
		do i=1,n
                	write(1,*) X(i)
		end do
        deallocate(X,X0)
        close(1)

end program newtonslau
