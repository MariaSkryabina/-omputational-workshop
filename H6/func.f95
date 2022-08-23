module func

	contains

	function Ffun(X) result(Ff)
	! Функция Ffun(X) получает вектор X, и возвращает вектор Ff той же размерности
	implicit none
	real, dimension(1:), intent(in) :: X
	real, dimension(1:size(X)) :: Ff
	integer :: i, n

	n=size(X)
	
	Ff=0
	Ff(1)=3*X(1)+5.0
	do i=2,n
	    Ff(i)=product(X(i-1:i))+X(i)
	! Функция придумана для примера
	    Ff(i-1)=Ff(i-1)*3*X(n)-sum(X(1:i))  
	enddo

	end function Ffun

end module func
