module lead

	contains

	subroutine leadfun(M,X,n)
	! Процедура строит матрицу для решения методом Гаусса, 
	!  но на каждом шаге ставит наибольший по модулю элемент ведущим, 
	!  чтобы избегать деления на близкий к нулю элемент
	implicit none
	character(6) :: choice='lead'
	real, dimension(1:n+1,n) :: M
	real, dimension(1:n) :: X
	integer :: i, j, k, n
	real, dimension(n+1) :: Help 
	integer, dimension(2,n) :: Trans  

	! Тот же цикл, что и в методе Гаусса, но с перестановкой строк и столбцов
	do k=1,n-1 
	    if ( abs( maxval(M(k:n,k:n)) ) > abs( minval(M(k:n,k:n)) ) ) then 
		Trans(1:2,k)=(k-1)+maxloc(M(k:n,k:n))                        
	    else
		Trans(1:2,k)=(k-1)+minloc(M(k:n,k:n))                         
	    endif

	    Help(1:n+1)=M(1:n+1,Trans(2,k)) 
	    M(1:n+1,Trans(2,k))=M(1:n+1,k) 
	    M(1:n+1,k)=Help(1:n+1)         
	    Help(1:n)=M(Trans(1,k),1:n)   
	    M(Trans(1,k),1:n)=M(k,1:n)
	    M(k,1:n)=Help(1:n)

	! Оба "forall" строят матрицу с переставленными элементами "по Гауссу"
	    forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k) 
	    forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k)
	enddo
	M(n:n+1,n)=M(n:n+1,n)/M(n,n) 
	
	call solution(choice,M,X,n)
	
	! После выдачи процедурой solution массива из иксов, нужно его привести в порядок, потому что мы делали перестановки
	do k=n-1,1,-1
	    Help(k)=X(Trans(1,k))
	    X(Trans(1,k))=X(k)
	    X(k)=Help(k)
	enddo

	end subroutine leadfun

	subroutine solution(choice,M,X,n)
	! Процедура высчитывает массив решений, используя уже преобразованную расширенную матрицу
	implicit none
	integer :: i, n
	real, dimension(1:n+1,n) :: M
	real, dimension(1:n) :: X
	character(6) :: choice 

	! Для метода Жордана счет иксов ведется по-другому
	if (choice /= 'jordan') then 
	    do i=n,1,-1
		X(i)=M(n+1,i)-dot_product(M(i+1:n,i),X(i+1:n))
	    enddo
	else
	    do i=1,n
		X(i)=M(n+1,i)
	    enddo
	endif

	end subroutine solution

end module lead
