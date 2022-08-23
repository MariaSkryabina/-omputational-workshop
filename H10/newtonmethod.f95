	module newtonmethod
	use mprecision
	use lead

	contains

	subroutine newton(X0,nummax,X,f)
	! Процедура контролирует получение вектора решения на каком-то шаге итераций
	implicit none
	real(mp), dimension(1:), intent(in) :: X0 ! X0 - начальное приближение
	real(mp), dimension(1:size(X0)), intent(out) :: X
	real(mp), dimension(1:size(X0)) :: Xnew ! Xnew - вектор X, получаемый при новом шаге итераций
	real(mp) :: eps=0.1_mp**4
	integer, intent(in) :: nummax
	integer :: i, n
	interface
	function F(t, X)
	use mprecision
	real(mp), intent(in) :: t
	real(mp), dimension(1:), intent(in) :: X
	real(mp), dimension(1:size(X)) :: F
	end function F
	end interface

	n=size(X0)

	X=X0; i=nummax; Xnew=X0+1.0_mp ! Число выбрано случайно для устранения возможного совпадения с X0
	do while (sum(abs(X-Xnew))>eps .and. i>=1)
		X=Xnew
		call solve(X,Xnew,F)
		i=i-1
	enddo
	X=Xnew ! Здесь X присваивает значение найденного решения


	end subroutine newton

	subroutine solve(X,Xnew,F)
	! Процедура получает новый вектор решений по методу Ньютона (с помощью итераций)
	implicit none
	real(mp), dimension(1:), intent(in) :: X
	real(mp), dimension(1:size(X)), intent(out) :: Xnew
	real(mp), dimension(1:size(X),1:size(X)) :: df ! df - матрица Якоби функции f
	real(mp), dimension(1:size(X)+1,1:size(X)) :: M ! M - расширенная матрица системы f+sum(df*(Xnewk-Xk))=0
	integer :: i, n
	interface
	function F(t, X)
	use mprecision
	real(mp), intent(in) :: t
	real(mp), dimension(1:), intent(in) :: X
	real(mp), dimension(1:size(X)) :: F
	end function F
	end interface

	n=size(X)

	call yakobmatrix(X,df,F)
	M(n+1,1:n)=-F(0.0, X)
	M(1:n,1:n)=df(1:n,1:n)
	call leadfun(M,Xnew,n)
	Xnew=Xnew+X ! (Так как при решении системы, был получен вектор Xnew-X)


	end subroutine solve

	subroutine yakobmatrix(X0,df,F)
	! Процедура создает матрицу Якоби для f(x) в точке X0
        real(mp) :: x0(:), h
        real(mp) :: df(size(x0),size(x0))
        real(mp), dimension(size(x0)) :: fx0, fx1, fx2, x1, x2
        integer :: i, l, n
        interface
		function F(t, X)
		use mprecision
		real(mp), intent(in) :: t
		real(mp), dimension(1:), intent(in) :: X
		real(mp), dimension(1:size(X)) :: F
		end function F
        end interface
        n = size(X0)
        h = 0.000001_mp
        print*, x0, 'ssss'
        do l=1,n
            x1 = x0
            x2 = x0
            x1(l) = x0(l)+h
            x2(l) = x0(l)+2*h
            fx0 = f(0.0, x0)
            fx1 = f(0.0, x1)
            fx2 = f(0.0, x2)
            df(:,l) = -1.5_mp*fx0(:)/h + 2.0_mp*fx1(:)/h - 0.5_mp*fx2(:)/h
        enddo
	end subroutine yakobmatrix

	end module newtonmethod
