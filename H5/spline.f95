module spline_module
    use my_prec ! mp
    use run_module ! run_solver
    implicit none
    contains
    subroutine f(x, y, weigth, x_res, y_res)
        real(mp), intent(inout), dimension(0:) :: x, y, weigth, x_res, y_res
        integer :: n, n_res
        real(mp), dimension(size(x), 3) :: A, B, BQ6, B_t, QB_t
        real(mp), dimension(size(x), 5) :: product_all, A_pad 
        real(mp), dimension(0:size(x)-1) :: S, R, BY6
        real(mp) :: h, t, start, finish, x_i, h_res, tmp
        integer :: i, index
        n = size(x) - 1

        ! compute A, B
        A = 0
        B = 0
        A(1, 1) = 2*(x(1) - x(0))
        A(2, 2) = 2*(x(2) - x(0))
        A(2, 3) = x(2) - x(1)
        A(n, 1) = x(n-1) - x(n-2)
        A(n, 2) = 2*(x(n) - x(n-2))
        A(n+1, 3) = 2*(x(n) - x(n-1))

        do i=2, n
            B(i, 1) = 1 / (x(i-1) - x(i-2))
            B(i, 2) = -1 / (x(i-1) - x(i-2)) - 1/(x(i) - x(i-1))
            B(i, 3) = 1 / (x(i) - x(i-1))

            if (i == 2 .or. i == n) then 
                continue
            else 
                A(i, 1) = x(i-1) - x(i-2)
                A(i, 2) = 2*(x(i) - x(i-2))
                A(i, 3) = x(i) - x(i-1)
            endif
        enddo

        ! compute BQ
        ! когда мы умножаем справа на диагональную матрицу, мы должны каждый столбец B домножить на число из Q
        ! с соотвествующим индексом (матрицы Q нет, вместо нее просто делим на weigth(индексация weigth с 0))
        BQ6 = B
        do i=2, n
            BQ6(i, 1) = B(i, 1) / weigth(i-2)
            BQ6(i, 2) = B(i, 2) / weigth(i-1)
            BQ6(i, 3) = B(i, 3) / weigth(i)
        enddo
        BQ6 = 6*BQ6  

        !compute B_T
        ! для траспонирования B нужно поменять местами только значения (1,2)<->(2,1) и (n,3)<->(n+1,2)
        B_t = B
        B_t(1, 2) = B_t(2, 1)
        B_t(2, 1) = 0
        B_t(n+1, 2) = B_t(n, 3)
        B_t(n, 3) = 0        
        B_t(n+1, 1:2) = B_t(n+1, 2:)     
 
        ! compute BY6
        call triagonal_mul_v(B, y, BY6)
        BY6 = 6*BY6

        !compute product_all
        call triagonal_mul(BQ6, B_t, product_all)

        ! add padding in A, 
        A_pad = 0
        A_pad(1, 1:3) = A(1, :) 
        A_pad(2, 1:3) = A(2, :)
        A_pad(n+1, 3:5) = A(n+1, :)
        A_pad(n, 3:5) = A(n, :)
        
        do i=3, n-1
            A_pad(i, 2:4) = A(i, :) 
        end do
        
        product_all = A_pad + product_all

        ! compute S
        call run_solver(product_all, BY6, S)


        QB_t = B
        do i=2, n
            QB_t(i, 1) = QB_t(i, 1) / weigth(i-2)
            QB_t(i, 2) = QB_t(i, 2) / weigth(i-1)
            QB_t(i, 3) = QB_t(i, 3) / weigth(i)
        enddo
        tmp = QB_t(1, 2)
        QB_t(1, 2) = QB_t(2, 1)
        QB_t(2, 1) = tmp
        tmp = QB_t(n+1, 2)
        QB_t(n+1, 2) = QB_t(n, 3)
        QB_t(n, 3) = tmp
        do i=2, n-1
            tmp = QB_t(i, 3)
            QB_t(i, 3) = QB_t(i+1, 1)
            QB_t(i+1, 1) = tmp
        enddo

        call triagonal_mul_v(QB_t, S, R)
        R = y - R
!         print*, R

        ! make x_res and f(x_res) =: y_res
        start = x(0)
        finish = x(n)
        n_res = 100*n 
        h_res = (finish - start) / n_res
        x_res(0) = x(0)

        do i=1, n_res
            x_res(i) = x_res(i-1) + h_res
        enddo
        y_res(0) = R(0)

        ! compute y_res
        do i=0, n_res-1
            ! ищем индекс начала промежутка
            index =  sum(maxloc(x, mask=x_res(i)>=x))  - 1
            h = x(index + 1) - x(index)
            x_i = x(index)
            t = (x_res(i) - x_i)/h
            y_res(i) = R(index)*(1-t) + R(index+1)*t - h**2*t*(1-t)/6*((2-t)*S(index) + (1+t)*S(index+1))
        enddo

        return 
    end subroutine f


    subroutine triagonal_mul_v(A, x, y)
            !input  :: A - 3-diagonal matrix (n x 3)
            !       :: x - vector (n x 1)
            !return :: y - result of multiplication of A and y(n x 1)
        real(mp), intent(in), dimension(:,:) :: A
        real(mp), intent(in), dimension(0:) :: x
        real(mp), intent(out), dimension(0:) :: y
        integer :: i, n
        n = size(x) - 1
        ! compute y
        y(0) = A(1, 1)*x(0) + A(1, 2)*x(1)
        y(n) = A(n+1, 1)*x(n-2) + A(n+1, 2)*x(n-1) + A(n+1,1)*x(n)
        do i=1, n-1
            y(i) = A(i+1, 1)*x(i-1) + A(i+1, 2)*x(i) + A(i+1, 3)*x(i+1)
        enddo

        return
        end subroutine triagonal_mul_v


        subroutine triagonal_mul(A, B, X)
            !input  :: A, B - 3-diagonal matrix (n x 3)
            !return :: X - result of multiplication of A by B(n x 5)
        implicit none
        real(mp), intent(in), dimension(:,:) :: A, B
        real(mp), intent(out), dimension (:,:) :: X
   
        integer :: i, j, n
        n = size(A, dim=1)

        !initilize with zero
        X = 0

        X(1, 1) = A(1,1)*B(1,1) + A(1,2)*B(2,1) 
        X(2, 1) = A(2,1)*B(1,1) + A(2,2)*B(2,1)


        X(1, 2) = A(1,1)*B(1,2) + A(1,2)*B(2,2)
        X(2, 2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,1) 
        do i=3,n
            X(i, 1) = A(i,1)*B(i-1,1)
            X(i, 2) = A(i,1)*B(i-1,2) + A(i,2)*B(i,1)
            X(i, 3) = A(i,1)*B(i-1,3) + A(i,2)*B(i,2) + A(i,3)*B(i+1,1)
        enddo

        X(1, 3) = A(1,2)*B(2,3)
        X(2, 3) = A(2,2)*B(2,3) + A(2,3)*B(3,2)


        X(2, 4) = A(2,3)*B(3,3)
        X(n-1,4) = A(n-1,2)*B(n-1,3) + A(n-1,3)*B(n,2)
        do i=3,n-2
            X(i, 4) = A(i,2)*B(i,3) + A(i,3)*B(i+1,2)
            X(i,5) = A(i,3)*B(i+1,3)
        enddo

        X(n-1, 2:) = X(n-1, :4)
        X(n-1, :1) = 0

        X(n, 3:) = X(n, :3)
        X(n, :2) = 0

        return
        end subroutine triagonal_mul
end module spline_module
