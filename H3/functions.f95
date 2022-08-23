module funs
contains
    subroutine gauss(C, x)
        real(8), intent(inout), dimension(:, :) :: C
        real(8), intent(inout), dimension(:) :: x
        real(8), parameter :: eps = 1e-8
        integer :: i, j, k, n
        n = size(C, dim=1)

        do k=1,n
            if (abs(C(k,k)) < eps) then 
                write(*,*) 'WARNING!! Ведущий элемент близок к нулю!!'
            endif
            forall (j=k:n+1) 
                C(k,j) = C(k,j)/C(k,k) 
            end forall   
            forall (i=k+1:n, j=k:n+1) 
                C(i,j) = C(i,j) - C(k,j)*C(i,k)
            end forall
        enddo
        do i=n,1,-1
            x(i) = C(i, n+1) - dot_product(C(i, i+1:n), x(i+1:n))
        enddo
    end subroutine gauss

  subroutine jordan(C, x)
        real(8), intent(inout), dimension(:, :) :: C
        real(8), intent(inout), dimension(:) :: x
        real(8), parameter :: eps = 1e-8
        integer :: i, j, k, n

        n = size(C, dim=1)
        do k=1,n
            if (abs(C(k,k)) < eps) then 
	        write(*,*) 'WARNING!! Ведущий элемент близок к нулю!!'
            endif
            forall (j=k:n+1) 
                C(k,j) = C(k,j)/C(k,k) 
            end forall   
            forall (i=1:n, j=k:n+1, i/=k) 
                C(i,j) = C(i,j) - C(k,j)*C(i,k)
            end forall
        enddo
        do i=n,1,-1
            x(i) = C(i, n+1)
        enddo  
  end subroutine jordan

  subroutine choice(C, x)
      real(8), intent(inout), dimension(:, :) :: C
      real(8), intent(inout), dimension(:) :: x   
      real(8), dimension(size(C, dim=1)+1) :: temp_row
      real(8), dimension(size(C, dim=1)) :: temp_col
      real(8) :: temp
      integer :: i, j, k, n, indexx(2), rem(size(C, dim=1)), max_row_index, max_col_index
      n = size(C, dim=1)  

! rem- remembers transpositions
! temp_col temp_row нужны вследствие того, что мы меняем местами строки и столбцы. 
!rem запоминает все эти перестановки, чтобы потом восстановить x
! temp нужен для того, чтобы поменять местами 2 элемента в x

      do k=1,n
        	indexx = maxloc(abs(C(k:n, k:n)))
        	max_row_index = indexx(1)+k-1
        	max_col_index = indexx(2)+k-1
        	rem(k) = max_col_index 

        	temp_col = C(:, k)
        	C(:, k) = C(:, max_col_index) 
        	C(:, max_col_index) = temp_col        

        	temp_row = C(k,:)
        	C(k,:) = C(max_row_index,:) 
        	C(max_row_index,:) = temp_row

        	forall (j=k:n+1) 
            	C(k,j) = C(k,j)/C(k,k) 
        	end forall   

        	forall (i=k+1:n, j=k:n+1) 
           	    C(i,j) = C(i,j) - C(k,j)*C(i,k)
       	    end forall
      enddo

      do i=n,1,-1
          x(i) = C(i, n+1) - dot_product(C(i, i+1:n), x(i+1:n))
      enddo

      do i=n, 1, -1
        temp = x(i)
        x(i) = x(rem(i))
        x(rem(i)) = temp
      enddo


  end subroutine choice

end module
