module fun
    implicit none
    contains
    subroutine Lagrange(N,X,Y,arbitrary_x,output)
        real(4), intent(out) :: output
        real(4), intent(in) :: X(:), Y(:)
	real, intent(in) :: arbitrary_x
	real :: L
        integer :: n,i,j
	
	output=0.00
	do i=1,N+1
		L=1.00
		do j=1,N+1
			if (i==j) then
				L=L*1.00
			else 
				L=L*(arbitrary_x-X(j))/(X(i)-X(j))
			endif
		enddo
	output=output+L*Y(i)
	enddo
	end subroutine Lagrange
end module fun
	
	
