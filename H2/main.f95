program h2
	use fun
	implicit none
	!Требуется написать программу, вычисляющую значения интерполяционного
	!полинома для функции в произвольных промежуточных точках.
	integer   :: N, i
	real      :: a, b, t, L, arbitrary_x, pi=3.141592
	character :: value
	real, allocatable :: X(:),Y(:)

	call getarg(1,value)
	select case (value)
		case("u")
			open(1,file="uniform.dat") 
    		case("c")
			open(1,file="chebyshev.dat") 
	end select
   	read(1,'(2x,I6)') N
	read(1,*) a,b
	allocate(Y(N+1))
        allocate(X(N+1))
        read(1,*) Y
	close(1)
            select case(value)
 		case("u")
    		    do i=0,N
			print*, 'first step'
		        X(i)=(b-a)*(1.*i/N)+a
   		    enddo
		    open(1,file="res_uniform.dat")
		        do i=0,100*N
			print*, 'second step'
			    arbitrary_x=(b-a)*(1.*i/(100.*N))+a
			    call Lagrange(N,X,Y,arbitrary_x,L)
				print*, 'third step'
			    write(1,*) arbitrary_x, L
			    print*, 'fifth step'
		        enddo
		    close(1)
    	        case("c")
		    t=(abs(b-a))/2
   			do i=0,N
			    X(i)=a+t-t*cos(((2*i+1.)/(2.*N+2))*pi)
    		        enddo
			open(1,file="res_chebyshev.dat")
			do i=0,100*N
				arbitrary_x=a+t-t*cos(((2*i+1.)/(200.*N+2))*pi)
				call Lagrange(N,X,Y,arbitrary_x,L)
				write(1,*) arbitrary_x, L
			end do
			close(1)
	    end select
    deallocate(X,Y)
end program h2


