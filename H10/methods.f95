module methods
use ffunction
use newtonmethod
use lead
use mprecision
use adams
use rungekut
implicit none
contains

subroutine Rosenbroke()
	real(mp) alpha,beta,gammma, t
	real(mp), dimension(size(Xdata)) :: results,B,help
	real(mp), dimension(size(Xdata), size(Xdata)) :: A, J, E
	real(mp), dimension(size(Xdata), size(Xdata)+1) :: M
	real(mp), dimension(size(Xdata)+1, size(Xdata)) :: M_t
	integer i,n
		open(100, action='write', file='rosen.dat')
		n=size(Xdata)
		alpha=1.077
		beta=-0.372
		gammma=-0.577
		E=0
		results=Xdata
		forall(i=1:n)E(i,i)=1
		do while (t<tdata+h)
        	print*, results
		call yakobmatrix(results,J,f)
		A=(E-alpha*h*J-beta*h*h*matmul(J,J))/h
        help=gammma*h*f(0.0,results)+results
        B=matmul(A,results)+f(0.0,help)
        M = reshape([A, B], shape=(/n,n+1/))
        M_t = transpose(M)
	    call leadfun(M_t,results,n)
        print*, M(1, 2), 'aa', results
        t=t+h
        write(100,*)t,results
    enddo
    close(100)
  end subroutine
  
subroutine pred_cor_method()
            
    real(mp), dimension(size(Xdata)) :: X, Xnew
    real(mp) :: Xae(0:Nextradams, size(Xdata)), Xai(0:Nextradams-2, size(Xdata))
    real(mp), dimension(size(Xdata)) :: predict_correct, predictions
    real(mp) :: t_e, t_i, h_e, h_i, t
    integer :: i, j, n

    open(1, file='precor.dat', action='write')

    n = size(Xdata)
    t = 0.0_mp
    predict_correct = Xdata
    write(1,*) t, predict_correct

    h_e = h / Nextradams
    h_i = h / Nextradams 

    t_i = t  

    do while (t < tdata)  
        t_e = t_i
        
        Xae(0, :) = predict_correct
        do i=0, Nextradams-2
            call rk4(Xae(i,:), t_e, Xae(i+1,:), h_e)
            t_e = t_e + h_e
        enddo
        X = Xae(Nextradams-1,:)
        do while (t_e < t + h)
                call extradams(X, t_e, Xnew, Xae, h_e)
                forall (j=0:Nextradams-2) Xae(j,:)=Xae(j+1,:)
                Xae(Nextradams-1,:)=Xnew
                t_e = t_e + h_e
                X=Xnew
        enddo

        predictions = Xnew




        t_i = t_i + h_i
        call rk4(predict_correct, t_i, Xai(0,:), h_i)
        do i = 0, Nextradams-3
            t_i = t_i + h_i
            call rk4(Xai(i,:), t_i, Xai(i+1,:), h_i)
        enddo

        X = Xai(Nextradams-2, :)
        do while (t_i < t + h)
            call interadams(X, t_i, Xnew, Xai, h_i, predictions)
            forall (j=0:Nextradams-3) Xai(j,:) = Xai(j+1,:)
            Xai(Nextradams-2,:) = Xnew
            t_i = t_i + h_i
            X = Xnew
        enddo

        predict_correct = Xnew

        t = t + h
        write(1,*) t_i, predict_correct
        enddo
    close(1)
    end subroutine pred_cor_method

end module 


