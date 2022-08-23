results: main
	./main

main: ffunction.f90 integral.f90
	gfortran -fbounds-check $^  -o main
	rm -f *.mod

tablichki: bernul.f90 lead.f90 problems.f90 quad.f90
	gfortran -fbounds-check $^  -o tablichki
	rm -f *.mod
