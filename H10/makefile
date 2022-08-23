results: main
	rm -f *.mod

main: myprec.f95 lead.f95 newtonmethod.f95 ffunction.f95 ABfunc.f95  rungekut.f95 adams.f95  methods.f95 main.f95  
	gfortran myprec.f95 lead.f95 newtonmethod.f95 ffunction.f95 ABfunc.f95  rungekut.f95 adams.f95  methods.f95 main.f95 -o main

