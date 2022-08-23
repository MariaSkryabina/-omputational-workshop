results: main
	rm -f *.mod
#!  make rk - Метод Рунге-Кутты        
#!  make ae - Экстраполяционный метод Адамса  
#!  make ai - Интерполяционный метод Адамса   


main:  
	gfortran myprec.f90 lead.f90 newtonmethod.f90 ffunction.f90 ABfunc.f90  rungekut.f90 adams.f90 main.f90  -o main

rk: main
	./main rk

ae: main
	./main ae

ai: main
	./main ai
