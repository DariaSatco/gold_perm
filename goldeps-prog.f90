program goldeps

use gold

implicit none

real, parameter:: pi=3.14159265359
real(8), parameter :: kb=8.617e-5  ! eV/K kb=1.38e-23 J/K  Boltzman constant
real(8), parameter :: h=6.585e-16  ! eV*s h=1.054e-34J*s  Plank constant

integer,parameter:: n1=310 !rows number in file
integer,parameter:: nw=500 !number of Matsubara frequencies

real(8):: matrixAu(n1,3) , x(nw), eps_Drude(nw), eps_Mar(nw), eps_Gen(nw), epsAur(nw), eps_BB(nw)
!counters
integer:: i,k,j
!integration results
real(8):: integralA1,integralA2,integralA3
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funvalAu(n1),SpcoefA(3,n1), exintA(n1), workA(n1), funvalAu1(n1)
!variables for Kr-K integral in qags subroutine
integer(4)::nevalA,ierA
real(8):: eV,T


print*, "Process..."
open(unit=15, file='resultAu.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

open(unit=16, file='resultAu_eV.txt', status='replace')

eV=1.519e15 !rad/s
200 format(3(f15.3))

do i=1,n1
read(15,*), (matrixAu(n1+1-i,j),j=1,3)
end do

do i=1,n1
matrixAu(i,1)=matrixAu(i,1)/eV
write(16,200), (matrixAu(i,j),j=1,3)
end do

T=300._8
!Kramers-Kronig for 'freq'-vector

	do i=1,nw
	
!	x(i)=2*pi*kb*T*i/h
!   x(i)=x(i)/eV
    x(i)=0.2_8+0.01*(i-1)
    freqA=x(i)

	!integrates the Drude-model approximation
	call qnc79(drude_to_int, 0.0_8, matrixAu(1,1), 1.0e-3_8, integralA1, ierA, nevalA)

	!builds vector to integrate from the Palik data
		do k=1,n1
		funvalAu(k)=matrixAu(k,3)*matrixAu(k,1)/(matrixAu(k,1)**2+freqA**2)
		end do
	!integrates the vector of data
 	call cspint(n1, matrixAu(1:n1,1), funvalAu, matrixAu(1,1), matrixAu(n1,1), SpcoefA, exintA, workA, integralA2)

	!Kr-Kr formula for permittivity 
	epsAur(i)=(integralA1+integralA2)*2/pi+1

	end do

open(unit=18, file='golddata.txt', status='replace')

do i=1,nw
eps_Drude(i)=Drude(x(i),parAproxIm)
eps_Mar(i)=epsMar(x(i))
eps_Gen(i)=Gen_Plasma(x(i),gAu,gammaAu,wAu)
eps_BB(i)=Brendel_Bormann(x(i), parRakic, fRb*parRakic(1)**2, gammaRb, wRb, sigmaRb)
end do

100 format(7(f15.3))
150 format(6(A15))

!create the top of the table
write(18,150) 'freq (eV)', 'K-K (Drude)', 'Drude', 'Marachevsky', 'Gen_plasma', 'Bren_Borm'
!write values in the table
do i=1,nw
write(18,100) x(i), epsAur(i), eps_Drude(i), eps_Mar(i), eps_Gen(i), eps_BB(i)
end do

close(15)
close(16)
close(18)

print*, "Done!"

end program
