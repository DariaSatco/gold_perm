program goldeps

use silicon
use gold
use dopcasimir
use spline
use silicon
use quadpack

implicit none

integer,parameter:: n1=310 !rows number in file
!Palik data table
real(8):: matrixAu(n1,3) , x(100), y1(100), y2(100)
!permettivity vector
real(8), allocatable:: epsAur(:)
!counters
integer:: i,k,j
!integration results
real(8):: integralA1,integralA2,integralA3
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funvalAu(n1),SpcoefA(3,n1), exintA(n1), workA(n1)
!variables for Kr-K integral in qags subroutine
real(4):: abserrA,abserrA2
integer(4)::nevalA,ierA,nevalA2,ierA2


open(unit=15, file='resultAu.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,n1
read(15,*), (matrixAu(n1+1-i,j),j=1,3)
end do

200 format(' ', e15.3, f15.7, f15.7)

!do i=1,n1
!write(*,200), (matrixAu(i,j),j=1,3)
!end do

!Kramers-Kronig for 'freq'-vector

allocate(epsAur(100))
freqA=20._8

	do i=1,100
	
	!integrates the Drude-model approximation
	call qnc79(drude_to_int, 0.0_8, matrixAu(1,1), 1.0e-3_8, integralA1, ierA, nevalA)
	PRINT*, 'INTEGRAL1=', integralA1
	!print*, 'ier=', ierA

	
	!builds vector to integrate from the Palik data
		do k=1,n1
		funvalAu(k)=matrixAu(k,3)*matrixAu(k,1)/(matrixAu(k,1)**2+freqA**2)
		end do
	!integrates the vector of data
 	call cspint(n1, matrixAu(1:n1,1), funvalAu, matrixAu(1,1), matrixAu(n1,1), SpcoefA, exintA, workA, integralA2)
	print*, 'integral2=', integralA2

	!integrates the oscillator-function
	call monte_carlo(oscil_to_int, matrixAu(n1,1), matrixAu(n1,1)*1.0e4_8, 1000000, integralA3)
	print*, 'integral3=', integralA3
	!print*, 'ier3=', ierA2

	
	!Kr-Kr formula for permittivity 
	epsAur(i)=(integralA1+integralA2+integralA3)*2/pi+1

	x(i)=freqA
	freqA=freqA*10**(0.2_8)

	end do

open(unit=18, file='golddata1.txt', status='replace')

do i=1,100
y1(i)=Drude(x(i),parAproxIm)
y2(i)=epsMar(x(i))
end do

do i=1,100
write(18,*) x(i), epsAur(i), y1(i), y2(i)
end do

deallocate(epsAur)
close(15)
close(18)

end program
