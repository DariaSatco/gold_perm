program goldeps

use gold

implicit none

real, parameter:: pi=3.14159265359
real(8), parameter :: kb=8.617e-5  ! eV/K kb=1.38e-23 J/K  Boltzman constant
real(8), parameter :: h=6.585e-16  ! eV*s h=1.054e-34J*s  Plank constant
real(8), parameter :: eV=1.519e15 !rad/s

integer,parameter:: n1=663 !rows number in file
integer,parameter:: nw=663 !number of Matsubara frequencies

real(8):: matrixAu(n1,3)
real(8):: x(nw), eps_Drude(nw), eps_Mar(nw), eps_Gen(nw)
real(8):: epsAur(nw), eps_BB(nw), eps_LD(nw)
!counters
integer:: i,k,j
!integration results
real(8):: integralA1,integralA2
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funvalAu(n1),SpcoefA(3,n1), exintA(n1), workA(n1)
!variables for Kr-K integral in qags subroutine
integer(4)::nevalA,ierA


real(8):: T
real(8):: Lorentz_oscil_re(nw), Lorentz_oscil_im(nw)
real(8):: Mostep_oscil_re(nw), Mostep_oscil_im(nw)
real(8):: BB_oscil_re(nw), BB_oscil_im(nw)
real(8):: LD_model_re(nw), BB_mod_re(nw), Gen_pl_mod(nw)
real(8):: LD_model_im(nw), BB_mod_im(nw)
real(8):: epsMar_re(nw), epsMar_im(nw)
real(8):: Drude_re(nw), Drude_im(nw)
real(8):: force_LD(6), force_BB(6)

print*, "Process..."
open(unit=15, file='gold_eps_im_re_ev+Olmon.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

read(15,*)

do i=1,n1
read(15,*), (matrixAu(i,j),j=1,3)
end do

T=300._8
!Kramers-Kronig calculation for eps(iw)

	do i=1,nw
	
!	x(i)=2*pi*kb*T*i/h
!   x(i)=x(i)/eV
    x(i)=matrixAu(i,1)
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
	epsAur(i)=(integralA1+integralA2)*2/pi+1.0

	end do

open(unit=18, file='golddata.txt', status='replace')
open(unit=19, file='gold_oscillators_real.txt', status='replace')
open(unit=20, file='gold_oscillators_imag.txt', status='replace')
open(unit=21, file='gold_drude_plasma_re.txt', status='replace')
open(unit=22, file='gold_drude_im.txt', status='replace')
open(unit=25, file='eps(w)_re.txt', status='replace')
open(unit=26, file='eps(w)_im.txt', status='replace')

!calculate eps(iw) in different analytical models

!calculate oscillator forces from Rakic article
force_LD=fR*wpR**2
force_BB=fRb*wpR**2

do i=1,nw
eps_Drude(i)=Drude(x(i),parAproxIm)
eps_Mar(i)=epsMar(x(i))
eps_Gen(i)=Gen_Plasma(x(i),parMost, gAu,gammaAu,wAu)
eps_BB(i)=Brendel_Bormann(x(i), parRakic_BB, force_LD, gammaRb, wRb, sigmaRb)
eps_LD(i)=Lorentz_Drude(x(i), parRakic_LD, force_BB, gammaR, wR)
end do

!calculate real parts of oscillators in different models and parameters
do i=1,nw
Lorentz_oscil_re(i)=oscil_gold_re_im(x(i), force_LD, gammaR, wR,6, 0)
Mostep_oscil_re(i)=oscil_gold_re_im(x(i), gAu, gammaAu, wAu,6, 0)
BB_oscil_re(i)=Brendel_Bormann_oscil_Re_Im(x(i), force_BB, gammaRb, wRb, sigmaRb,6, 0)
end do

!calculate imaginary parts of oscillators in different models and parameters
do i=1,nw
Lorentz_oscil_im(i)=oscil_gold_re_im(x(i), force_LD, gammaR, wR, 6, 1)
Mostep_oscil_im(i)=oscil_gold_re_im(x(i), gAu, gammaAu, wAu, 6, 1)
BB_oscil_im(i)=Brendel_Bormann_oscil_Re_Im(x(i), force_BB, gammaRb, wRb, sigmaRb, 6, 1)
end do

!calculate real parts Drude/Plasma models
do i=1,nw
LD_model_re(i)=Drude_Re_Im(x(i), parRakic_LD, 0)
Gen_pl_mod(i)=Plasma(x(i), parMost, 0)
BB_mod_re(i)=Drude_Re_Im(x(i), parRakic_BB, 0)
end do

!calculate imaginary parts Drude models
do i=1,nw
LD_model_im(i)=Drude_Re_Im(x(i), parRakic_LD, 1)
BB_mod_im(i)=Drude_Re_Im(x(i), parRakic_BB, 1)
end do

100 format(7(f15.3))
150 format(7(A15))

!golddata.txt - eps(iw)
!create the top of the table
write(18,150) 'freq (eV)', 'K-K (Drude)', 'Drude', 'Marachevsky', 'Gen_plasma', 'Bren_Borm', 'Lor-Dr'
!write values in the table
do i=1,nw
write(18,100) x(i), epsAur(i), eps_Drude(i), eps_Mar(i), eps_Gen(i), eps_BB(i), eps_LD(i)
end do

200 format(6(f15.3))
250 format(6(A15))

!gold_oscillators_real.txt
!create the top of the table
write(19,250) 'freq (eV)', 'Lorentz', 'Gen_Plasma', 'Bren_Borm'
!write values in the table
do i=1,nw
write(19,200) x(i), Lorentz_oscil_re(i), Mostep_oscil_re(i), BB_oscil_re(i)
end do

!gold_oscillators_imag.txt
!create the top of the table
write(20,250) 'freq (eV)', 'Lorentz', 'Gen_Plasma', 'Bren_Borm'
!write values in the table
do i=1,nw
write(20,200) x(i), Lorentz_oscil_im(i), Mostep_oscil_im(i), BB_oscil_im(i)
end do

!gold_drude_plasma_re.txt
!create the top of the table
write(21,250) 'freq (eV)', 'Lorentz', 'Gen_Plasma', 'Bren_Borm'
!write values in the table
do i=1,nw
write(21,200) x(i), LD_model_re(i), Gen_pl_mod(i), BB_mod_re(i)
end do

!gold_drude_im.txt
!create the top of the table
write(22,250) 'freq (eV)', 'Lorentz', 'Bren_Borm'
!write values in the table
do i=1,nw
write(22,200) x(i), LD_model_im(i), BB_mod_im(i)
end do

!eps(w)_re.txt
write(25,250) 'freq (eV)', 'Lorentz', 'Gen_Plasma', 'Bren_Borm', 'Marachevsky', 'Drude'

do i=1,nw

!Marachevsky
epsMar_re(i)=epsMar_re_im(x(i), 0)
Drude_re(i)=Drude_Re_Im(x(i), parAproxIm, 0)

write(25,200) x(i), LD_model_re(i)+Lorentz_oscil_re(i), Gen_pl_mod(i)+Mostep_oscil_re(i), BB_mod_re(i)+BB_oscil_re(i), &
epsMar_re(i), Drude_re(i)
end do

!eps(w)_im.txt
write(26,250) 'freq (eV)', 'Lorentz', 'Gen_Plasma', 'Bren_Borm', 'Marachevsky', 'Drude'

do i=1,nw

!Marachevsky
epsMar_im(i)=epsMar_re_im(x(i), 1)
Drude_im(i)=Drude_Re_Im(x(i), parAproxIm, 1)

write(26,200) x(i), LD_model_im(i)+Lorentz_oscil_im(i), Mostep_oscil_im(i), BB_mod_im(i)+BB_oscil_im(i),&
epsMar_im(i), Drude_im(i)
end do

!experimental data - only oscillation part
open(unit=23, file='eps_minus_drude_re_im.txt', status='replace')

!eps(w)_im.txt
write(23,250) 'freq (eV)', 're_eps', 'im_eps'

do i=1,nw
write(23,200) x(i), matrixAu(i,2)-Drude_re(i), matrixAu(i,3)-Drude_im(i)
end do


close(15)
close(18)
close(19)
close(20)
close(21)
close(22)
close(23)
close(25)
close(26)

print*, "Done!"

end program
