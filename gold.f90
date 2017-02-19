module gold

interface
    SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
    real(8):: xi, yi, u,v
    logical::flag
    end subroutine
end interface

real(8), public:: freqA
!parameters (Drude model)
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter:: parMost=(/9.0, 0.035/) !eV
!Rakic
real(8), dimension(2), parameter:: parRakic=(/7.87, 0.053/) !eV
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter:: parAproxIm=(/6.62, 0.089/) !eV
!approximation data of the real part (Palik)
real(8), dimension(2), parameter:: parAproxRe=(/7.93, 0.08/) !eV

!oscillator model, Mostepanenko parameters
real(8), dimension(6), parameter:: gAu=(/7.091, 41.46, 2.7, 154.7, 44.55, 309.6/) !eV
real(8), dimension(6), parameter:: gammaAu=(/0.75, 1.85, 1.0, 7.0, 6.0, 9.0/) !eV
real(8), dimension(6), parameter:: wAu=(/3.05, 4.15, 5.4, 8.5, 13.5, 21.5/) !eV

!Lorentz-Drude model, Rakic parameters
real(8), dimension(6), parameter:: fR=(/0.024, 0.01, 0.071, 0.601, 4.384, 0.0/) !eV
! gR=fR*parRakic(1)**2
real(8), dimension(6), parameter:: gammaR=(/0.241, 0.345, 0.87, 2.494, 2.214, 0.0/) !eV
real(8), dimension(6), parameter:: wR=(/0.415, 0.83, 2.969, 4.304, 13.32, 0.0/) !eV

!Brendel-Bormann model, Rakic parameters
real(8), dimension(6), parameter:: fRb=(/0.054, 0.05, 0.312, 0.719, 1.648, 0.0/)
! gRb=fRb*parRakic(1)**2
real(8), dimension(6), parameter:: gammaRb=(/0.074, 0.035, 0.083, 0.125, 0.179, 0.0/)
real(8), dimension(6), parameter:: wRb=(/0.218, 2.885, 4.069, 6.137, 27.97, 0.0/)
real(8), dimension(6), parameter:: sigmaRb=(/0.742, 0.349, 0.830, 1.246, 1.795, 1.0/)


contains

function Drude(x, param)
!Drude-model permettivity function eps(iw)
!param - vector of parameters (wp, damp)

real(8):: x, Drude, param(2)

	Drude=1+param(1)**2/(x**2+param(2)*x);

end function

!----------------------------------------------------

function epsMar(x)
!Marachevsky formula fo eps(iw)
real(8):: x,epsMar,num,denum, eV
real(8), dimension(2,4):: paramMar

eV=1.519e15
!Marachevsky model parameters (wl1,wl2//gl1,gl2//gt1,gt2//wt2,0)
paramMar(1,1)=exp(-0.96)*1e16
paramMar(2,1)=exp(0.2866)*1e16
paramMar(1,2)=exp(-2.536)*1e16
paramMar(2,2)=exp(1.255)*1e16
paramMar(1,3)=exp(-4.7922)*1e16
paramMar(2,3)=exp(-0.957)*1e16
paramMar(1,4)=exp(-0.8359)*1e16
paramMar(2,4)=0.0_8

paramMar=paramMar/eV

num=(paramMar(1,1)**2+x**2+paramMar(1,2)*x)*(paramMar(2,1)**2+x**2+paramMar(2,2)*x)
denum=(x**2+paramMar(1,3)*x)*(paramMar(1,4)**2+x**2+paramMar(2,3)*x)

epsMar=num/denum

end function

!---------------------------------------------------

function oscil_gold_im(x,a,b,c)
!imaginary part of oscillator approximation (Mostepanenko/Racik/other...)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: oscil_gold_im, x, a(6), b(6), c(6)

oscil_gold_im=0.0_8

do i=1,size(a)
oscil_gold_im=oscil_gold_im-a(i)*b(i)*x/((c(i)**2-x**2)**2+(b(i)*x)**2)
end do

end function
!--------------------------------------------------

function oscil_gold_re(x,a,b,c)
!real part of oscillator approximation (Mostepanenko/Racik/other...)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: oscil_gold_re, x, a(6), b(6), c(6)

oscil_gold_re=0.0_8

do i=1,size(a)
oscil_gold_re=oscil_gold_re+a(i)*(c(i)**2-x**2)/((c(i)**2-x**2)**2+(b(i)*x)**2)
end do

end function
!--------------------------------------------------

function oscil_gold(x,a,b,c)
!Lorentz function (oscillators) eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: oscil_gold, x, a(6), b(6), c(6)

oscil_gold=0.0_8

do i=1,size(a)
oscil_gold=oscil_gold+a(i)/(c(i)**2+x**2-b(i)*x)
end do

end function

!--------------------------------------------------

!function oscil_to_int(x)
!!KK expression under integral for oscillators
!
!real(8):: x, oscil_to_int
!
!oscil_to_int=x/(x**2+freqA**2)*oscil_gold(x)
!
!end function

!--------------------------------------------------

function drude_to_int(x)
!KK expression under integral for Drude model

real(8):: x, drude_to_int, param(2)

param=parAproxIm
drude_to_int=param(1)**2*param(2)/((x**2+param(2)**2)*(x**2+freqA**2))

end function


!-------------------------------------------------

function Gen_Plasma(x,a,b,c)
!Generalized plasma-like dielectric permettivity function eps(iw)

real(8):: x, Gen_Plasma
real(8):: a(6), b(6), c(6)

Gen_Plasma=1._8+(parMost(1)/x)**2+oscil_gold(x,a,-b,c)

end function

!-------------------------------------------------

function Lorentz_Drude(x,a,b,c)
!Lorentz-Drude model (Rakic) eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: x, Lorentz_Drude, a(6), b(6), c(6)

Lorentz_Drude=Drude(x, parRakic)+oscil_gold(x,a,b,c)

end function

!-------------------------------------------------

function Brendel_Bormann_oscil(x,a,b,c,d)
!Brendel-Bormann model, oscillation part in eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!d - vector of sigma (dispersion)

real, parameter:: pi=3.14159265359

real(8):: x, Brendel_Bormann_oscil, a(6),b(6),c(6),d(6)
real(8):: WvalueRe1, WvalueIm1, WvalueRe2, WvalueIm2
logical:: logvar
real(8):: Re,Im, const
!argW=Re+i*Im

Brendel_Bormann_oscil=0.0D+00

do i=1,size(a)

!calculate argument
Re=c(i)/(sqrt(2.0)*d(i))
Im=sqrt(x**2+x*b(i))/(sqrt(2.0)*d(i))
!call Faddeeva function w(z)
call wofz(Re, Im, WvalueRe1, WvalueIm1, logvar) !first summand
call wofz(-Re, Im, WvalueRe2, WvalueIm2, logvar) !second summand


!calculate prefactor
const=sqrt(pi)*a(i)/(2*sqrt(2.0)*d(i)*sqrt(x**2+x*b(i)))

if (logvar.eqv. .True.) then
print*, 'Attention, error in wofz subroutine'
else
Brendel_Bormann_oscil=Brendel_Bormann_oscil+(WvalueRe1+WvalueRe2)*const
end if

end do

end function

!-------------------------------------------------
function Brendel_Bormann(x,param,a,b,c,d)
!Brendel-Bormann model, eps(iw)

real(8):: Brendel_Bormann, x
real(8):: a(6),b(6),c(6),d(6), param(2)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!d - vector of sigma (dispersion)

Brendel_Bormann=Drude(x, param)+Brendel_Bormann_oscil(x, a, b, c, d)

end function

!-------------------------------------------------

function epsRe(x)
!real part of dielectric perm. function, Drude model

real(8):: x, epsRe
real(8):: param(2)

param=parAproxRe

epsRe=1._8-param(1)**2/(x**2+param(2)**2)

end function


!--------------------------------------------------


end module



