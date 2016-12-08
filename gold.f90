module gold

real(8), public:: freqA
!parameters (Drude model)
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter:: parMost=(/9.0, 0.035/) !eV
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter:: parAproxIm=(/6.62, 0.089/) !eV
!approximation data of the real part (Palik)
real(8), dimension(2), parameter:: parAproxRe=(/7.93, 0.08/) !eV

!oscillator model, Mostepanenko parameters
real(8), dimension(6), parameter:: gAu=(/7.091, 41.46, 2.7, 154.7, 44.55, 309.6/) !*(1.519e15)**2
real(8), dimension(6), parameter:: gammaAu=(/0.75, 1.85, 1.0, 7.0, 6.0, 9.0/) !*1.519e15
real(8), dimension(6), parameter:: wAu=(/3.05, 4.15, 5.4, 8.5, 13.5, 21.5/) !*1.519e15

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

function oscil_gold(x)
!imaginary part of oscillator approximation

real(8):: oscil_gold, x

oscil_gold=0.0_8

do i=1,6
oscil_gold=oscil_gold+gAu(i)*gammaAu(i)*x/((wAu(i)**2-x**2)**2+(gammaAu(i)*x)**2)
end do

end function

!--------------------------------------------------

function oscil_to_int(x)
!KK expression under integral fo oscillators

real(8):: x, oscil_to_int

oscil_to_int=x/(x**2+freqA**2)*oscil_gold(x)

end function

!--------------------------------------------------

function drude_to_int(x)
!KK expression under integral for Drude model

real(8):: x, drude_to_int, param(2)

param=parAproxRe
drude_to_int=param(1)**2*param(2)/((x**2+param(2)**2)*(x**2+freqA**2))

end function


!-------------------------------------------------

!Generalized plasma-like dielectric permettivity function

function Gen_Plasma(x)

real(8):: x, Gen_Plasma

Gen_Plasma=1._8+(parMost(1)/x)**2

do i=1,6
Gen_Plasma=Gen_Plasma+gAu(i)/(wAu(i)**2+x**2-gammaAu(i)*x)
end do

end function
!-------------------------------------------------


end module



