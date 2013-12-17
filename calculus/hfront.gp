Kw=1e-14;
dd(pH)=10**(-pH)-Kw/(10**(-pH));


U(pH,pH_0,pH_inf) = (dd(pH)-dd(pH_inf))/(dd(pH_0)-dd(pH_inf));
ModU(pH,pH_0,pH_inf) = inverf(1-U(pH,pH_0,pH_inf));

V(pH,pH_0,pH_inf) = (10**(-pH)-10**(-pH_inf))/(10**(-pH_0)-10**(-pH_inf));
ModV(pH,pH_0,pH_inf) = inverf(1-V(pH,pH_0,pH_inf));

FreeOverChem(pH,pH_0,pH_inf) = ModV(pH,pH_0,pH_inf)/ModU(pH,pH_0,pH_inf);

Kernel(U,pH_0,pH_inf) = 10**(-pH_inf) + (10**(-pH_0)-10**(-pH_inf)) * erfc(U);

#plot [2.1:9.9] Coeff(10.0**(-x)) w l



