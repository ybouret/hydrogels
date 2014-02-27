1;
M0=87500;
Merr=12500;
Mup=M0+Merr;
Mlo=M0-Merr;

fac0   = 1e6/M0;
fac_up = 1e6/Mlo;
fac_lo = 1e6/Mup;

X    = [1 2 4 8]';
Xfac = X./(100-X);

C0   = fac0 * Xfac;
Clo  = fac_lo * Xfac;
Cup  = fac_up * Xfac;
Cerr = 0.5*(Cup-Clo);
res  = [ X C0 Cerr -log10(1e-3*C0) ];
