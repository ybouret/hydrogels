1;
ag   = [ 0.5 1 2 4];
D_av = [11.8 11.8 11.8 11.8];
D_lo = D_av - 0.8;
D_hi = D_av + 0.8;

D_ag  = [6.84 6.34 6.07 6.51 ];
D_err = [ 0.01 0.01 0.02 0.01];
D_ag_lo = D_ag - D_err;
D_ag_hi = D_ag + D_err;

rho_av = sqrt(D_ag./D_av); 
rho_lo = sqrt(D_ag_lo ./ D_hi);
rho_hi = sqrt(D_ag_hi ./ D_lo);

eta_av = (2*rho_av) ./ (1+rho_av);
eta_lo = (2*rho_lo) ./ (1+rho_hi);
eta_hi = (2*rho_hi) ./ (1+rho_lo);

eta_err = (eta_hi-eta_lo)/2;
eta_plus = eta_hi - eta_av;
eta_minus = eta_lo - eta_av;

[ ag' eta_av' eta_lo' eta_hi' eta_minus' eta_plus' eta_err' ]
