N=1000;

alpha=1;

function Cstar(t)
return 0.5;
end

D = 2.0e-9; -- m^2/s

R_ini = 0.1 * 1e-3; -- initial size, in m
R_end = 10  * 1e-3; -- finial  size, in m

t_run  = 120;
dt     = 1e-3;
dt_sav = 1;