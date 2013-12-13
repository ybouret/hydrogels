volumes=600; -- #volumes
length=0.02; -- #in meters
alpha=0.2;   -- #max dt D/dx
Dh  = 9.3046e-9;
-- Dw  = 5.2291e-9;
Dw = Dh;

ftol = 1e-3; -- #tolerance for Newton/ODE

t_run  =10;
dt_save=0.1;

iter=1;
