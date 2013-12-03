Dh  = 9.3046e-9;

species =
{
    { "H+", 0, Dh } -- no charge => no electroneutrality
};

eqs =
{
    
}

ini_left =
{
    { 2.0, { 1, "H+" } }
}

ini_core =
{
    { 1.0, { 1, "H+" } }
};

right_wall = 1;

volumes = 100;
length  = 0.01; -- in meters

alpha = 0.4;  -- dt Dmax/dx_min^2
Tmax  = 10;  -- run time in seconds
dt    = 0.05; -- required dt
save  = 0.1;  -- in seconds
