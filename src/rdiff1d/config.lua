
Dh  = 9.3046e-9;
Dw  = 5.2291e-9;
DNa = 1.3328e-9;
DCl = 2.0313e-9;
DK  = 1.9568e-9;

-- database of species
species =
{
    { "H+",   1, Dh },
    { "HO-", -1, Dw },
    --{ "InH",  0, 0  },
    --{ "In-", -1, 0  },
    { "Na+",  1, DNa},
    { "Cl-", -1, DCl}
};

-- equations at stake
eqs =
{
    { "water",  1e-14,     { 1, "H+"}, { 1, "HO-" } },
    --{ "color", 10^(-4.8),  { 1, "H+"}, { 1, "In-" }, { -1, "InH" } },
};

Csalt = 0.0;

-- some constraints
Na    = { Csalt, {1,"Na+"} };
Cl    = { Csalt, {1,"Cl-" } };
Indic = { 0, {1,"InH"}, {1,"In-" } };

-- boundary/initial conditions
pH_left  = 5;
pH_core  = 9;
pH_right = 6;

-- initialize: set of side constraints

-- corrected by HCl => Na is given, Cl- is variable
ini_left =
{
   -- Indic,
    Na,
    { 10^(-pH_left), {1, "H+" } }
};

-- initialize: set of core constraints
-- corrected by NaOH => Cl- is given, Na+ is variable
ini_core =
{
   -- Indic,
    Cl,
    { 10^(-pH_core), {1, "H+" } }
};


right_wall = 1;

ini_right =
{
    Na,
    { 10^(-pH_right), {1,"H+"} }
}


volumes = 50;
length  = 0.001; -- in meters

alpha = 0.4;  -- dt Dmax/dx_min^2
Tmax  = 100;  -- run time in seconds
dt    = 0.01; -- required dt
save  = 0.1;  -- in seconds

