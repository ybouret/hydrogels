
D_h = 9.3046e-9;
D_w = 5.2291e-9;

D_CO2       = 1.9e-9;
D_Bicarb    = 1.4532e-9;
D_Carbonate = 0.9222e-9;

D_Acetate = 1.0888e-09;

D_Na = 1.3328e-09;
D_K  = 1.9568e-09;
D_Cl = 2.0313e-09;


species =
{ 
    { "H+",     1, D_h  },
    { "HO-",   -1, D_w  },
    { "Cl-",   -1, D_Cl },
    { "Na+",    1, D_Na },
    { "K+",     1, D_K  },
    { "CO2",    0, D_CO2 },
    { "HCO3-", -1, D_Bicarb },
    { "CO3--", -2, D_Carbonate}
};


ftol = 1e-7;

chemsys =
{
    { "water",  1e-14,     {  1, "H+" }, { 1, "HO-" } },
    { "K1",     4.6e-7,    {  1, "H+" }, { 1, "HCO3-" }, { -1, "CO2" } },
    { "K2",     4.69e-11,  {  1, "H+" }, { 1, "CO3--" }, { -1, "HCO3-" }}
};

ini_K  = { 5e-3,   { 1, "K+"  } }
ini_Na = { 140e-3, { 1, "Na+" } }
ini_Cl = { 100e-3, { 1, "Cl-" } }

ini_skin =
{
    { 10^(-7.4), { 1, "H+" } },
    ini_K,
    ini_Na,
    ini_Cl
}

ini_bath =
{
    { 10^(-6), { 1, "H+" } },
    ini_K,
    ini_Na,
    ini_Cl
}

-- Geometry

xvolumes = 30;
yvolumes = 20;
Lx       = 5e-3;
Ly       = 10e-3;
alpha    = 0.02;
t_run    = 5;
dt_save  = 0.05;