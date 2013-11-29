Dh    = 9.3046e-9;
Dw    = 5.2291e-9;
DNa   = 1.3328e-9;
DCl   = 2.0313e-9;
DK    = 1.9568e-9;
DCO2  = 1.9000e-9;
DHCO3 = 1.4532e-9;
DCO3  = 0.9222e-9;
DAH   = 1.2100e-9;
DAm   = DAH;

-- database of species
species =
{
    { "H+",     1,  Dh  },
    { "HO-",   -1,  Dw  },
    { "AH",     0,  DAH  },
    { "A-",    -1, DAm  },
    { "Na+",    1,  DNa  },
    { "Cl-",   -1,  DCl  },
    { "CO2",    0,  DCO2 },
    { "HCO3-", -1, DHCO3 }
};

-- equations at stake
eqs =
{
    { "water",  1e-14,     { 1, "H+"}, { 1, "HO-" } },
    { "acetic", 10^(-4.76), { 1, "H+" }, { 1, "A-"}, {-1, "AH" } },
    { "K1",     4.6e-7,     { 1, "H+" }, { 1, "HCO3-" }, { -1, "CO2" } }
};

C_acid   = 1e-5;
C_bicarb = 1e-5;

NoNa = { 0, {1, "Na+"  } };
NoCl = { 0, {1, "Cl-"  } };


NoCarbonate = { 0, { 1, "CO2" }, { 1, "HCO3-" } }
NoAcid      = { 0, { 1, "AH" },  { 1, "A-"    } }

ini_left =
{
    { C_acid, {1, "AH"}, {1,"A-" } },
    NoNa,
    NoCl,
    NoCarbonate
};

ini_core =
{
    { 1e-7, { 1, "H+" } },
    NoNa,
    NoCl,
    NoCarbonate
}

right_wall = 0;

ini_right =
{
    NoAcid,
    { C_bicarb, {1,"Na+"} },
    NoCl,
    { C_bicarb, {1, "CO2"}, {1,"HCO3-"} }
}

volumes = 100;
Tmax    = 800;
alpha   = 0.1;
dt      = 0.1;
save    = 1;
length  = 0.01;

