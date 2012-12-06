
D_h = 9.3046e-9;
D_w = 5.2291e-9;

D_CO2       = 1.9e-9;
D_Bicarb    = 1.4532e-9;
D_Carbonate = 0.9222e-9;

D_Acetate = 1.0888e-09;

D_Na = 1.3328e-09;
D_K  = 1.9568e-09;
D_Cl = 2.0313e-09;

D_In = 1e-9;

species =
{ 
    { "H+",   1, D_h  },
    { "HO-", -1, D_w  },
    { "Cl-", -1, D_Cl },
    { "Na+",  1, D_Na },
    { "InH",  0, D_In },
    { "In-", -1, D_Cl }
};


ftol = 1e-7;

chemsys =
{
    { "water",  1e-14,     { 1,  "H+" }, { 1, "HO-" } },
    { "indic",  10^(-3.5), { -1, "InH"}, { 1, "H+"  }, { 1, "In-" } }
};

C_In = 0;
indic_total = { C_In, { 1, "InH" }, {1, "In-"} };

ini_bulk =
{
    { 1e-2, { 1, "H+"  } }, -- pH=2
    { 0.1,  { 1, "Na+" } },
    indic_total
};

ini_core = 
{
    { 1e-10, { 1, "H+"}  }, -- pH=10
    { 0.1,   { 1, "Cl-"} },
    indic_total
};


volumes = 1000;
length  = 2e-3;
alpha   = 0.02;
t_run   = 1;
dt_save = 0.01;

-- pH_front = 3.39
pH_front = 6
fit      = true;

