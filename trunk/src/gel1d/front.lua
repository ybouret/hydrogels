
-- -----------------------------------------------------------------------------
-- diffusion coefficients
-- -----------------------------------------------------------------------------

D_h = 9.3046e-9;
D_w = 5.2291e-9;


D_Na = 1.3328e-09;
D_K  = 1.9568e-09;
D_Cl = 2.0313e-09;

D_color = 1e-9;
-- -----------------------------------------------------------------------------
-- species
-- -----------------------------------------------------------------------------
species =
{
	{ "H+",     1, D_h },
 	{ "HO-",   -1, D_w },
    { "Na+",    1, D_Na },
    { "Cl-",   -1, D_Cl },
    { "AH",     0, D_color },
    { "A-",    -1, D_color }
};


Ki = 10^(-3.39)

equilibria =
{
    { "water",  1e-14,     { 1, "H+" }, { 1, "HO-" } },
    { "indic",  Ki,        { 1, "H+" }, { 1, "A-"}, {-1,"AH" }}
};

Ca = 1e-4;
-- the program adds the electroneutrality
ini_left =
{
    { 0.1,      { 1, "Na+"} },
    { 0.1+1e-2, { 1, "Cl-" } },
    { Ca,       { 1, "AH"}, {1,"A-"}}
}

ini_core =
{
    { 0.1,     { 1, "Na+"} },
    { 0.1+1e-5,     { 1, "Cl-"} },
    { Ca,    { 1, "AH"}, {1,"A-"}}
}

ini_right = ini_core;

ntop        = 440;     -- mesh: 0..ntop
gel_length  = 0.022;   -- in meters
noRightFlux = true;    -- boundary condition
pH_front    = 3.5;

