
Dh  = 9.3046e-9;
Dw  = 5.2291e-9;
DNa = 1.3328e-9;
DCl = 2.0313e-9;
DK  = 1.9568e-9;

DInH  = 0.6140e-9;
DInm  = DInH;

rescale = 1.0

-- database of species
species =
{
    { "H+",     1, Dh   * rescale   },
    { "HO-",   -1, Dw   * rescale   },
    { "InH",    0, DInH * rescale   },
    { "In-",   -1, DInm * rescale   },
    { "Na+",    1, DNa  * rescale   },
    { "Cl-",   -1, DCl  * rescale   },
    { "GelH2+", 1, 0     },
    { "GelH",   0, 0     },
    { "Gel-",  -1, 0     }
};

-- equations at stake

-- type A Gelatine
pKa1 = 5.2;
pKa2 = 11.5;

--type B Gelatine
-- pKa1 = 3.6
-- pKa2 = 7.8

eqs =
{
    { "water",     1e-14,       { 1, "H+"}, { 1, "HO-" } },
    { "color",     10^(-3.39),  { 1, "H+"}, { 1, "In-" }, { -1, "InH" } },
    { "gelatine1", 10^(-pKa1),  { 1, "H+"}, { 1, "GelH"}, { -1, "GelH2+" } },
    { "gelatine2", 10^(-pKa2),  { 1, "H+"}, { 1, "Gel-"}, { -1, "GelH"   } }
};

Csalt = 0.1;
--Cgel  = 0.0000;

-- some constraints
Na       = { Csalt, {1,"Na+"} };
Cl       = { Csalt, {1,"Cl-" } };
Indic    = { 1e-4,  {1,"InH"}, {1,"In-" } };
Gelatine = { Cgel,  {1, "GelH2+"}, {1,"GelH"}, {1,"Gel-"} };

-- boundary/initial conditions
pH_left  = 2;
pH_core  = 10;
pH_right = pH_core;

-- initialize: set of side constraints

-- corrected by HCl => Na is given, Cl- is variable
ini_left =
{
    Indic,
    Na,
    { 10^(-pH_left), {1, "H+" } },
    Gelatine
};

-- initialize: set of core constraints
-- corrected by NaOH => Cl- is given, Na+ is variable
ini_core =
{
    Indic,
    Cl,
    { 10^(-pH_core), {1, "H+" } },
    Gelatine
};


right_wall = 1;

search_output = string.format("Cgel%.2g.dat", Cgel);

