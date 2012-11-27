species =
{ 
    { "H+", 1, 0.01},
    { "HO-", -1, 0.002},
    { "Cl-", -1, 0.000},
    { "Na+",  1, 0.000}
};


ftol = 1e-5;

chemsys =
{
    { "water",  1e-14,     { 1, "H+" }, { 1, "HO-" } }
};

ini_bulk =
{
    { 1e-5, { 1, "H+" } },
    { 0, { 1, "Na+" } }
};

ini_core = 
{
    { 1e-8, { 1, "H+"} },
    { 0, { 1, "Cl-"} }
};


volumes=200;
length =2;