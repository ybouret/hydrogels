species =
{ 
    { "H+", 1, 0.01},
    { "HO-", -1, 0.002},
    { "Cl-", -1, 0.000}
};


ftol = 1e-5;

chemsys =
{
    { "water",  1e-14,     { 1, "H+" }, { 1, "HO-" } }
};

ini_left =
{
    { 1e-5, { 1, "H+" } }
};



volumes=100;
length =2;