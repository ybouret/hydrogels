
-- database of species
species =
{
    { "H+",   1, 1},
    { "HO-", -1, 1},
    { "InH",  0, 1},
    { "In-", -1, 1},
    { "Na+",  1, 1},
    { "Cl-", -1, 1}
};

-- equations at stake
eqs =
{
    { "water",  1e-14,     { 1, "H+"}, { 1, "HO-" } },
    { "color", 10^(-4.8),  { 1, "H+"}, { 1, "In-" }, { -1, "InH" } },
};


-- some constraints
Na    = { 0.1, {1,"Na+"} };
Cl    = { 0.1, {1,"Cl-" } };
Indic = { 0, {1,"InH"}, {1,"In-" } };

-- boundary/initial conditions
pH_side = 2;
pH_core = 10;

-- initialize: set of side constraints

-- corrected by HCl => Na is given, Cl- is variable
ini_side =
{
    Indic,
    Na,
    { 10^(-pH_side), {1, "H+" } }
};

-- initialize: set of core constraints
-- corrected by NaOH => Cl- is given, Na+ is variable
ini_core =
{
    Indic,
    Cl,
    { 10^(-pH_core), {1, "H+" } }
};
