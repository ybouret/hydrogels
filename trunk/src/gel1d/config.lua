
species =
{
	{ "H+", 1, 1 },
	{ "HO-", -1, 1 },
    { "AcH", 0, 1 },
    { "Ac-", -1, 1}
};


equilibria =
{
    { "water",  1e-14,     { 1, "H+" }, { 1, "HO-" } },
    { "acetic", 10^(-4.7), { 1, "H+"},  { 1, "Ac-"}, { -1, "AcH"} }
};


ini_left =
{
    { 10^(-5), { 1, "H+" } },
    { 0.0,     { 1, "AcH" }, { 1, "Ac-" } }
}

ini_right =
{
    { 10^(-7), { 1, "H+" } },
    { 0.0,     { 1, "AcH" }, { 1, "Ac-" } }
}

ini_core =
{
    { 10^(-7), { 1, "H+" } },
    { 0.0,     { 1, "AcH" }, { 1, "Ac-" } }
}

ntop       = 10;  -- mesh: 0..ntop
gel_length = 1e-2; -- in meters
