#include "config.h"
#include "EstimateXFe.h"

struct StructAverageElementMass{
	 char Name[100];
	 double  Value;
};
const int AverageElementMassSize = 86;
struct StructAverageElementMass AverageElementMass[] = {
	{.Name = "H",
	 .Value = 1.0000200000000001,},
	{.Name = "He",
	 .Value = 3.999834,},
	{.Name = "Li",
	 .Value = 6.9241,},
	{.Name = "Be",
	 .Value = 9.0,},
	{.Name = "B",
	 .Value = 10.801,},
	{.Name = "C",
	 .Value = 12.011061999999999,},
	{.Name = "N",
	 .Value = 14.002290000000002,},
	{.Name = "O",
	 .Value = 16.004379000000004,},
	{.Name = "F",
	 .Value = 19.0,},
	{.Name = "Ne",
	 .Value = 20.138910000000003,},
	{.Name = "Na",
	 .Value = 23.0,},
	{.Name = "Mg",
	 .Value = 24.3202,},
	{.Name = "Al",
	 .Value = 27.0,},
	{.Name = "Si",
	 .Value = 28.108604,},
	{.Name = "P",
	 .Value = 31.0,},
	{.Name = "S",
	 .Value = 32.0942,},
	{.Name = "Cl",
	 .Value = 35.4844,},
	{.Name = "Ar",
	 .Value = 36.3086,},
	{.Name = "K",
	 .Value = 39.135889999999996,},
	{.Name = "Ca",
	 .Value = 40.11563,},
	{.Name = "Sc",
	 .Value = 45.0,},
	{.Name = "Ti",
	 .Value = 47.9183,},
	{.Name = "V",
	 .Value = 50.9975,},
	{.Name = "Cr",
	 .Value = 52.05541,},
	{.Name = "Mn",
	 .Value = 55.00000000000001,},
	{.Name = "Fe",
	 .Value = 54.53854000000001,},
	{.Name = "Fe",
	 .Value = 1.37139,},
	{.Name = "Co",
	 .Value = 59.0,},
	{.Name = "Ni",
	 .Value = 58.759575,},
	{.Name = "Cu",
	 .Value = 63.616600000000005,},
	{.Name = "Zn",
	 .Value = 65.4682,},
	{.Name = "Ga",
	 .Value = 69.79784000000001,},
	{.Name = "Ge",
	 .Value = 72.6905,},
	{.Name = "As",
	 .Value = 75.0,},
	{.Name = "Se",
	 .Value = 79.0421,},
	{.Name = "Br",
	 .Value = 79.9862,},
	{.Name = "Kr",
	 .Value = 83.88083999999999,},
	{.Name = "Rb",
	 .Value = 85.58312,},
	{.Name = "Sr",
	 .Value = 87.71136299999999,},
	{.Name = "Y",
	 .Value = 89.0,},
	{.Name = "Zr",
	 .Value = 91.31840000000001,},
	{.Name = "Nb",
	 .Value = 93.0,},
	{.Name = "Mo",
	 .Value = 96.05436999999999,},
	{.Name = "Ru",
	 .Value = 101.15980000000002,},
	{.Name = "Rh",
	 .Value = 103.0,},
	{.Name = "Pd",
	 .Value = 12.626000000000001,},
	{.Name = "Pd",
	 .Value = 93.8851,},
	{.Name = "Ag",
	 .Value = 107.96322,},
	{.Name = "Cd",
	 .Value = 112.50800000000001,},
	{.Name = "In",
	 .Value = 114.91420000000001,},
	{.Name = "Sn",
	 .Value = 118.8077,},
	{.Name = "Sb",
	 .Value = 121.85579999999999,},
	{.Name = "Te",
	 .Value = 127.6984,},
	{.Name = "I",
	 .Value = 127.0,},
	{.Name = "Xe",
	 .Value = 131.2881,},
	{.Name = "Cs",
	 .Value = 133.0,},
	{.Name = "Ba",
	 .Value = 137.42162000000002,},
	{.Name = "La",
	 .Value = 138.99909000000002,},
	{.Name = "Ce",
	 .Value = 140.20986000000002,},
	{.Name = "Pr",
	 .Value = 206.9724,},
	{.Name = "Nd",
	 .Value = 144.33496,},
	{.Name = "Sm",
	 .Value = 150.4481,},
	{.Name = "Eu",
	 .Value = 152.0438,},
	{.Name = "Gd",
	 .Value = 157.3281,},
	{.Name = "Tb",
	 .Value = 159.0,},
	{.Name = "Dy",
	 .Value = 0.23746,},
	{.Name = "Dy",
	 .Value = 162.33407,},
	{.Name = "Ho",
	 .Value = 165.0,},
	{.Name = "Er",
	 .Value = 167.32707,},
	{.Name = "Tm",
	 .Value = 169.0,},
	{.Name = "Yb",
	 .Value = 173.1335,},
	{.Name = "Lu",
	 .Value = 175.028205,},
	{.Name = "Hf",
	 .Value = 178.54163,},
	{.Name = "Ta",
	 .Value = 180.99988000000002,},
	{.Name = "W",
	 .Value = 183.89069999999998,},
	{.Name = "Re",
	 .Value = 186.28676000000002,},
	{.Name = "Os",
	 .Value = 190.29070000000002,},
	{.Name = "Ir",
	 .Value = 192.254,},
	{.Name = "Pt",
	 .Value = 129.14108,},
	{.Name = "Au",
	 .Value = 197.0,},
	{.Name = "Hg",
	 .Value = 200.6297,},
	{.Name = "Tl",
	 .Value = 204.40952,},
	{.Name = "Pb",
	 .Value = 207.34285,},
	{.Name = "Bi",
	 .Value = 209.0,},
	{.Name = "Th",
	 .Value = 231.99999999999997,},
	{.Name = "U",
	 .Value = 237.27134,},
};



struct StructSolarAbundancePattern{
	 char Name[100];
	 double  Value;
};

const int SolarAbundancePatternSize = 83;
struct StructSolarAbundancePattern SolarAbundancePattern[] = {
	{.Name = "H",
	 .Value = 12.0,},
	{.Name = "He",
	 .Value = 10.93,},
	{.Name = "Li",
	 .Value = 1.05,},
	{.Name = "Be",
	 .Value = 1.38,},
	{.Name = "B",
	 .Value = 2.7,},
	{.Name = "C",
	 .Value = 8.43,},
	{.Name = "N",
	 .Value = 7.83,},
	{.Name = "O",
	 .Value = 8.69,},
	{.Name = "F",
	 .Value = 4.56,},
	{.Name = "Ne",
	 .Value = 7.93,},
	{.Name = "Na",
	 .Value = 6.24,},
	{.Name = "Mg",
	 .Value = 7.6,},
	{.Name = "Al",
	 .Value = 6.45,},
	{.Name = "Si",
	 .Value = 7.51,},
	{.Name = "P",
	 .Value = 5.41,},
	{.Name = "S",
	 .Value = 7.12,},
	{.Name = "Cl",
	 .Value = 5.5,},
	{.Name = "Ar",
	 .Value = 6.4,},
	{.Name = "K",
	 .Value = 5.03,},
	{.Name = "Ca",
	 .Value = 6.34,},
	{.Name = "Sc",
	 .Value = 3.15,},
	{.Name = "Ti",
	 .Value = 4.95,},
	{.Name = "V",
	 .Value = 3.93,},
	{.Name = "Cr",
	 .Value = 5.64,},
	{.Name = "Mn",
	 .Value = 5.43,},
	{.Name = "Fe",
	 .Value = 7.5,},
	{.Name = "Co",
	 .Value = 4.99,},
	{.Name = "Ni",
	 .Value = 6.22,},
	{.Name = "Cu",
	 .Value = 4.19,},
	{.Name = "Zn",
	 .Value = 4.56,},
	{.Name = "Ga",
	 .Value = 3.04,},
	{.Name = "Ge",
	 .Value = 3.65,},
	{.Name = "As",
	 .Value = -10.0,},
	{.Name = "Se",
	 .Value = -10.0,},
	{.Name = "Br",
	 .Value = -10.0,},
	{.Name = "Kr",
	 .Value = 3.25,},
	{.Name = "Rb",
	 .Value = 2.52,},
	{.Name = "Sr",
	 .Value = 2.87,},
	{.Name = "Y",
	 .Value = 2.21,},
	{.Name = "Zr",
	 .Value = 2.58,},
	{.Name = "Nb",
	 .Value = 1.46,},
	{.Name = "Mo",
	 .Value = 1.88,},
	{.Name = "Ru",
	 .Value = 1.75,},
	{.Name = "Rh",
	 .Value = 0.91,},
	{.Name = "Pd",
	 .Value = 1.57,},
	{.Name = "Ag",
	 .Value = 0.94,},
	{.Name = "Cd",
	 .Value = -10.0,},
	{.Name = "In",
	 .Value = 0.8,},
	{.Name = "Sn",
	 .Value = 2.04,},
	{.Name = "Sb",
	 .Value = -10.0,},
	{.Name = "Te",
	 .Value = -10.0,},
	{.Name = "I",
	 .Value = -10.0,},
	{.Name = "Xe",
	 .Value = 2.2,},
	{.Name = "Cs",
	 .Value = -10.0,},
	{.Name = "Ba",
	 .Value = 2.18,},
	{.Name = "La",
	 .Value = 1.1,},
	{.Name = "Ce",
	 .Value = 1.58,},
	{.Name = "Pr",
	 .Value = 0.72,},
	{.Name = "Nd",
	 .Value = 1.42,},
	{.Name = "Sm",
	 .Value = 0.96,},
	{.Name = "Eu",
	 .Value = 0.52,},
	{.Name = "Gd",
	 .Value = 1.07,},
	{.Name = "Tb",
	 .Value = 0.3,},
	{.Name = "Dy",
	 .Value = 1.1,},
	{.Name = "Ho",
	 .Value = 0.48,},
	{.Name = "Er",
	 .Value = 0.92,},
	{.Name = "Tm",
	 .Value = 0.1,},
	{.Name = "Yb",
	 .Value = 0.84,},
	{.Name = "Lu",
	 .Value = 0.1,},
	{.Name = "Hf",
	 .Value = 0.85,},
	{.Name = "Ta",
	 .Value = -10.0,},
	{.Name = "W",
	 .Value = 0.85,},
	{.Name = "Re",
	 .Value = -10.0,},
	{.Name = "Os",
	 .Value = 1.4,},
	{.Name = "Ir",
	 .Value = 1.38,},
	{.Name = "Pt",
	 .Value = -10.0,},
	{.Name = "Au",
	 .Value = 0.92,},
	{.Name = "Hg",
	 .Value = -10.0,},
	{.Name = "Tl",
	 .Value = 0.9,},
	{.Name = "Pb",
	 .Value = 1.75,},
	{.Name = "Bi",
	 .Value = -10.0,},
	{.Name = "Th",
	 .Value = 0.02,},
	{.Name = "U",
	 .Value = -10.0,},
};

struct StructFittingCoef{
	 char Name[100];
	 double  Coef[4];
};


const int FittingCoefSagaFullSize = 16;
struct StructFittingCoef FittingCoefSagaFull[] = {
	{.Name = "C",
	 .Coef[0] = 0.07244646097660129,
	 .Coef[1] = -0.38232437645098766,
	 .Coef[2] = -0.22640313008025711,
	 .Coef[3] = -0.04876283187993332,
	},
	{.Name = "N",
	 .Coef[0] = -1.9633904300083036,
	 .Coef[1] = -2.7255280230440073,
	 .Coef[2] = -0.8877989753263729,
	 .Coef[3] = -0.1067796779037367,
	},
	{.Name = "O",
	 .Coef[0] = 0.08154680418775746,
	 .Coef[1] = 0.02374188384530344,
	 .Coef[2] = -0.021999618808812082,
	 .Coef[3] = -0.013279367851578914,
	},
	{.Name = "Na",
	 .Coef[0] = 0.08154680418775746,
	 .Coef[1] = 0.02374188384530344,
	 .Coef[2] = -0.021999618808812082,
	 .Coef[3] = -0.013279367851578914,
	},
	{.Name = "Mg",
	 .Coef[0] = 0.03386681373050431,
	 .Coef[1] = -0.3517163631997052,
	 .Coef[2] = -0.14717761672705915,
	 .Coef[3] = -0.022912801484105363,
	},
	{.Name = "Si",
	 .Coef[0] = 0.09813587375019486,
	 .Coef[1] = -0.2860955675285689,
	 .Coef[2] = -0.09576470764350152,
	 .Coef[3] = -0.012500320944729491,
	},
	{.Name = "S",
	 .Coef[0] = 0.05605155220960065,
	 .Coef[1] = 0.06957786714485809,
	 .Coef[2] = 0.02743337687249802,
	 .Coef[3] = 0.0020124385155206857,
	},
	{.Name = "Ca",
	 .Coef[0] = 0.047820906064949285,
	 .Coef[1] = -0.23433298546460374,
	 .Coef[2] = -0.053973641693984366,
	 .Coef[3] = -0.003701738810780395,
	},
	{.Name = "Ni",
	 .Coef[0] = 0.05605155220960065,
	 .Coef[1] = 0.06957786714485809,
	 .Coef[2] = 0.02743337687249802,
	 .Coef[3] = 0.0020124385155206857,
	},
	{.Name = "Zn",
	 .Coef[0] = 0.09697378398242171,
	 .Coef[1] = -0.028910478769483078,
	 .Coef[2] = -0.03374891111412462,
	 .Coef[3] = -0.016920968370606568,
	},
	{.Name = "Ba",
	 .Coef[0] = -0.039927756642270315,
	 .Coef[1] = -0.03294114304357123,
	 .Coef[2] = 0.06699439407706859,
	 .Coef[3] = 0.033939500697402196,
	},
	{.Name = "Sr",
	 .Coef[0] = -0.5012114226749763,
	 .Coef[1] = -1.131553974529649,
	 .Coef[2] = -0.5142702192842246,
	 .Coef[3] = -0.057566197338619694,
	},
	{.Name = "Eu",
	 .Coef[0] = 0.03599806247486634,
	 .Coef[1] = -0.4155936192128074,
	 .Coef[2] = -0.1085800221335415,
	 .Coef[3] = -0.0074883828155753234,
	},
//#Note: Since there is no available data of Ne in SAGA database, the coefs of Mg are used.
	{.Name = "Ne",
	 .Coef[0] = 0.03386681373050431,
	 .Coef[1] = -0.3517163631997052,
	 .Coef[2] = -0.14717761672705915,
	 .Coef[3] = -0.022912801484105363,
	},
//#Note: Since there is (almost) no available data of P in SAGA database, the flat is adopted.
	{.Name = "P",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},
};

const int FittingCoefSagaCompactSize = 22;
struct StructFittingCoef FittingCoefSagaCompact[] = {
	{.Name = "C",
	 .Coef[0] = 0.116114,
	 .Coef[1] = -0.1074, 
	 .Coef[2] = -0.10963,
	 .Coef[3] = -0.0253474,
	},
	{.Name = "N", // Almost no data for [Fe/H] > -1
	 .Coef[0] = 0.302988, 
	 .Coef[1] = 0.199542, 
	 .Coef[2] = 0.134288,
	 .Coef[3] = 0.0111476,
	},
	{.Name = "O",
	 .Coef[0] = 0.153876, 
	 .Coef[1] = -0.801999, 
	 .Coef[2] = -0.48025, 
	 .Coef[3] = -0.093347,
	},
	{.Name = "Na",
	 .Coef[0] = 0.0576812, 
	 .Coef[1] = -0.0155756, 
	 .Coef[2] = -0.0264432, 
	 .Coef[3] = -0.00866148,
	},
	{.Name = "Mg",
	 .Coef[0] = 0.0464037, 
	 .Coef[1] = -0.342724, 
	 .Coef[2] = -0.143449, 
	 .Coef[3] = -0.0212719,
	},

	{.Name = "Al",
	 .Coef[0] = 0.0821379, 
	 .Coef[1] = -0.232196, 
	 .Coef[2] = -0.416064, 
	 .Coef[3] = -0.0838226,
	},
	{.Name = "Si",
	 .Coef[0] = 0.188843, 
	 .Coef[1] = -0.317837, 
	 .Coef[2] = -0.0722834, 
	 .Coef[3] = -0.00471598,
	},
	{.Name = "P", // The data points are not enough.
	 // .Coef[0] = 0.370647, 
	 // .Coef[1] = 0.272895, 
	 // .Coef[2] = 0.739934, 
	 // .Coef[3] = 0.333831
	 .Coef[0] = 0.0,
	 .Coef[1] = 0.0,
	 .Coef[2] = 0.0,
	 .Coef[3] = 0.0,
	},
	{.Name = "S",
	 .Coef[0] = 0.188843, 
	 .Coef[1] = -0.317837,
	 .Coef[2] = -0.0722834, 
	 .Coef[3] = -0.00471598,
	},
	{.Name = "K",
	 .Coef[0] = 0.189191, 
	 .Coef[1] = -0.197355, 
	 .Coef[2] = -0.0981321, 
	 .Coef[3] = -0.0211905,
	},
	{.Name = "Ca",
	 .Coef[0] = 0.0560119, 
	 .Coef[1] = -0.227, 
	 .Coef[2] = -0.0565123, 
	 .Coef[3] = -0.00455427,
	},
	{.Name = "Ni",
	 .Coef[0] = 0.0541504, 
	 .Coef[1] = 0.0686141, 
	 .Coef[2] = 0.0205764, 
	 .Coef[3] = 3.14348e-5,
	},
	{.Name = "Co",// Almost no data for [Fe/H] > -1
	 .Coef[0] = -0.00098173, 
	 .Coef[1] = -0.427932, 
	 .Coef[2] = -0.393026, 
	 .Coef[3] = -0.0888319,
	},
	{.Name = "Cr",
	 .Coef[0] = 0.0529467, 
	 .Coef[1] = -0.0300297, 
	 .Coef[2] = -0.0763834, 
	 .Coef[3] = -0.0112343,
	},
	{.Name = "Mn",// Almost no data for [Fe/H] > -1
	 .Coef[0] = 0.638244, 
	 .Coef[1] = 1.53493, 
	 .Coef[2] = 0.718141, 
	 .Coef[3] = 0.108057,
	},
	{.Name = "Sr",// Almost no data for [Fe/H] > -1
	 .Coef[0] = 0.198974, 
	 .Coef[1] = -0.196958, 
	 .Coef[2] = -0.396965, 
	 .Coef[3] = -0.108067,
	},
	{.Name = "Y",
	 .Coef[0] = 0.0132705, 
	 .Coef[1] = -0.0321223, 
	 .Coef[2] = -0.065291, 
	 .Coef[3] = -0.00988122,
	},
	{.Name = "Ti",
	 .Coef[0] = 0.0764876, 
	 .Coef[1] = -0.223838, 
	 .Coef[2] = -0.0764549, 
	 .Coef[3] = -0.00995204,
	},
	{.Name = "Zn",
	 .Coef[0] = 0.114945, 
	 .Coef[1] = -0.102028, 
	 .Coef[2] = -0.133357, 
	 .Coef[3] = -0.0482015,
	},
	{.Name = "Ba",// Almost no data for [Fe/H] > -1
	 .Coef[0] = -0.130867,
	 .Coef[1] = -0.223465, 
	 .Coef[2] = -0.283106, 
	 .Coef[3] = -0.0522968,
	},
	{.Name = "Eu",// Almost no data for [Fe/H] > -1
	 .Coef[0] = 0.258222, 
	 .Coef[1] = -0.177767, 
	 .Coef[2] = -0.183894, 
	 .Coef[3] = -0.0479506,
	},
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},
};

const int FittingCoefGalahSize = 23;
struct StructFittingCoef FittingCoefGalah[] = {
	{.Name = "Li",
	 .Coef[0] = 0.635258, 
	 .Coef[1] = -1.05084, 
	 .Coef[2] = -1.65439, 
	 .Coef[3] = -0.749118,
	},
	{.Name = "C",
	 .Coef[0] = -0.0505651, 
	 .Coef[1] = -0.0478904,
	 .Coef[2] = -0.244249,
	 .Coef[3] = -0.112594,
	},
	{.Name = "O",
	 .Coef[0] = 0.0323147, 
	 .Coef[1] = -0.566089,
	 .Coef[2] = 0.204004,
	 .Coef[3] = 0.233069,
	},
	{.Name = "Na",
	 .Coef[0] = 0.0582172,
	 .Coef[1] = -0.0495983,
	 .Coef[2] = 0.11599,
	 .Coef[3] = 0.115027,
	},
	{.Name = "Mg",
	 .Coef[0] = 0.0808625,
	 .Coef[1] = -0.156105,
	 .Coef[2] = 0.305877,
	 .Coef[3] = 0.201449,
	},
	{.Name = "Al",
	 .Coef[0] = 0.000705132,
	 .Coef[1] = -0.0445434,
	 .Coef[2] = 0.32374,
	 .Coef[3] = 0.192471,
	},
	{.Name = "Si",
	 .Coef[0] = 0.0277335,
	 .Coef[1] = -0.0646924,
	 .Coef[2] = 0.359036,
	 .Coef[3] = 0.19908,
	},
	{.Name = "K",
	 .Coef[0] = 0.144495,
	 .Coef[1] = -0.540476,
	 .Coef[2] = -0.220491,
	 .Coef[3] = 0.0358996,
	},
	{.Name = "Ca",
	 .Coef[0] = 0.0631838, 
	 .Coef[1] = -0.220376,
	 .Coef[2] = 0.101468,
	 .Coef[3] = 0.114314,
	},
	{.Name = "Sc",
	 .Coef[0] = 0.0826671,
	 .Coef[1] = -0.19594,
	 .Coef[2] = -0.0955347,
	 .Coef[3] = 0.0127929,
	},
	{.Name = "Ti",
	 .Coef[0] = 0.0315109,
	 .Coef[1] = -0.157515,
	 .Coef[2] = 0.188432,
	 .Coef[3] = 0.124102,
	},
	{.Name = "V",
	 .Coef[0] = 0.234487,
	 .Coef[1] = 0.00829416, 
	 .Coef[2] = 0.202581,
	 .Coef[3] = 0.157872,
	},
	{.Name = "Cr",
	 .Coef[0] = 0.0202264, 
	 .Coef[1] = 0.0395926, 
	 .Coef[2] = -0.218872,
	 .Coef[3] = -0.108495,
	},
	{.Name = "Mn",
	 .Coef[0] = 0.0132889,
	 .Coef[1] = 0.247854,
	 .Coef[2] = -0.302151, 
	 .Coef[3] = -0.146052,
	},
	{.Name = "Co",
	 .Coef[0] = -0.0942132, 
	 .Coef[1] = 0.111552,
	 .Coef[2] = 0.362292,
	 .Coef[3] = 0.138045,
	},
	{.Name = "Ni",
	 .Coef[0] = 0.152981, 
	 .Coef[1] = 0.165406, 
	 .Coef[2] = 0.260255, 
	 .Coef[3] = 0.144395,
	},
	{.Name = "Cu",
	 .Coef[0] = 0.0388541, 
	 .Coef[1] = 0.167546, 
	 .Coef[2] = 0.0296788, 
	 .Coef[3] = 0.0629747,
	},
	{.Name = "Zn",
	 .Coef[0] = 0.0188398,
	 .Coef[1] = -0.0659812,
	 .Coef[2] = 0.272219,
	 .Coef[3] = 0.168769,
	},
	{.Name = "Y",
	 .Coef[0] = 0.0847561,
	 .Coef[1] = -0.118752,
	 .Coef[2] = -0.267887,
	 .Coef[3] = -0.0966515,
	},
	{.Name = "Ba",
	 .Coef[0] = 0.0733723, 
	 .Coef[1] = -0.2918,
	 .Coef[2] = -0.546874,
	 .Coef[3] = -0.220023,
	},
	{.Name = "La",
	 .Coef[0] = 0.0877738,
	 .Coef[1] = -0.11566, 
	 .Coef[2] = -0.175408,
	 .Coef[3] = -0.0936792,
	},
	{.Name = "Eu",
	 .Coef[0] = 0.0956002,
	 .Coef[1] = -0.650759, 
	 .Coef[2] = -0.405441, 
	 .Coef[3] = -0.062431,
	},
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},

};

const int FittingCoefGalahApogeeSize = 27;
struct StructFittingCoef FittingCoefGalahApogee[] = {
	{.Name = "Li",
	 .Coef[0] = 0.635258, 
	 .Coef[1] = -1.05084, 
	 .Coef[2] = -1.65439, 
	 .Coef[3] = -0.749118,
	},
	{.Name = "C",
	 .Coef[0] = -0.0505651, 
	 .Coef[1] = -0.0478904,
	 .Coef[2] = -0.244249,
	 .Coef[3] = -0.112594,
	},
    {.Name = "N", // from APOGEE
        .Coef[0] = 0.248618, 
        .Coef[1] = 0.190538, 
        .Coef[2] = 0.218319,
        .Coef[3] = 0.0575116,
    },
	{.Name = "O",
	 .Coef[0] = 0.0323147, 
	 .Coef[1] = -0.566089,
	 .Coef[2] = 0.204004,
	 .Coef[3] = 0.233069,
	},
	{.Name = "Na",
	 .Coef[0] = 0.0582172,
	 .Coef[1] = -0.0495983,
	 .Coef[2] = 0.11599,
	 .Coef[3] = 0.115027,
	},
	{.Name = "Mg",
	 .Coef[0] = 0.0808625,
	 .Coef[1] = -0.156105,
	 .Coef[2] = 0.305877,
	 .Coef[3] = 0.201449,
	},
	{.Name = "Al",
	 .Coef[0] = 0.000705132,
	 .Coef[1] = -0.0445434,
	 .Coef[2] = 0.32374,
	 .Coef[3] = 0.192471,
	},
	{.Name = "Si",
	 .Coef[0] = 0.0277335,
	 .Coef[1] = -0.0646924,
	 .Coef[2] = 0.359036,
	 .Coef[3] = 0.19908,
	},
    {.Name = "P", // from APOGEE
        .Coef[0] = -0.0453955, 
        .Coef[1] = -0.140561, 
        .Coef[2] = -0.218989, 
        .Coef[3] = -0.0798785,
    },
    {.Name = "S", // from APOGEE
        .Coef[0] = -0.00587101,
        .Coef[1] = -0.197565, 
        .Coef[2] = 0.282541, 
        .Coef[3] = 0.138067,
    },
	{.Name = "K",
	 .Coef[0] = 0.144495,
	 .Coef[1] = -0.540476,
	 .Coef[2] = -0.220491,
	 .Coef[3] = 0.0358996,
	},
	{.Name = "Ca",
	 .Coef[0] = 0.0631838, 
	 .Coef[1] = -0.220376,
	 .Coef[2] = 0.101468,
	 .Coef[3] = 0.114314,
	},
	{.Name = "Sc",
	 .Coef[0] = 0.0826671,
	 .Coef[1] = -0.19594,
	 .Coef[2] = -0.0955347,
	 .Coef[3] = 0.0127929,
	},
	{.Name = "Ti",
	 .Coef[0] = 0.0315109,
	 .Coef[1] = -0.157515,
	 .Coef[2] = 0.188432,
	 .Coef[3] = 0.124102,
	},
	{.Name = "V",
	 .Coef[0] = 0.234487,
	 .Coef[1] = 0.00829416, 
	 .Coef[2] = 0.202581,
	 .Coef[3] = 0.157872,
	},
	{.Name = "Cr",
	 .Coef[0] = 0.0202264, 
	 .Coef[1] = 0.0395926, 
	 .Coef[2] = -0.218872,
	 .Coef[3] = -0.108495,
	},
	{.Name = "Mn",
	 .Coef[0] = 0.0132889,
	 .Coef[1] = 0.247854,
	 .Coef[2] = -0.302151, 
	 .Coef[3] = -0.146052,
	},
	{.Name = "Co",
	 .Coef[0] = -0.0942132, 
	 .Coef[1] = 0.111552,
	 .Coef[2] = 0.362292,
	 .Coef[3] = 0.138045,
	},
	{.Name = "Ni",
	 .Coef[0] = 0.152981, 
	 .Coef[1] = 0.165406, 
	 .Coef[2] = 0.260255, 
	 .Coef[3] = 0.144395,
	},
	{.Name = "Cu",
	 .Coef[0] = 0.0388541, 
	 .Coef[1] = 0.167546, 
	 .Coef[2] = 0.0296788, 
	 .Coef[3] = 0.0629747,
	},
	{.Name = "Zn",
	 .Coef[0] = 0.0188398,
	 .Coef[1] = -0.0659812,
	 .Coef[2] = 0.272219,
	 .Coef[3] = 0.168769,
	},
	{.Name = "Y",
	 .Coef[0] = 0.0847561,
	 .Coef[1] = -0.118752,
	 .Coef[2] = -0.267887,
	 .Coef[3] = -0.0966515,
	},
	{.Name = "Sr", // Assume that Sr has the same trend of Y
	 .Coef[0] = 0.0847561,
	 .Coef[1] = -0.118752,
	 .Coef[2] = -0.267887,
	 .Coef[3] = -0.0966515,
	},
	{.Name = "Ba",
	 .Coef[0] = 0.0733723, 
	 .Coef[1] = -0.2918,
	 .Coef[2] = -0.546874,
	 .Coef[3] = -0.220023,
	},
	{.Name = "La",
	 .Coef[0] = 0.0877738,
	 .Coef[1] = -0.11566, 
	 .Coef[2] = -0.175408,
	 .Coef[3] = -0.0936792,
	},
	{.Name = "Eu",
	 .Coef[0] = 0.0956002,
	 .Coef[1] = -0.650759, 
	 .Coef[2] = -0.405441, 
	 .Coef[3] = -0.062431,
	},
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},

};

const int FittingCoefGalahSagaFullSize = 27;
struct StructFittingCoef FittingCoefGalahSagaFull[] = {
	{.Name = "Li",
	 .Coef[0] = 0.635258, 
	 .Coef[1] = -1.05084, 
	 .Coef[2] = -1.65439, 
	 .Coef[3] = -0.749118,
	},
	{.Name = "C",
	 .Coef[0] = -0.0505651, 
	 .Coef[1] = -0.0478904,
	 .Coef[2] = -0.244249,
	 .Coef[3] = -0.112594,
	},
	{.Name = "O",
	 .Coef[0] = 0.0323147, 
	 .Coef[1] = -0.566089,
	 .Coef[2] = 0.204004,
	 .Coef[3] = 0.233069,
	},
	{.Name = "Na",
	 .Coef[0] = 0.0582172,
	 .Coef[1] = -0.0495983,
	 .Coef[2] = 0.11599,
	 .Coef[3] = 0.115027,
	},
	{.Name = "Mg",
	 .Coef[0] = 0.0808625,
	 .Coef[1] = -0.156105,
	 .Coef[2] = 0.305877,
	 .Coef[3] = 0.201449,
	},
	{.Name = "Al",
	 .Coef[0] = 0.000705132,
	 .Coef[1] = -0.0445434,
	 .Coef[2] = 0.32374,
	 .Coef[3] = 0.192471,
	},
	{.Name = "Si",
	 .Coef[0] = 0.0277335,
	 .Coef[1] = -0.0646924,
	 .Coef[2] = 0.359036,
	 .Coef[3] = 0.19908,
	},
	{.Name = "K",
	 .Coef[0] = 0.144495,
	 .Coef[1] = -0.540476,
	 .Coef[2] = -0.220491,
	 .Coef[3] = 0.0358996,
	},
	{.Name = "Ca",
	 .Coef[0] = 0.0631838, 
	 .Coef[1] = -0.220376,
	 .Coef[2] = 0.101468,
	 .Coef[3] = 0.114314,
	},
	{.Name = "Sc",
	 .Coef[0] = 0.0826671,
	 .Coef[1] = -0.19594,
	 .Coef[2] = -0.0955347,
	 .Coef[3] = 0.0127929,
	},
	{.Name = "Ti",
	 .Coef[0] = 0.0315109,
	 .Coef[1] = -0.157515,
	 .Coef[2] = 0.188432,
	 .Coef[3] = 0.124102,
	},
	{.Name = "V",
	 .Coef[0] = 0.234487,
	 .Coef[1] = 0.00829416, 
	 .Coef[2] = 0.202581,
	 .Coef[3] = 0.157872,
	},
	{.Name = "Cr",
	 .Coef[0] = 0.0202264, 
	 .Coef[1] = 0.0395926, 
	 .Coef[2] = -0.218872,
	 .Coef[3] = -0.108495,
	},
	{.Name = "Mn",
	 .Coef[0] = 0.0132889,
	 .Coef[1] = 0.247854,
	 .Coef[2] = -0.302151, 
	 .Coef[3] = -0.146052,
	},
	{.Name = "Co",
	 .Coef[0] = -0.0942132, 
	 .Coef[1] = 0.111552,
	 .Coef[2] = 0.362292,
	 .Coef[3] = 0.138045,
	},
	{.Name = "Ni",
	 .Coef[0] = 0.152981, 
	 .Coef[1] = 0.165406, 
	 .Coef[2] = 0.260255, 
	 .Coef[3] = 0.144395,
	},
	{.Name = "Cu",
	 .Coef[0] = 0.0388541, 
	 .Coef[1] = 0.167546, 
	 .Coef[2] = 0.0296788, 
	 .Coef[3] = 0.0629747,
	},
	{.Name = "Zn",
	 .Coef[0] = 0.0188398,
	 .Coef[1] = -0.0659812,
	 .Coef[2] = 0.272219,
	 .Coef[3] = 0.168769,
	},
	{.Name = "Y",
	 .Coef[0] = 0.0847561,
	 .Coef[1] = -0.118752,
	 .Coef[2] = -0.267887,
	 .Coef[3] = -0.0966515,
	},
	{.Name = "Ba",
	 .Coef[0] = 0.0733723, 
	 .Coef[1] = -0.2918,
	 .Coef[2] = -0.546874,
	 .Coef[3] = -0.220023,
	},
	{.Name = "La",
	 .Coef[0] = 0.0877738,
	 .Coef[1] = -0.11566, 
	 .Coef[2] = -0.175408,
	 .Coef[3] = -0.0936792,
	},
	{.Name = "Eu",
	 .Coef[0] = 0.0956002,
	 .Coef[1] = -0.650759, 
	 .Coef[2] = -0.405441, 
	 .Coef[3] = -0.062431,
	},
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},
// from SAGA
	{.Name = "N",
	 .Coef[0] = -1.9633904300083036,
	 .Coef[1] = -2.7255280230440073,
	 .Coef[2] = -0.8877989753263729,
	 .Coef[3] = -0.1067796779037367,
	},
	{.Name = "Ne",
	 .Coef[0] = 0.03386681373050431,
	 .Coef[1] = -0.3517163631997052,
	 .Coef[2] = -0.14717761672705915,
	 .Coef[3] = -0.022912801484105363,
	},
	{.Name = "S",
	 .Coef[0] = 0.05605155220960065,
	 .Coef[1] = 0.06957786714485809,
	 .Coef[2] = 0.02743337687249802,
	 .Coef[3] = 0.0020124385155206857,
	},
	{.Name = "Sr",
	 .Coef[0] = -0.5012114226749763,
	 .Coef[1] = -1.131553974529649,
	 .Coef[2] = -0.5142702192842246,
	 .Coef[3] = -0.057566197338619694,
	},

};

const int FittingCoefApogeeSize = 20;
struct StructFittingCoef FittingCoefApogee[] = {
    {.Name = "C",
        .Coef[0] = -0.0443456,
        .Coef[1] = -0.0189738, 
        .Coef[2] = 0.105897, 
        .Coef[3] = 0.0954529,
    },
    {.Name = "Cl",
        .Coef[0] = -0.092206,
        .Coef[1] = -0.172721,
        .Coef[2] = 0.0983815,
        .Coef[3] = 0.0905185,
    },
    {.Name = "N",
        .Coef[0] = 0.248618, 
        .Coef[1] = 0.190538, 
        .Coef[2] = 0.218319,
        .Coef[3] = 0.0575116,
    },
    {.Name = "O",
        .Coef[0] = 0.0206403,
        .Coef[1] = -0.183682, 
        .Coef[2] = 0.131491, 
        .Coef[3] = 0.0761359,
    },
    {.Name = "Na",
        .Coef[0] = 0.0203534, 
        .Coef[1] = 0.26372, 
        .Coef[2] = 0.426456, 
        .Coef[3] = 0.0722008,
    },
    {.Name = "Mg",
        .Coef[0] = 0.0222359, 
        .Coef[1] = -0.121263, 
        .Coef[2] = 0.238401, 
        .Coef[3] = 0.119879,
    },
    {.Name = "Al",
        .Coef[0] = 0.0522939, 
        .Coef[1] = 0.0192999, 
        .Coef[2] = -0.107505, 
        .Coef[3] = -0.0282489,
    },
    {.Name = "Si",
        .Coef[0] = 0.00979022, 
        .Coef[1] = -0.132581, 
        .Coef[2] = 0.15076, 
        .Coef[3] = 0.0702173,
    },
    {.Name = "P",
        .Coef[0] = -0.0453955, 
        .Coef[1] = -0.140561, 
        .Coef[2] = -0.218989, 
        .Coef[3] = -0.0798785,
    },
    {.Name = "S",
        .Coef[0] = -0.00587101,
        .Coef[1] = -0.197565, 
        .Coef[2] = 0.282541, 
        .Coef[3] = 0.138067,
    },
    {.Name = "K",
        .Coef[0] = 0.0151212, 
        .Coef[1] = -0.102125,
        .Coef[2] = 0.106084,
        .Coef[3] = 0.0582489,
    },
    {.Name = "Ca",
        .Coef[0] = 0.00129988, 
        .Coef[1] = -0.0978388, 
        .Coef[2] = 0.12864, 
        .Coef[3] = 0.075861,
    },
    {.Name = "Ti",
        .Coef[0] = 0.0220352, 
        .Coef[1] = 0.0743516, 
        .Coef[2] = 0.150817,
        .Coef[3] = 0.0531193,
    },
    {.Name = "Tiii",
        .Coef[0] = 0.0449238, 
        .Coef[1] = -0.473386, 
        .Coef[2] = -0.0405162, 
        .Coef[3] = 0.0667414,
    },
    {.Name = "V",
        .Coef[0] = 0.0101649, 
        .Coef[1] = 0.1553,
        .Coef[2] = 0.1223, 
        .Coef[3] = 0.027049,
    },
    {.Name = "Cr",
        .Coef[0] = -0.13226,
        .Coef[1] = -0.213362, 
        .Coef[2] = 0.107713, 
        .Coef[3] = 0.0760279,
    },
    {.Name = "Mn",
        .Coef[0] = -0.0159211, 
        .Coef[1] = 0.250874,
        .Coef[2] = 0.0838985, 
        .Coef[3] = 0.00645244,
    },
    {.Name = "Co",
        .Coef[0] = 0.0586986, 
        .Coef[1] = 0.100623, 
        .Coef[2] = 0.246905, 
        .Coef[3] = 0.105207,
    },
    {.Name = "Ni",
        .Coef[0] = 0.0145689, 
        .Coef[1] = -0.02149, 
        .Coef[2] = 0.11953, 
        .Coef[3] = 0.0530565,
    },
	{.Name = "Fe",
	 .Coef[0] = 0.e0,
	 .Coef[1] = 0.e0,
	 .Coef[2] = 0.e0,
	 .Coef[3] = 0.e0,
	},
};

static int GetIndex(const char ElementName[], const int mode, const int DataType){

    if(mode == 0){ // Coef
        if(DataType == XFeFeH_SAGAFull){
            for(int i=0;i<FittingCoefSagaFullSize;i++){
                if(strncmp(ElementName,FittingCoefSagaFull[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        } else if(DataType == XFeFeH_SAGACompact){
            for(int i=0;i<FittingCoefSagaCompactSize;i++){
                if(strncmp(ElementName,FittingCoefSagaCompact[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        } else if(DataType == XFeFeH_GALAH){
            for(int i=0;i<FittingCoefGalahSize;i++){
                if(strncmp(ElementName,FittingCoefGalah[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        } else if(DataType == XFeFeH_GALAH_APOGEE){
            for(int i=0;i<FittingCoefGalahApogeeSize;i++){
                if(strncmp(ElementName,FittingCoefGalahApogee[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        } else if(DataType == XFeFeH_GALAH_SAGAFull){
            for(int i=0;i<FittingCoefGalahSagaFullSize;i++){
                if(strncmp(ElementName,FittingCoefGalahSagaFull[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        } else if(DataType == XFeFeH_APOGEE){
            for(int i=0;i<FittingCoefApogeeSize;i++){
                if(strncmp(ElementName,FittingCoefApogee[i].Name,MaxCharactersInLine-1) == 0){
                    return i;
                }
            }
        }
        return NONE;
    } else if(mode == 1){ // AverageMass
        for(int i=0;i<AverageElementMassSize;i++){
            if(strncmp(ElementName,AverageElementMass[i].Name,MaxCharactersInLine-1) == 0){
                return i;
            }
        }
        return NONE;
    } else if(mode == 2){ // SolarAbundance
        for(int i=0;i<SolarAbundancePatternSize;i++){
            if(strncmp(ElementName,SolarAbundancePattern[i].Name,MaxCharactersInLine-1) == 0){
                return i;
            }
        }
        return NONE;
    } else {
        fprintf(stderr,"Mode incorrect\n");
        assert(mode < 3);
        assert(-1 < mode);
    }
    return NONE;
}

static double FitModel(const int Index, const double FeH, const int DataType){

    assert(Index!=NONE);

    double LogFeH = log10(FeH);
    if(DataType == XFeFeH_SAGAFull){

        LogFeH = fmin(0.5,fmax(LogFeH,-3));

        double f = FittingCoefSagaFull[Index].Coef[0]+
                   FittingCoefSagaFull[Index].Coef[1]*LogFeH+
                   FittingCoefSagaFull[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefSagaFull[Index].Coef[3]*CUBE(LogFeH);
        if(strncmp("N",FittingCoefSagaFull[Index].Name,MaxCharactersInLine-1) == 0){
            return fmax(f,0.0);
        } else if(strncmp("Sr",FittingCoefSagaFull[Index].Name,MaxCharactersInLine-1) == 0){
            return fmax(f,0.0);
        } else {
            return f;
        }
    } else if(DataType == XFeFeH_SAGACompact){

        LogFeH = fmin(0.5,fmax(LogFeH,-3));

        double f = FittingCoefSagaCompact[Index].Coef[0]+
                   FittingCoefSagaCompact[Index].Coef[1]*LogFeH+
                   FittingCoefSagaCompact[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefSagaCompact[Index].Coef[3]*CUBE(LogFeH);

        if(LogFeH > -1) {
            if(strncmp("N",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0 ||
              strncmp("Co",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0 ||
              strncmp("Mn",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0 ||
              strncmp("Sr",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0 ||
              strncmp("Ba",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0 ||
              strncmp("Eu",FittingCoefSagaCompact[Index].Name,MaxCharactersInLine-1) == 0){
                return fmax(f,0.0);
            } else {
                return f;
            }
        } else {
            return f;
        }
    } else if(DataType == XFeFeH_GALAH){
        LogFeH = fmin(0.5,fmax(LogFeH,-2));

        double f = FittingCoefGalah[Index].Coef[0]+
                   FittingCoefGalah[Index].Coef[1]*LogFeH+
                   FittingCoefGalah[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefGalah[Index].Coef[3]*CUBE(LogFeH);
        return f;
    } else if(DataType == XFeFeH_APOGEE){

        LogFeH = fmin(0.5,fmax(LogFeH,-2));

        double f = FittingCoefApogee[Index].Coef[0]+
                   FittingCoefApogee[Index].Coef[1]*LogFeH+
                   FittingCoefApogee[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefApogee[Index].Coef[3]*CUBE(LogFeH);
        if(strncmp("Na",FittingCoefApogee[Index].Name,MaxCharactersInLine-1) == 0){
            LogFeH = fmin(0.5,fmax(LogFeH,-1));
            return FittingCoefApogee[Index].Coef[0]+
                   FittingCoefApogee[Index].Coef[1]*LogFeH+
                   FittingCoefApogee[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefApogee[Index].Coef[3]*CUBE(LogFeH);
        } else {
            return f;
        }
    } else if(DataType == XFeFeH_GALAH_APOGEE){
        LogFeH = fmin(0.5,fmax(LogFeH,-2));

        double f = FittingCoefGalahApogee[Index].Coef[0]+
                   FittingCoefGalahApogee[Index].Coef[1]*LogFeH+
                   FittingCoefGalahApogee[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefGalahApogee[Index].Coef[3]*CUBE(LogFeH);
        return f;
    } else if(DataType == XFeFeH_GALAH_SAGAFull){
        LogFeH = fmin(0.5,fmax(LogFeH,-2));

        double f = FittingCoefGalahSagaFull[Index].Coef[0]+
                   FittingCoefGalahSagaFull[Index].Coef[1]*LogFeH+
                   FittingCoefGalahSagaFull[Index].Coef[2]*SQ(LogFeH)+
                   FittingCoefGalahSagaFull[Index].Coef[3]*CUBE(LogFeH);
        if(strncmp("N",FittingCoefGalahSagaFull[Index].Name,MaxCharactersInLine-1) == 0){
            return fmax(f,0.0);
        } else if(strncmp("Sr",FittingCoefGalahSagaFull[Index].Name,MaxCharactersInLine-1) == 0){
            return fmax(f,0.0);
        } else {
            return f;
        }

    }
}

/*
 * This is the function which evaluates an adequate value of [X/Fe] of
 * a given element "ElementName" at the given metallicity "FeH".  The dependence
 * of "FeH" is estimated by using the pre-evaluated fitting function. The fitting
 * coefficients are obtained by using the SAGA database (Suda et al. 2008) and
 * APOGEE data. 
 */
double EstimateXFe(const char ElementName[], const double FeH, const int DataType){
    return FitModel(GetIndex(ElementName,0,DataType),FeH,DataType);
}

/*
 * This is the function which evaluates an adequate value of the number density of
 * a given element "ElementName" at the given metallicity "FeH".  The dependence
 * of "FeH" is estimated by using the pre-evaluated fitting function. Not that
 * "FeH" actually means that [Fe/H], the normalized fraction of iron normalized
 * by the solar values. The fitting coefficients are obtained by using the SAGA
 * database (Suda et al. 2008) and APOGEE data. 
 */
double EstimateNX(const char ElementName[], const double FeH, const int DataType){

    int Index = GetIndex(ElementName,0,DataType);

    // Get nB
    //    nB = Z*10.0^(GetSolarAbudance("Fe")-12.0)
    double nB = FeH*pow(10.0,SolarAbundancePattern[GetIndex("Fe",2,DataType)].Value-12.0);

    // Get nAnB_sun
    //nAnB_sun = 10.0^(GetSolarAbudance(Element)-12.0)/10.0^(GetSolarAbudance("Fe")-12.0)
    double nAnB_sun = pow(10.0,(SolarAbundancePattern[GetIndex(ElementName,2,DataType)].Value-12.0))/
                      pow(10.0,(SolarAbundancePattern[GetIndex("Fe",2,DataType)].Value-12.0));

    if(Index == NONE){
        return nB*nAnB_sun;
        //return NONE;
    } else {
        //return nB*nAnB_sun*10.0^(XFe)
        return nB*nAnB_sun*pow(10.0,FitModel(GetIndex(ElementName,0,DataType),FeH,DataType));
    }
}

/*
 * This is the function which evaluates an adequate value of the mass of a given
 * element "ElementName" at the given metallicity "Fe".  The dependence of "Fe" is
 * estimated by using the pre-evaluated fitting function. The fitting
 * coefficients are obtained by using the SAGA database (Suda et al. 2008). 
 */
double EstimateMassX(const char ElementName[], const double FeH, const int DataType){

    int Index = GetIndex(ElementName,0,DataType);
    return AverageElementMass[GetIndex(ElementName,1,DataType)].Value*EstimateNX(ElementName,FeH,DataType);
}


static double MassFractionA09[3] = {0.7381,0.2485,0.0134};
void EstimateAbudancePattern(double Elements[], const double FeH, const int DataType){
    // Assume that this is the standard CELib 13 elements.
    // H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe,Ni,Eu

    char labels[13][MaxCharactersInLine] = {"H","He","C","N","O","Ne","Mg","Si","S","Ca","Fe","Ni","Eu"};
    double FeHCurrent = FeH;

    // FeH
    //[X/Fe] = log10(nX/nFe)-log10(nX/nFe)_sun
    //[X/Fe] = log10({nX/nFe}/{nX/nFe}_sun)
    //10^[X/Fe] * {nX/nFe}_sun = nX/nFe
    // nX = nFe*(nX/nFe)_sun*10^[X/Fe]

    double nXnFe_sun[13] = {0.e0};
    double nX[13] = {0.e0};
    double MassX[13] = {0.e0};
    for(int i=0;i<13;i++){
        if(i>1){
            Elements[i] = EstimateMassX(labels[i],FeHCurrent,DataType);
        } else {
            Elements[i] = AverageElementMass[GetIndex(labels[i],1,DataType)].Value
                *pow(10.0,(SolarAbundancePattern[GetIndex(labels[i],2,DataType)].Value-12.0));
        }
    }

    double Zsum = 0.e0;
    for(int i=2;i<13;i++)
        Zsum += Elements[i];

    double Y = MassFractionA09[1];
    // vprint(Y/Elements[1]);
    double correction = Y/Elements[1];
    // Y = Elements[1]
    // Z:Y = Elements[2>=]:Elemtns[1]
    // Z = Y*Elements[2>=]/Elements[1]
    double Z = Zsum*correction;

    // X+Y+Z=1
    // X'+Y'+Z'=C
    // C(X+Y+Z)=C
    // CX =X'
    // Y'=CY
    // Z'=CZ
    // C=Y'/Y = Elements[1]/Y;
    // vprint(Elements[1]/Y);
    double X = 1-Y-Z;
    // vprint(X);
    // vprint(Y);
    // vprint(Z);
    // vprint(X+Y+Z);
    Elements[0] = 1-Y-Z; 
    Elements[1] = Y;
    for(int i=2;i<13;i++)
        Elements[i] *= correction;

    double sum = 0.e0;
    for(int i=0;i<13;i++)
        sum += Elements[i];
    // vprint(sum);

    return ;
}



void EstimateAbudancePatternForMass(double Elements[], const double FeH, const double Mass, const int DataType){
    // Assume that this is the standard CELib 13 elements.
    // H,He,C,N,O,Ne,Mg,Si,S,Ca,Fe,Ni,Eu

    char labels[13][MaxCharactersInLine] = {"H","He","C","N","O","Ne","Mg","Si","S","Ca","Fe","Ni","Eu"};

    double FeHCurrent = FeH;

    // FeH
    //[X/Fe] = log10(nX/nFe)-log10(nX/nFe)_sun
    //[X/Fe] = log10({nX/nFe}/{nX/nFe}_sun)
    //10^[X/Fe] * {nX/nFe}_sun = nX/nFe
    // nX = nFe*(nX/nFe)_sun*10^[X/Fe]

    double nXnFe_sun[13] = {0.e0};
    double nX[13] = {0.e0};
    double MassX[13] = {0.e0};
    for(int i=0;i<13;i++){
        if(i>1){
            Elements[i] = EstimateMassX(labels[i],FeHCurrent,DataType);
        } else {
            Elements[i] = AverageElementMass[GetIndex(labels[i],1,DataType)].Value
                *pow(10.0,(SolarAbundancePattern[GetIndex(labels[i],2,DataType)].Value-12.0));
        }
    }

    double Zsum = 0.e0;
    for(int i=2;i<13;i++)
        Zsum += Elements[i];

    double Y = MassFractionA09[1];
    double correction = Y/Elements[1];
    double Z = Zsum*correction;
    double X = 1-Y-Z;   

    Elements[0] = 1-Y-Z; 
    Elements[1] = Y;
    for(int i=2;i<13;i++)
        Elements[i] *= correction;

    for(int i=0;i<13;i++)
        Elements[i] *= Mass;

    return ;
}



void EstimateAbundancePatternTest(const int DataType){

    FILE *fp;
    char fname[MaxCharactersInLine];
    Snprintf(fname,"test.%02d.dat",DataType);
    FileOpen(fp,fname,"w");
    
    double logFeHmin = -3;
    double logFeHmax = 0.5;
    int Nbin = 100;

    fprintf(fp,"#[Fe/H] [C/Fe] [N/Fe] [O/Fe] [Ne/Fe] [Mg/Fe] [Si/Fe] [S/Fe] [Ca/Fe] [Ni/Fe] [Eu/Fe]\n");

    double dlogFeH = (logFeHmax-logFeHmin)/Nbin;
    for(int i=0;i<Nbin;i++){
        double logFeH = logFeHmin+i*dlogFeH;          
        double FeH = pow(10.0,logFeH);
        vprint(logFeH);

        double Elements[13] = {0.e0};
        EstimateAbudancePattern(Elements,pow(10.0,logFeH),DataType);

        double nC_sun  = pow(10.0,SolarAbundancePattern[GetIndex("C",2,DataType)].Value-12.0);
        double nN_sun  = pow(10.0,SolarAbundancePattern[GetIndex("N",2,DataType)].Value-12.0);
        double nO_sun  = pow(10.0,SolarAbundancePattern[GetIndex("O",2,DataType)].Value-12.0);
        double nNe_sun = pow(10.0,SolarAbundancePattern[GetIndex("Ne",2,DataType)].Value-12.0);
        double nMg_sun = pow(10.0,SolarAbundancePattern[GetIndex("Mg",2,DataType)].Value-12.0);
        double nSi_sun = pow(10.0,SolarAbundancePattern[GetIndex("Si",2,DataType)].Value-12.0);
        double nS_sun  = pow(10.0,SolarAbundancePattern[GetIndex("S",2,DataType)].Value-12.0);
        double nCa_sun = pow(10.0,SolarAbundancePattern[GetIndex("Ca",2,DataType)].Value-12.0);
        double nNi_sun = pow(10.0,SolarAbundancePattern[GetIndex("Ni",2,DataType)].Value-12.0);
        double nEu_sun = pow(10.0,SolarAbundancePattern[GetIndex("Eu",2,DataType)].Value-12.0);

        double nFe_sun = pow(10.0,SolarAbundancePattern[GetIndex("Fe",2,DataType)].Value-12.0);
        double nH_sun = 1.0;
        double nFe = Elements[10]/AverageElementMass[GetIndex("Fe",1,DataType)].Value;

        double nH = Elements[0]/AverageElementMass[GetIndex("H",1,DataType)].Value;
        // vprint(nH);
        // vprint(Elements[0]);
        // vprint(nFe_sun);

#if 1
        fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g\n",
                log10((Elements[10]/AverageElementMass[GetIndex("Fe",1,DataType)].Value)/nH)  -log10(nFe_sun/nH_sun),
                log10((Elements[2]/AverageElementMass[GetIndex("C",1,DataType)].Value)/nFe)   -log10(nC_sun/nFe_sun),
                log10((Elements[3]/AverageElementMass[GetIndex("N",1,DataType)].Value)/nFe)   -log10(nN_sun/nFe_sun),
                log10((Elements[4]/AverageElementMass[GetIndex("O",1,DataType)].Value)/nFe)   -log10(nO_sun/nFe_sun),
                log10((Elements[5]/AverageElementMass[GetIndex("Ne",1,DataType)].Value)/nFe)  -log10(nNe_sun/nFe_sun),
                log10((Elements[6]/AverageElementMass[GetIndex("Mg",1,DataType)].Value)/nFe)  -log10(nMg_sun/nFe_sun),
                log10((Elements[7]/AverageElementMass[GetIndex("Si",1,DataType)].Value)/nFe)  -log10(nSi_sun/nFe_sun),
                log10((Elements[8]/AverageElementMass[GetIndex("S",1,DataType)].Value)/nFe)   -log10(nS_sun/nFe_sun),
                log10((Elements[9]/AverageElementMass[GetIndex("Ca",1,DataType)].Value)/nFe)  -log10(nCa_sun/nFe_sun),
                log10((Elements[11]/AverageElementMass[GetIndex("Ni",1,DataType)].Value)/nFe) -log10(nNi_sun/nFe_sun),
                log10((Elements[12]/AverageElementMass[GetIndex("Eu",1,DataType)].Value)/nFe) -log10(nEu_sun/nFe_sun)
                );
#else 
        fprintf(fp,"%g %g %g %g %g %g %g %g\n",
                Elements[10],AverageElementMass[GetIndex("Fe",1,DataType)].Value,nH,nFe_sun,nH_sun,
                log10((Elements[10]/AverageElementMass[GetIndex("Fe",1,DataType)].Value)/nH),log10(nFe_sun/nH_sun),
                log10((Elements[10]/AverageElementMass[GetIndex("Fe",1,DataType)].Value)/nH)-log10(nFe_sun/nH_sun));
#endif
    }
    fclose(fp);

    return ;
}



