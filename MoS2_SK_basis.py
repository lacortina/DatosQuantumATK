# =============================================================================
# MoS2_SK_basis.py
# Slater-Koster BasisSet para bicapa MoS2
# Parametros de: Zahid et al., "A generic tight-binding model for monolayer,
# bilayer and bulk MoS2"
# Base: s, p, d para Mo y S
# =============================================================================

from atkpython import *

# =============================================================================
# PSEUDOELEMENTO S1 (azufres de la capa 2, mismas propiedades que S)
# =============================================================================
S1 = PseudoElement(
    symbol        = 'S1',
    atomic_number = 16,
    mass          = 32.06 * amu,
    color         = [1.0, 0.8, 0.0],
)

# =============================================================================
# BASIS SETS (s, p, d para Mo y S/S1)
# =============================================================================
Mo_basis = SlaterKosterBasisSet(
    element  = Molybdenum,
    orbitals = [
        SphericalHarmonic(n=5, l=0, m= 0),   # s
        SphericalHarmonic(n=5, l=1, m=-1),   # px
        SphericalHarmonic(n=5, l=1, m= 0),   # pz
        SphericalHarmonic(n=5, l=1, m= 1),   # py
        SphericalHarmonic(n=4, l=2, m=-2),   # d_{x2-y2}
        SphericalHarmonic(n=4, l=2, m=-1),   # d_{xz}
        SphericalHarmonic(n=4, l=2, m= 0),   # d_{z2}
        SphericalHarmonic(n=4, l=2, m= 1),   # d_{yz}
        SphericalHarmonic(n=4, l=2, m= 2),   # d_{xy}
    ]
)

S_basis = SlaterKosterBasisSet(
    element  = Sulfur,
    orbitals = [
        SphericalHarmonic(n=3, l=0, m= 0),   # s
        SphericalHarmonic(n=3, l=1, m=-1),   # px
        SphericalHarmonic(n=3, l=1, m= 0),   # pz
        SphericalHarmonic(n=3, l=1, m= 1),   # py
        SphericalHarmonic(n=3, l=2, m=-2),   # d_{x2-y2}
        SphericalHarmonic(n=3, l=2, m=-1),   # d_{xz}
        SphericalHarmonic(n=3, l=2, m= 0),   # d_{z2}
        SphericalHarmonic(n=3, l=2, m= 1),   # d_{yz}
        SphericalHarmonic(n=3, l=2, m= 2),   # d_{xy}
    ]
)

S1_basis = SlaterKosterBasisSet(
    element  = S1,
    orbitals = [
        SphericalHarmonic(n=3, l=0, m= 0),
        SphericalHarmonic(n=3, l=1, m=-1),
        SphericalHarmonic(n=3, l=1, m= 0),
        SphericalHarmonic(n=3, l=1, m= 1),
        SphericalHarmonic(n=3, l=2, m=-2),
        SphericalHarmonic(n=3, l=2, m=-1),
        SphericalHarmonic(n=3, l=2, m= 0),
        SphericalHarmonic(n=3, l=2, m= 1),
        SphericalHarmonic(n=3, l=2, m= 2),
    ]
)

# =============================================================================
# PARAMETROS ONSITE (eV)
# =============================================================================
Mo_onsite = OnSiteParameters(
    element = Molybdenum,
    s       =  5.5994 * eV,
    p       =  6.7128 * eV,
    d       =  2.6429 * eV,
)

S_onsite = OnSiteParameters(
    element = Sulfur,
    s       =  7.6595 * eV,
    p       = -2.1537 * eV,
    d       =  8.7689 * eV,
)

S1_onsite = OnSiteParameters(
    element = S1,
    s       =  7.6595 * eV,
    p       = -2.1537 * eV,
    d       =  8.7689 * eV,
)

# =============================================================================
# SPIN-ORBIT COUPLING (eV)
# =============================================================================
Mo_soc = SpinOrbitCouplingParameters(
    element   = Molybdenum,
    p_orbital = 1.0675 * eV,
    d_orbital = 1.0675 * eV,
)

S_soc = SpinOrbitCouplingParameters(
    element   = Sulfur,
    p_orbital = 0.2129 * eV,
    d_orbital = 0.2129 * eV,
)

S1_soc = SpinOrbitCouplingParameters(
    element   = S1,
    p_orbital = 0.2129 * eV,
    d_orbital = 0.2129 * eV,
)

# =============================================================================
# DISTANCIAS Y CUTOFFS (Angstrom)
# Los cutoffs se definen como el punto medio entre la distancia del vecino
# considerado y el siguiente vecino no deseado, con un margen de seguridad
# =============================================================================
d_MoMo      = 3.1604   # distancia Mo-Mo primeros vecinos
d_MoS       = 2.41764  # distancia Mo-S  primeros vecinos
d_SS_intra  = 3.1604   # distancia S-S / S1-S1 intracapa
d_SS_inter  = 3.4903   # distancia S-S1 intercapa

# Cutoffs: distancia_vecino + margen (~15% de la distancia)
cutoff_MoMo     = d_MoMo     * 1.15   # ~ 3.63 Ang
cutoff_MoS      = d_MoS      * 1.15   # ~ 2.78 Ang
cutoff_SS_intra = d_SS_intra * 1.15   # ~ 3.63 Ang
cutoff_SS_inter = d_SS_inter * 1.15   # ~ 4.01 Ang

# =============================================================================
# TABLAS SK: MO-MO
# =============================================================================
MoMo_H_table = NumericalSlaterKosterTable(
    distances = [d_MoMo, cutoff_MoMo] * Angstrom,
    ss_sigma  = [ 0.1768,  0.0],
    sp_sigma  = [ 1.0910,  0.0],
    pp_sigma  = [-0.3842,  0.0],
    pp_pi     = [ 0.5203,  0.0],
    sd_sigma  = [-0.5635,  0.0],
    pd_sigma  = [-0.2316,  0.0],
    pd_pi     = [ 0.0582,  0.0],
    dd_sigma  = [ 0.3602,  0.0],
    dd_pi     = [ 0.0432,  0.0],
    dd_delta  = [ 0.1008,  0.0],
)

MoMo_S_table = NumericalSlaterKosterTable(
    distances = [d_MoMo, cutoff_MoMo] * Angstrom,
    ss_sigma  = [-0.0575,  0.0],
    sp_sigma  = [ 0.0057,  0.0],
    pp_sigma  = [ 0.0296,  0.0],
    pp_pi     = [ 0.0946,  0.0],
    sd_sigma  = [-0.1082,  0.0],
    pd_sigma  = [ 0.0212,  0.0],
    pd_pi     = [-0.0448,  0.0],
    dd_sigma  = [-0.0216,  0.0],
    dd_pi     = [-0.0285,  0.0],
    dd_delta  = [ 0.0432,  0.0],
)

MoMo_SK = SlaterKosterTablePair(
    element_pair      = (Molybdenum, Molybdenum),
    hamiltonian_table = MoMo_H_table,
    overlap_table     = MoMo_S_table,
)

# =============================================================================
# TABLAS SK: MO-S
# Notacion del paper: Vsss = s(Mo)-s(S), Vsps = s(Mo)-p(S),
#                     Vpss = p(Mo)-s(S), Vsds = s(Mo)-d(S), etc.
# =============================================================================
MoS_H_table = NumericalSlaterKosterTable(
    distances = [d_MoS, cutoff_MoS] * Angstrom,
    ss_sigma  = [-0.0917,  0.0],   # Vsss
    sp_sigma  = [-1.6515,  0.0],   # Vsps  s(Mo)-p(S)
    ps_sigma  = [ 0.6656,  0.0],   # Vpss  p(Mo)-s(S)
    pp_sigma  = [ 1.4008,  0.0],   # Vpps
    pp_pi     = [-0.4812,  0.0],   # Vppp
    sd_sigma  = [-1.0654,  0.0],   # Vsds  s(Mo)-d(S)
    ds_sigma  = [ 0.2177,  0.0],   # Vdss  d(Mo)-s(S)
    pd_sigma  = [ 2.1898,  0.0],   # Vpds  p(Mo)-d(S)
    dp_sigma  = [-2.9732,  0.0],   # Vdps  d(Mo)-p(S)
    pd_pi     = [-1.9408,  0.0],   # Vpdp
    dp_pi     = [ 0.7739,  0.0],   # Vdpp
    dd_sigma  = [-3.1425,  0.0],   # Vdds
    dd_pi     = [ 2.4975,  0.0],   # Vddp
    dd_delta  = [-0.3703,  0.0],   # Vddd
)

MoS_S_table = NumericalSlaterKosterTable(
    distances = [d_MoS, cutoff_MoS] * Angstrom,
    ss_sigma  = [ 0.0294,  0.0],
    sp_sigma  = [ 0.1765,  0.0],   # Ssps
    ps_sigma  = [ 0.1042,  0.0],   # Spss
    pp_sigma  = [-0.1865,  0.0],
    pp_pi     = [ 0.0303,  0.0],
    sd_sigma  = [-0.1432,  0.0],   # Ssds
    ds_sigma  = [-0.0480,  0.0],   # Sdss
    pd_sigma  = [ 0.2002,  0.0],   # Spds
    dp_sigma  = [ 0.0942,  0.0],   # Sdps
    pd_pi     = [-0.2435,  0.0],   # Spdp
    dp_pi     = [ 0.0132,  0.0],   # Sdpp
    dd_sigma  = [ 0.0273,  0.0],
    dd_pi     = [ 0.1940,  0.0],
    dd_delta  = [ 0.1261,  0.0],
)

MoS_SK = SlaterKosterTablePair(
    element_pair      = (Molybdenum, Sulfur),
    hamiltonian_table = MoS_H_table,
    overlap_table     = MoS_S_table,
)

# Mo-S1: mismos parametros que Mo-S
MoS1_SK = SlaterKosterTablePair(
    element_pair      = (Molybdenum, S1),
    hamiltonian_table = MoS_H_table,   # reutilizamos la misma tabla
    overlap_table     = MoS_S_table,
)

# =============================================================================
# TABLAS SK: S-S INTRACAPA (y S1-S1, mismos parametros)
# =============================================================================
SS_intra_H_table = NumericalSlaterKosterTable(
    distances = [d_SS_intra, cutoff_SS_intra] * Angstrom,
    ss_sigma  = [ 0.3093,  0.0],
    sp_sigma  = [-0.9210,  0.0],
    pp_sigma  = [ 0.7132,  0.0],
    pp_pi     = [-0.1920,  0.0],
    sd_sigma  = [-0.2016,  0.0],
    pd_sigma  = [-0.5204,  0.0],
    pd_pi     = [-0.1203,  0.0],
    dd_sigma  = [ 0.8347,  0.0],
    dd_pi     = [ 0.7434,  0.0],
    dd_delta  = [-0.1919,  0.0],
)

SS_intra_S_table = NumericalSlaterKosterTable(
    distances = [d_SS_intra, cutoff_SS_intra] * Angstrom,
    ss_sigma  = [-0.0532,  0.0],
    sp_sigma  = [ 0.0240,  0.0],
    pp_sigma  = [ 0.0478,  0.0],
    pp_pi     = [-0.0104,  0.0],
    sd_sigma  = [ 0.0946,  0.0],
    pd_sigma  = [ 0.0724,  0.0],
    pd_pi     = [ 0.0772,  0.0],
    dd_sigma  = [ 0.1849,  0.0],
    dd_pi     = [-0.0429,  0.0],
    dd_delta  = [-0.0333,  0.0],
)

SS_intra_SK = SlaterKosterTablePair(
    element_pair      = (Sulfur, Sulfur),
    hamiltonian_table = SS_intra_H_table,
    overlap_table     = SS_intra_S_table,
)

S1S1_intra_SK = SlaterKosterTablePair(
    element_pair      = (S1, S1),
    hamiltonian_table = SS_intra_H_table,   # reutilizamos
    overlap_table     = SS_intra_S_table,
)

# =============================================================================
# TABLAS SK: S-S1 INTERCAPA
# =============================================================================
SS1_inter_H_table = NumericalSlaterKosterTable(
    distances = [d_SS_inter, cutoff_SS_inter] * Angstrom,
    ss_sigma  = [ 0.3207,  0.0],
    sp_sigma  = [-0.1302,  0.0],
    pp_sigma  = [ 0.7053,  0.0],
    pp_pi     = [-0.0980,  0.0],
    sd_sigma  = [ 0.1164,  0.0],
    pd_sigma  = [-0.0334,  0.0],
    pd_pi     = [-0.0370,  0.0],
    dd_sigma  = [-0.2300,  0.0],
    dd_pi     = [ 0.0050,  0.0],
    dd_delta  = [-0.1104,  0.0],
)

SS1_inter_S_table = NumericalSlaterKosterTable(
    distances = [d_SS_inter, cutoff_SS_inter] * Angstrom,
    ss_sigma  = [-0.1430,  0.0],
    sp_sigma  = [ 0.0196,  0.0],
    pp_sigma  = [-0.0486,  0.0],
    pp_pi     = [ 0.0117,  0.0],
    sd_sigma  = [ 0.0297,  0.0],
    pd_sigma  = [-0.0087,  0.0],
    pd_pi     = [-0.0031,  0.0],
    dd_sigma  = [ 0.0060,  0.0],
    dd_pi     = [-0.0378,  0.0],
    dd_delta  = [ 0.0007,  0.0],
)

SS1_inter_SK = SlaterKosterTablePair(
    element_pair      = (Sulfur, S1),
    hamiltonian_table = SS1_inter_H_table,
    overlap_table     = SS1_inter_S_table,
)

# =============================================================================
# ENSAMBLAJE FINAL DEL BASISSET
# =============================================================================
MoS2_SK_BasisSet = SlaterKosterHamiltonianParameterSet(
    basis_sets        = [Mo_basis, S_basis, S1_basis],
    onsite_parameters = [Mo_onsite, S_onsite, S1_onsite],
    soc_parameters    = [Mo_soc, S_soc, S1_soc],
    pair_parameters   = [
        MoMo_SK,
        MoS_SK,
        MoS1_SK,
        SS_intra_SK,
        S1S1_intra_SK,
        SS1_inter_SK,
    ],
)
