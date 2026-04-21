# =============================================================================
# MoS2_bilayer_bands.py
# Estructura de bandas para bicapa MoS2 con modelo SK de Zahid et al.
# =============================================================================

from atkpython import *

# Importar el BasisSet
exec(open('MoS2_SK_basis.py').read())

# =============================================================================
# PARAMETROS DE RED
# Bicapa MoS2 apilamiento AA' (el mas comun)
# a = 3.1604 Ang (igual que distancia Mo-Mo primeros vecinos)
# c = suficientemente grande para incluir las dos capas + vacio
# =============================================================================
a     = 3.1604   # Angstrom, parametro de red hexagonal
dz_S  = 1.5690   # Angstrom, desplazamiento vertical S respecto a Mo
                 # (distancia Mo-S proyectada en z, aprox sqrt(d_MoS^2 - (a/sqrt(3))^2))
gap   = 3.4903   # Angstrom, separacion S-S1 intercapa (= d_SS_inter)
c     = 2 * (2 * dz_S) + gap + 20.0  # c grande para evitar interaccion entre imagenes

# Posiciones fractionales (z en unidades de c):
# Capa 1:   S(arriba) - Mo - S(abajo)
# gap van der Waals
# Capa 2:   S1(arriba) - Mo - S1(abajo)

# Definimos z=0 en el centro de la bicapa
z_Mo1  =  (gap/2 + dz_S) / c + 0.5
z_S1u  =  (gap/2 + 2*dz_S) / c + 0.5   # S superior capa 1
z_S1d  =   gap/2 / c + 0.5              # S inferior capa 1  (el que acopla con S1)

z_Mo2  = -(gap/2 + dz_S) / c + 0.5
z_S2u  = -gap/2 / c + 0.5               # S1 superior capa 2 (el que acopla con S)
z_S2d  = -(gap/2 + 2*dz_S) / c + 0.5   # S1 inferior capa 2

lattice = HexagonalLattice(
    a = a * Angstrom,
    c = c * Angstrom,
)

# Posiciones atomicas
# Mo en 1/3, 2/3 y S en 2/3, 1/3 (red hexagonal, 2 atomos por celda en plano)
bulk_config = BulkConfiguration(
    bravais_lattice = lattice,
    elements        = [Molybdenum, Sulfur,  Sulfur,    # capa 1: Mo, S_up, S_down
                       Molybdenum, S1,      S1],        # capa 2: Mo, S1_up, S1_down
    fractional_coordinates = [
        [1/3, 2/3, z_Mo1],   # Mo capa 1
        [2/3, 1/3, z_S1u],   # S  superior capa 1
        [2/3, 1/3, z_S1d],   # S  inferior capa 1
        [2/3, 1/3, z_Mo2],   # Mo capa 2  (apilamiento AA': mismo x,y que S)
        [1/3, 2/3, z_S2u],   # S1 superior capa 2
        [1/3, 2/3, z_S2d],   # S1 inferior capa 2
    ],
)

# =============================================================================
# CALCULADOR TB
# =============================================================================
calculator = TightBindingCalculator(
    hamiltonian_parameters = MoS2_SK_BasisSet,
    spin                   = SpinOrbit,
    k_point_sampling       = MonkhorstPackGrid(15, 15, 1),
    occupations_method     = FermiDirac(300 * Kelvin),
)

bulk_config.setCalculator(calculator)
bulk_config.update()

# =============================================================================
# ESTRUCTURA DE BANDAS
# Ruta: Gamma -> K -> M -> Gamma (puntos de alta simetria hexagonal)
# =============================================================================
bandstructure = BandStructure(
    configuration      = bulk_config,
    route              = [
        ('G', [0.000, 0.000, 0.000]),
        ('K', [1/3,  1/3,   0.000]),
        ('M', [0.500, 0.000, 0.000]),
        ('G', [0.000, 0.000, 0.000]),
    ],
    points_per_segment = 60,
)

bandstructure.evaluate()
nlsave('MoS2_bilayer_TB.hdf5', bandstructure)
nlprint(bandstructure)

# =============================================================================
# DOS (opcional, comentar si no se necesita)
# =============================================================================
dos = DensityOfStates(
    configuration  = bulk_config,
    energies       = numpy.linspace(-10, 10, 1000) * eV,
    broadening     = 0.05 * eV,
)

dos.evaluate()
nlsave('MoS2_bilayer_TB.hdf5', dos)
