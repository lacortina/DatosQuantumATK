import numpy as np
import matplotlib.pyplot as plt
import re

# =========================
# 1. LECTURA DEL FICHERO
# =========================
def read_band_file(filename):
    bands = []
    current_band = []

    fermi_level = None

    with open(filename, 'r') as f:
        lines = f.readlines()

    for line in lines:
        # =========================
        # 1. Leer Fermi level (robusto)
        # =========================
        if "Fermi level" in line:
            match = re.search(r'Fermi level\s*=\s*([-+0-9.eE]+)', line)
            if match:
                fermi_level = float(match.group(1))

        # =========================
        # 2. Detectar nueva banda
        # =========================
        if line.strip().startswith("# Band"):
            if current_band:
                bands.append(np.array(current_band))
                current_band = []
            continue

        # =========================
        # 3. Leer datos numéricos
        # =========================
        parts = line.split()
        if len(parts) >= 2:
            try:
                # Solo necesitamos la energía
                energy = float(parts[1])
                current_band.append(energy)
            except ValueError:
                continue

    # Añadir última banda
    if current_band:
        bands.append(np.array(current_band))

    # =========================
    # 4. Validación
    # =========================
    if fermi_level is None:
        raise ValueError(f"No se encontró el nivel de Fermi en {filename}")

    if len(bands) == 0:
        raise ValueError(f"No se encontraron bandas en {filename}")

    return bands, fermi_level



# =========================
# 2. FILTRADO DE BANDAS
# =========================
def filter_bands(bands, emin=-15, emax=15):
    filtered = []

    for band in bands:
        e0 = band[0]

        if e0 < emin:
            continue
        if e0 > emax:
            continue

        filtered.append(band)

    return filtered


# =========================
# 3. RENORMALIZACIÓN
# =========================
def renormalize_bands(bands, fermi_file, fermi_global):
    shift = fermi_file - fermi_global
    return [band + shift for band in bands]


# =========================
# 4. PLOT
# =========================
def plot_bands(all_band_sets, styles=None, labels=None, save = True, filename = 'ComLAmpliado', dpi=300):
    """
    all_band_sets: lista de listas de bandas
    styles: lista de dicts -> {'color': 'k', 'linestyle': '-'}
    """

    plt.figure(figsize=(6, 5))

    n_k = len(all_band_sets[0][0])
    k_axis = np.linspace(0, 1, n_k)

    # Puntos de alta simetría
    ticks = [0, 20/(n_k-1), 40/(n_k-1), 60/(n_k-1)]
    tick_labels = [r'$\Gamma$', 'M', 'K', r'$\Gamma$']

    # Dibujar líneas verticales
    for t in ticks:
        plt.axvline(t, color='gray', linewidth=0.5)

    # Dibujar nivel de Fermi global
    plt.axhline(0, color='red', linestyle='--', linewidth=1)

    # Plot de bandas
    for i, bands in enumerate(all_band_sets):
        style = styles[i] if styles else {}

        for band in bands:
            plt.plot(k_axis, band,
                     color=style.get('color', 'black'),
                     linestyle=style.get('linestyle', '-'),
                     linewidth=1)

    plt.xticks(ticks, tick_labels, fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel("Energy (eV)", fontsize = 14)
    #plt.xlim(0, 1)
    plt.ylim(-12, 12)

    if labels:
        for i, label in enumerate(labels):
            plt.plot([], [], label=label,
                     color=styles[i].get('color', 'black'),
                     linestyle=styles[i].get('linestyle', '-'))
        plt.legend(fontsize=12)

    plt.tight_layout()
    if save:
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"Figura guardada en: {filename}")

    
    plt.show()

# =========================
# 5. PIPELINE COMPLETO
# =========================
def process_file(filename, fermi_global):
    bands, fermi_file = read_band_file(filename)
    bands = filter_bands(bands)
    bands = renormalize_bands(bands, fermi_file, fermi_global)
    return bands

# =========================
# EJEMPLO DE USO
# =========================

fermi_global = -5.04949  # <-- tú decides este valor

files = [
    "./Datos/Bandas9.txt"
    #"./Datos/Bandas32.txt",
    #"./Datos/Bandas21.txt",
    #"./Datos/Bandas13.txt",
    #"./Datos/Bandas9.txt",
     
]

styles = [
    {'color': 'black', 'linestyle': '-'},
    #{'color': 'red', 'linestyle': '--'},
    #{'color': 'blue', 'linestyle': '--'},
    #{'color': 'purple', 'linestyle': '--'},
    #{'color': 'brown', 'linestyle': '--'}
    
]

labels = ["9"]#, "32", "21", "13", "9"]

all_band_sets = []

for f in files:
    bands = process_file(f, fermi_global)
    all_band_sets.append(bands)

plot_bands(all_band_sets, styles=styles, labels=labels)