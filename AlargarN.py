import numpy as np

# --------------------------------------------------
# 1. Leer CIF
# --------------------------------------------------
def read_cif(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    cell = {}
    atoms = []

    read_atoms = False

    for line in lines:
        if "_cell_length_a" in line:
            cell['a'] = float(line.split()[1])
        elif "_cell_length_b" in line:
            cell['b'] = float(line.split()[1])
        elif "_cell_length_c" in line:
            cell['c'] = float(line.split()[1])

        elif "_atom_site_label" in line:
            read_atoms = True
            continue

        elif read_atoms:
            if line.strip() == "" or line.startswith("_"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                atoms.append({
                    "label": parts[0],
                    "fract": np.array([float(parts[1]),
                                       float(parts[2]),
                                       float(parts[3])])
                })

    return cell, atoms, lines

# --------------------------------------------------
# 2. Fraccional → real
# --------------------------------------------------
def frac_to_cart(frac, cell):
    return np.array([
        frac[0] * cell['a'],
        frac[1] * cell['b'],
        frac[2] * cell['c']
    ])


def cart_to_frac(cart, cell):
    return np.array([
        cart[0] / cell['a'],
        cart[1] / cell['b'],
        cart[2] / cell['c']
    ])

# --------------------------------------------------
# 3. Crear set por rango en X
# --------------------------------------------------
def create_set(atoms, cell, xmin, xmax):
    selected = []

    for atom in atoms:
        r = frac_to_cart(atom["fract"], cell)
        if xmin <= r[0] <= xmax:
            selected.append({
                "label": atom["label"],
                "cart": r,
                "atom_ref": atom
            })

    return selected

# --------------------------------------------------
# 4. Emparejar y calcular diferencias en X
# --------------------------------------------------
def match_and_compute(set1, set2, tol=5e-2):
    diffs = []

    used = set()

    for i, a1 in enumerate(set1):
        y1, z1 = a1["cart"][1], a1["cart"][2]

        for j, a2 in enumerate(set2):
            if j in used:
                continue

            y2, z2 = a2["cart"][1], a2["cart"][2]

            if abs(y1 - y2) < tol and abs(z1 - z2) < tol:
                dx = a1["cart"][0] - a2["cart"][0]
                diffs.append(dx)
                used.add(j)
                break

    diffs = np.array(diffs)

    mean = np.mean(diffs)
    var = np.var(diffs)

    return diffs, mean, var


# --------------------------------------------------
# 5. Desplazar set en X
# --------------------------------------------------
def shift_set(set_atoms, shift_x):
    for atom in set_atoms:
        atom["cart"][0] += shift_x



def build_new_atoms(set3, set4, cell_new):
    """
    Combina set3 + set4 (ya desplazado) y recalcula coordenadas fraccionales
    """
    new_atoms = []

    for atom in set3 + set4:
        cart = atom["cart"]

        fract = np.array([
            cart[0] / cell_new['a'],
            cart[1] / cell_new['b'],
            cart[2] / cell_new['c']
        ])

        new_atoms.append({
            "label": atom["label"],  # Mo o S
            "fract": fract
        })

    return new_atoms


def write_new_cif(output_file, cell, new_atoms):
    with open(output_file, "w") as f:

        f.write("data_global\n")
        f.write(f"_cell_length_a {cell['a']}\n")
        f.write(f"_cell_length_b {cell['b']}\n")
        f.write(f"_cell_length_c {cell['c']}\n")
        f.write("_cell_angle_alpha 90\n")
        f.write("_cell_angle_beta 90\n")
        f.write("_cell_angle_gamma 90\n")
        f.write("_symmetry_space_group_name_H-M 'P -1'\n")

        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")

        for atom in new_atoms:
            f.write(f"{atom['label']} "
                    f"{atom['fract'][0]:.6f} "
                    f"{atom['fract'][1]:.6f} "
                    f"{atom['fract'][2]:.6f}\n")

# --------------------------------------------------
# 7. PIPELINE PRINCIPAL
# --------------------------------------------------

def main():
    N = 1  # número de iteraciones

    cell, atoms, lines = read_cif("MoS2_gb_60_zz_bulk_opt.cif_cut.cif")

    # -------------------------
    # Paso 1: calcular mean inicial
    # -------------------------
    set1 = create_set(atoms, cell, xmin=6.5, xmax=9.5)
    set2 = create_set(atoms, cell, xmin=9.5, xmax=12.5)

    diffs, mean, var = match_and_compute(set2, set1)

    print("Media ΔX:", mean)
    print("Varianza ΔX:", var)

    # -------------------------
    # Inicialización
    # -------------------------
    current_atoms = atoms.copy()
    cell_current = cell.copy()

    # -------------------------
    # Iteraciones
    # -------------------------
    for n in range(1, N+1):

        print(f"\nIteración {n}")

        shift_total = n * mean

        # límites dinámicos
        xmin1 = 0
        xmax1 = 9.5 + (n-1)*mean

        xmin2 = 6.5 + (n-1)*mean
        xmax2 = 100

        xmin3 = 0
        xmax3 = 29 + (n-1)*mean

        xmin4 = 25.8 + (n-1)*mean
        xmax4 = 100

        # crear sets
        setA = create_set(current_atoms, cell_current, xmin=xmin1, xmax=xmax1)
        setB = create_set(current_atoms, cell_current, xmin=xmin2, xmax=xmax2)

        # desplazar setB SOLO por mean (no acumulado aquí)
        shift_set(setB, mean)

        # actualizar celda
        cell_current['a'] += mean

        # construir nueva lista de átomos
        new_atoms = build_new_atoms(setA, setB, cell_current)

        # actualizar para siguiente iteración
        current_atoms = new_atoms

        setC = create_set(current_atoms, cell_current, xmin=xmin3, xmax=xmax3)
        setD = create_set(current_atoms, cell_current, xmin=xmin4, xmax=xmax4) 

        shift_set(setD, mean)

        cell_current['a'] += mean

        new_atoms = build_new_atoms(setC, setD, cell_current)
        

        current_atoms = new_atoms   
    # -------------------------
    # Guardar resultado final
    # -------------------------
    write_new_cif("output.cif", cell_current, current_atoms)

if __name__ == "__main__":
    main()
