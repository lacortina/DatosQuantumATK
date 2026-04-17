import numpy as np

def filter_cif_by_z(input_cif, output_cif, zmin, zmax):
    """
    Filtra átomos de un CIF según coordenada fraccionaria z

    Parámetros:
    - input_cif: archivo de entrada
    - output_cif: archivo de salida
    - zmin, zmax: rango en coordenadas fraccionarias
    """

    with open(input_cif, 'r') as f:
        lines = f.readlines()

    new_lines = []
    atoms = []

    read_atoms = False
    atom_start_idx = None

    # -------------------------
    # Leer CIF
    # -------------------------
    for i, line in enumerate(lines):
        if "_atom_site_label" in line:
            read_atoms = True
            atom_start_idx = i + 1
            new_lines.extend(lines[:i+1])
            continue

        if read_atoms:
            if line.strip() == "" or line.startswith("_") or line.startswith("loop_"):
                continue

            parts = line.split()
            if len(parts) >= 4:
                label = parts[0]
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])

                atoms.append((label, x, y, z))

    # -------------------------
    # Filtrar átomos
    # -------------------------
    filtered_atoms = []

    for atom in atoms:
        label, x, y, z = atom

        # aplicar PBC opcional
        z_mod = z % 1.0

        if zmin <= z_mod <= zmax:
            filtered_atoms.append((label, x, y, z_mod))

    # -------------------------
    # Escribir nuevo CIF
    # -------------------------
    with open(output_cif, 'w') as f:

        # escribir cabecera original (sin átomos)
        for line in lines:
            if "_atom_site_label" in line:
                f.write(line)
                f.write("_atom_site_fract_x\n")
                f.write("_atom_site_fract_y\n")
                f.write("_atom_site_fract_z\n")
                break
            f.write(line)

        # escribir átomos filtrados
        for atom in filtered_atoms:
            label, x, y, z = atom
            f.write(f"{label} {x:.6f} {y:.6f} {z:.6f}\n")


filter_cif_by_z(
    "MoS2_gb_60L4_zz_bulk.cif",
    "MoS2_gb_60L4_zz_bulkMono.cif",
    zmin=0,
    zmax=0.6
)
