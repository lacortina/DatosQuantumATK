def cut_cif_by_y_window(input_path, output_path, y_min, y_max):
    """
    Filtra átomos en un CIF manteniendo solo aquellos con:
        y_min < fract_y < y_max

    input_path: ruta del CIF original
    output_path: ruta del CIF filtrado
    y_min, y_max: umbrales en coordenadas fraccionales (0-1)
    """

    assert 0.0 <= y_min < y_max <= 1.0, "Los límites deben cumplir 0 <= y_min < y_max <= 1"
    
    with open(input_path, 'r') as f:
        lines = f.readlines()

    header = []
    atom_lines = []
    in_atom_loop = False
    atom_loop_started = False

    for line in lines:
        stripped = line.strip()

        # Detectar inicio de loop
        if stripped.startswith("loop_"):
            in_atom_loop = True
            atom_loop_started = False
            header.append(line)
            continue

        # Columnas de átomos
        if in_atom_loop and stripped.startswith("_atom_site"):
            atom_loop_started = True
            header.append(line)
            continue

        # Datos de átomos
        if in_atom_loop and atom_loop_started:
            if stripped == "" or stripped.startswith("_") or stripped.startswith("loop_"):
                in_atom_loop = False
                header.append(line)
            else:
                atom_lines.append(line)
            continue

        header.append(line)

    # 🔹 Filtrado
    filtered_atoms = []
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 4:
            continue

        y = float(parts[2])  # fract_y

        if y_min < y < y_max:
            filtered_atoms.append(line)

    # 🔹 Escritura
    with open(output_path, 'w') as f:
        for line in header:
            f.write(line)
        for line in filtered_atoms:
            f.write(line)

    print(f"Archivo guardado en: {output_path}")
    print(f"Átomos originales: {len(atom_lines)}")
    print(f"Átomos filtrados: {len(filtered_atoms)}")

cut_cif_by_y_window("MoS2_gb_60_zz_bulk_opt.cif", "MoS2_gb_60_zz_bulk_opt.cif_cut.cif", y_min=0.21, y_max = 0.55)
