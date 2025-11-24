import csv


def main():
    """
    Read a CSV of spheres and print MCNP cell and surface cards.

    Expected CSV header:
        type,x_cm,y_cm,z_cm,radius_cm

    Example row:
        Be,-0.17320508,-0.17320508,-0.17320508,0.05
    """

    cell_map     = {"Li4SiO4":  3,    "Be":  4, }
    density_map  = {"Li4SiO4": -2.17, "Be": -1.80, }
    surface_map  = {"Li4SiO4": 13,    "Be": 12, }


    cell_lines = []
    surf_lines = []
    
    CSV_PATH = "sphere_packing.csv" # _vol01

    with open(CSV_PATH, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        rows.reverse()

    # ---- Cell cards ----
    cell_lines.append(f"c       -----------------------------------------------------------------------")
    cell_lines.append(f"c       Cell cards generated from '{CSV_PATH}'")
    cell_lines.append(f"c       -----------------------------------------------------------------------")
    cell_lines.append(f" 3      300  -2.17  -93                                        imp:n=1  tmp=7.76e-8     u=30  $ Li4SiO4 univ vol=5.2360e-4")
    cell_lines.append(f" 4      400  -1.80  -92                                        imp:n=1  tmp=7.76e-8     u=40  $ Be univ      vol=3.3510e-4")


    current, complements = None, []
    for _, row in enumerate(rows):
        type = row["type"].strip()
        if type in ['Li4SiO4', 'Be']:
            cell = cell_map[type]
            surf = surface_map[type]

            if type != current: idx = 0

            cell_id = cell*1000 + 1 + idx
            # surf_id = surf_start_id + idx

            x = float(row["x_cm"])
            y = float(row["y_cm"])
            z = float(row["z_cm"])
            r = float(row["radius_cm"])

            line = f"{cell_id:5d}     0         -{surf}  trcl=({x:+.6f} {y:+.6f} {z:+.6f})  imp:n=1  tmp=7.76e-8  fill={cell}0  $ {type} No.{idx+1}"
            cell_lines.append(line)

            current = type
            idx +=1 
            complements.append(cell_id)



    # -----------------------------
    # Print output blocks
    # -----------------------------
    print("\n".join(cell_lines))

    for i in range(0, len(complements), 10):
        chunk = complements[i:i+10]
        line = " ".join(f"#{n}" for n in chunk)
        print(line)


if __name__ == "__main__":
    main()
