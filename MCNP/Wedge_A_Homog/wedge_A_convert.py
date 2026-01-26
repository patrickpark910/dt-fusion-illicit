import re

def convert_openmc_to_mcnp(input_path, mat_id="2", suffix=".02c"):

    """
    Extracts density and converts nuclides to an MCNP material card.
    """

    with open(input_path, 'r') as f:
        input_file = f.read()

    # Atom numbers
    z_map = {
        'He': 2, 'Li': 3, 'Be': 4, 'C': 6, 'N': 7, 'O': 8, 
        'Si': 14, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'W': 74, 'U': 92
    }

    # Extract density
    d_match = re.search(r'density[^>]+value="([\d.]+)"', input_file)
    density = float(d_match.group(1)) if d_match else None

    # Extract nuclides: finds ao="val" and name="SymA" (e.g., Li6)
    nuclides = re.findall(r'nuclide\s+ao="([\d.eE+-]+)"\s+name="([a-zA-Z]+)(\d+)"', input_file)

    # 3. Process data for mass fraction calculation (Mass = Atom Fraction * A)
    processed = []
    for ao_str, symbol, a_str in nuclides:
        ao, a = float(ao_str), int(a_str)
        z = z_map.get(symbol, 0)
        zaid = f"{z * 1000 + a}{suffix}"
        # print(zaid, ao, f"{symbol}-{a}")
        processed.append({'zaid': zaid, 'ao': ao, 'name': f"{symbol}-{a}"})

    # Format into MCNP card
    lines = []
    for i, item in enumerate(processed):
        prefix = f"m{mat_id}" if i == 0 else ""
        # Format: ID (8s), ZAID (10s), -Fraction (15.11f), Comment
        line = f"{prefix:<8} {item['zaid']:>10}  {item['ao']:14.11e}  $ {item['name']}"
        lines.append(line)

    return density, "\n".join(lines)