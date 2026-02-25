import pandas as pd
import numpy as np
from numpy.linalg import norm
from spglib import get_magnetic_symmetry, get_magnetic_spacegroup_type_from_symmetry, get_magnetic_symmetry_from_database
from spinspg import get_spin_symmetry
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser


def read_mcif(mcif_path):
    """Read mcif file from path"""
    with open(mcif_path, 'r') as f:
        string_data = f.read()

    try:
        # Try using pymatgen's built-in mcif parsing
        stru = Structure.from_str(string_data, "cif")
        lattice = stru.lattice.matrix
        positions = stru.frac_coords
        numbers = np.array(stru.atomic_numbers)

        # Check if magnetic moments were parsed
        if 'magmom' in stru.site_properties:
            magmoms = np.array([list(stru.site_properties['magmom'][i].moment) for i in range(stru.num_sites)])
        else:
            magmoms = np.zeros((len(stru), 3))

        return lattice, positions, numbers, magmoms
    except:
        # Fallback to manual parsing
        parser = CifParser.from_str(string_data)
        stru = parser.parse_structures()[0]

        lattice = stru.lattice.matrix
        positions = stru.frac_coords
        numbers = np.array(stru.atomic_numbers)

        cif_dict = parser.as_dict()
        data_block = list(cif_dict.values())[0]

        if '_atom_site_moment.label' in data_block:
            labels = data_block.get('_atom_site_moment.label', [])
            mx = data_block.get('_atom_site_moment.crystalaxis_x', [])
            my = data_block.get('_atom_site_moment.crystalaxis_y', [])
            mz = data_block.get('_atom_site_moment.crystalaxis_z', [])

            mag_map = {}
            for i in range(len(labels)):
                # Parse values with uncertainties like '-3.6(2)' -> -3.6
                def parse_value(val):
                    val_str = str(val)
                    if '(' in val_str:
                        return float(val_str.split('(')[0])
                    return float(val_str)

                mag_map[labels[i]] = [parse_value(mx[i]), parse_value(my[i]), parse_value(mz[i])]

            magmoms_list = []
            for site in stru:
                symbol = site.specie.symbol
                magmoms_list.append(mag_map.get(symbol, [0.0, 0.0, 0.0]))

            magmoms = np.array(magmoms_list)
        else:
            magmoms = np.zeros((len(stru), 3))

        return lattice, positions, numbers, magmoms


def find_spinspacegroup(mcif_path):
    lattice, positions, numbers, magmoms = read_mcif(mcif_path)
    sog, rotations, translations, spin_rotations = get_spin_symmetry(lattice, positions, numbers, magmoms)
    return sog, rotations, translations, spin_rotations


def find2_msg_wo_soc(mcif_path):
    lattice, positions, numbers, magmoms = read_mcif(mcif_path)
    sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_path)
    symmetry = {'rotations':[],'translations':[],'time_reversals':[]}
    for i in range(len(rotations)):
        magmoms_rot = []
        for j in range(len(magmoms)):
            rotmag = np.dot(spin_rotations[i],magmoms[j])
            magmoms_rot.append(rotmag)
        if np.all(magmoms-magmoms_rot<0.001):
            symmetry['rotations'].append(rotations[i])
            symmetry['translations'].append(translations[i])
            symmetry['time_reversals'].append(False)
        elif np.all(magmoms+magmoms_rot<0.001):
            symmetry['rotations'].append(rotations[i])
            symmetry['translations'].append(translations[i])
            symmetry['time_reversals'].append(True)
        else:
            continue
    rots = symmetry['rotations']
    trans = symmetry['translations']
    time = symmetry['time_reversals']
    msg_wo_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
    return msg_wo_soc


def find_msg_w_soc(mcif_path):
    lattice, positions, numbers, magmoms = read_mcif(mcif_path)
    symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
    rots = symmetry['rotations']
    trans = symmetry['translations']
    time = symmetry['time_reversals']
    msg_w_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
    return msg_w_soc


def is_compensated_mag(magmoms):
    net_mag = np.array([0.0, 0.0, 0.0])
    for magmom in magmoms:
        net_mag = net_mag + magmom
    if norm(net_mag) < 0.001:
        return True
    else:
        return False


def test_mcif_file(mcif_path):
    print(f"\n{'='*80}")
    print(f"Testing: {mcif_path}")
    print(f"{'='*80}\n")

    # Read structure
    lattice, positions, numbers, magmoms = read_mcif(mcif_path)
    print(f"Number of atoms: {len(positions)}")
    print(f"Magnetic moments:\n{magmoms}\n")

    # Get spin space group
    sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_path)
    print(f"Spin-only group type: {sog.spin_only_group_type}")
    print(f"Number of symmetry operations: {len(rotations)}\n")

    # Load MSG database
    df = pd.read_csv("msg_list.cvs", dtype=str)
    df = df.set_index(['UNI_NUM'])

    # MSG without SOC
    msg_wo_soc = find2_msg_wo_soc(mcif_path)
    bns_symbol = df.loc[df['BNS_NUM']==msg_wo_soc.bns_number].iloc[0]['BNS_SYM']
    print(f"MSG without SOC:")
    print(f"  BNS Symbol: {bns_symbol}")
    print(f"  BNS Number: {msg_wo_soc.bns_number}")
    print(f"  UNI Number: {msg_wo_soc.uni_number}")
    print(f"  MSG Type: {msg_wo_soc.type}\n")

    # MSG with SOC
    msg_w_soc = find_msg_w_soc(mcif_path)
    bns_symbol = df.loc[df['BNS_NUM']==msg_w_soc.bns_number].iloc[0]['BNS_SYM']
    print(f"MSG with SOC:")
    print(f"  BNS Symbol: {bns_symbol}")
    print(f"  BNS Number: {msg_w_soc.bns_number}")
    print(f"  UNI Number: {msg_w_soc.uni_number}")
    print(f"  MSG Type: {msg_w_soc.type}\n")

    # Determine material type
    is_compensated = is_compensated_mag(magmoms)
    print(f"Compensated magnet: {is_compensated}")
    print(f"Spin structure type: {sog.spin_only_group_type}")


if __name__ == "__main__":
    # Test the two mcif files
    mcif_files = [
        "../0.222_CuMnAs.mcif",
        "../1.6_NiO.mcif"
    ]

    for mcif_file in mcif_files:
        try:
            test_mcif_file(mcif_file)
        except Exception as e:
            print(f"Error testing {mcif_file}: {e}")
            import traceback
            traceback.print_exc()
