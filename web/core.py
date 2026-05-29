"""Pure analysis functions — no Streamlit dependency."""

import numpy as np
from numpy.linalg import det, norm
from spglib import (
    get_magnetic_symmetry,
    get_magnetic_spacegroup_type_from_symmetry,
    get_magnetic_symmetry_from_database,
    get_symmetry,
)
from spinspg import get_spin_symmetry
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser


def read_mcif(mcif_file):
    """Parse an mcif/cif file. Accepts a file-like object (with .getvalue() or
    .read()) or a plain string of CIF content."""
    try:
        if hasattr(mcif_file, 'getvalue'):
            bytes_data = mcif_file.getvalue()
            if isinstance(bytes_data, bytes):
                try:
                    string_data = bytes_data.decode("utf-8")
                except UnicodeDecodeError:
                    string_data = bytes_data.decode("latin-1")
            else:
                string_data = bytes_data
        elif hasattr(mcif_file, 'read'):
            content = mcif_file.read()
            if isinstance(content, bytes):
                try:
                    string_data = content.decode("utf-8")
                except UnicodeDecodeError:
                    string_data = content.decode("latin-1")
            else:
                string_data = content
        else:
            string_data = str(mcif_file)
    except Exception as e:
        raise ValueError(f"Could not read MCIF file: {e}")

    try:
        stru = Structure.from_str(string_data, "cif")
        lattice = stru.lattice.matrix
        positions = stru.frac_coords
        numbers = np.array(stru.atomic_numbers)
        if 'magmom' in stru.site_properties:
            magmoms = np.array(
                [list(stru.site_properties['magmom'][i].moment)
                 for i in range(stru.num_sites)]
            )
        else:
            magmoms = np.zeros((stru.num_sites, 3))
        return lattice, positions, numbers, magmoms

    except Exception:
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

            def parse_value(val):
                val_str = str(val)
                if '(' in val_str:
                    return float(val_str.split('(')[0])
                return float(val_str)

            mag_map = {}
            for i in range(len(labels)):
                mag_map[labels[i]] = [parse_value(mx[i]), parse_value(my[i]), parse_value(mz[i])]
            magmoms_list = []
            for site in stru:
                symbol = site.specie.symbol
                magmoms_list.append(mag_map.get(symbol, [0.0, 0.0, 0.0]))
            magmoms = np.array(magmoms_list)
        else:
            magmoms = np.zeros((len(stru), 3))
        return lattice, positions, numbers, magmoms


def find_spinspacegroup(mcif_file):
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    sog, rotations, translations, spin_rotations = get_spin_symmetry(
        lattice, positions, numbers, magmoms
    )
    return sog, rotations, translations, spin_rotations


def find_msg_wo_soc(mcif_file):
    sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
    symmetry = {'rotations': [], 'translations': [], 'time_reversals': []}
    for i in range(len(rotations)):
        symmetry['rotations'].append(rotations[i])
        symmetry['translations'].append(translations[i])
        symmetry['time_reversals'].append(det(spin_rotations[i]) == -1)
    return get_magnetic_spacegroup_type_from_symmetry(
        symmetry['rotations'], symmetry['translations'], symmetry['time_reversals']
    )


def find2_msg_wo_soc(mcif_file):
    """Improved MSG-without-SOC finder with better non-magnet handling."""
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
    symmetry = {'rotations': [], 'translations': [], 'time_reversals': []}
    for i in range(len(rotations)):
        magmoms_rot = [np.dot(spin_rotations[i], magmoms[j]) for j in range(len(magmoms))]
        if np.all(np.array(magmoms) - np.array(magmoms_rot) < 0.001):
            symmetry['rotations'].append(rotations[i])
            symmetry['translations'].append(translations[i])
            symmetry['time_reversals'].append(False)
        elif np.all(np.array(magmoms) + np.array(magmoms_rot) < 0.001):
            symmetry['rotations'].append(rotations[i])
            symmetry['translations'].append(translations[i])
            symmetry['time_reversals'].append(True)
    return get_magnetic_spacegroup_type_from_symmetry(
        symmetry['rotations'], symmetry['translations'], symmetry['time_reversals']
    )


def find_msg_w_soc(mcif_file):
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
    return get_magnetic_spacegroup_type_from_symmetry(
        symmetry['rotations'], symmetry['translations'], symmetry['time_reversals']
    )


def is_Centrosymmetric(mcif_file):
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    symmetry = get_symmetry((lattice, positions, numbers))
    rots = symmetry['rotations']
    inversion = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    for r in rots:
        if np.all(r == inversion):
            return True
    return False


def has_ThetaI(mcif_file):
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
    rots = symmetry['rotations']
    time = symmetry['time_reversals']
    for r, t in zip(rots, time):
        if np.all(r + np.eye(3, dtype=int) == 0) and t:
            return True
    return False


def is_compensated_mag(magmoms):
    net_mag = np.sum(np.array(magmoms), axis=0)
    return norm(net_mag) < 0.001


SST_DESCRIPTIONS = {
    'SST-1': 'antiferromagnet (no NRSS, no RSS)',
    'SST-2': 'antiferromagnet (no NRSS, no RSS)',
    'SST-3': 'antiferromagnet (no NRSS, yes RSS)',
    'SST-4x': 'Γ-degenerate NRSS antiferromagnet',
    'SST-4y': 'Γ-split NRSS antiferromagnet',
    'SST-5': 'ferromagnet or ferrimagnet',
    'SST-6': 'centrosymmetric nonmagnet',
    'SST-7': 'non-centrosymmetric nonmagnet',
}


def classify_sst(mcif_file):
    """Return (sst_key, sog) for the material in mcif_file."""
    lattice, positions, numbers, magmoms = read_mcif(mcif_file)
    sog, _, _, _ = find_spinspacegroup(mcif_file)
    msg_wo_soc = find2_msg_wo_soc(mcif_file)
    msg_w_soc = find_msg_w_soc(mcif_file)

    centrosym = is_Centrosymmetric(mcif_file)
    theta_i = has_ThetaI(mcif_file)
    compensated = is_compensated_mag(magmoms)

    if msg_w_soc.type == 2 and centrosym:
        sst_key = "SST-6"
    elif msg_w_soc.type == 2 and not centrosym:
        sst_key = "SST-7"
    elif msg_wo_soc.type == 3 and theta_i:
        sst_key = "SST-1"
    elif msg_wo_soc.type == 4 and theta_i:
        sst_key = "SST-2"
    elif msg_wo_soc.type == 4 and not theta_i:
        sst_key = "SST-3"
    elif msg_wo_soc.type == 3 and not theta_i:
        sst_key = "SST-4x"
    elif msg_wo_soc.type == 1 and compensated:
        sst_key = "SST-4y"
    elif msg_wo_soc.type == 1 and not compensated:
        sst_key = "SST-5"
    else:
        sst_key = "UNKNOWN"

    return sst_key, sog, msg_wo_soc, msg_w_soc, centrosym, theta_i, compensated
