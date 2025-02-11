import streamlit as st
import pandas as pd
from io import StringIO
import numpy as np
from numpy.linalg import det,norm
from spglib import get_magnetic_symmetry, get_magnetic_spacegroup_type_from_symmetry, get_magnetic_symmetry_from_database, get_symmetry
from spinspg import get_spin_symmetry
from pymatgen.core import Structure

def read_mcif(mcif_file):
	stringio = StringIO(mcif_file.getvalue().decode("utf-8"))
	string_data = stringio.read()
	stru = Structure.from_str(string_data,"cif")
	lattice = stru.lattice.matrix
	positions = stru.frac_coords
	numbers = np.array(stru.atomic_numbers)
	magmoms = np.array([list(stru.site_properties['magmom'][i].moment) for i in range(stru.num_sites)])
	return lattice, positions, numbers, magmoms

def find_spinspacegroup(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	sog, rotations, translations, spin_rotations = get_spin_symmetry(lattice, positions, numbers, magmoms)
	#st.write(f"Spin-only group: {sog}")  # COLLINEAR(axis=[0. 0. 1.])
	return sog, rotations, translations, spin_rotations

def find_msg_wo_soc(mcif_file):
	sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
	symmetry = {'rotations':[],'translations':[],'time_reversals':[]}
	for i in range(len(rotations)):
		symmetry['rotations'].append(rotations[i])
		symmetry['translations'].append(translations[i])
		symmetry['time_reversals'].append(det(spin_rotations[i])==-1)
	# Find MSG without SOC
	rots = symmetry['rotations']
	trans = symmetry['translations']
	time = symmetry['time_reversals']
	msg_wo_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
	#st.write(f"MSG without SOC: {msg_wo_soc}")
	return msg_wo_soc


def find2_msg_wo_soc(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
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
	# Find MSG without SOC
	rots = symmetry['rotations']
	trans = symmetry['translations']
	time = symmetry['time_reversals']
	msg_wo_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
	#st.write(f"MSG without SOC: {msg_wo_soc}")
	return msg_wo_soc

def find_msg_w_soc(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
	rots = symmetry['rotations']
	trans = symmetry['translations']
	time = symmetry['time_reversals']
	msg_w_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
	#st.write(f"MSG with SOC: {msg_w_soc}")
	return msg_w_soc

def find2_msg_w_soc(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
	symmetry = {'rotations':[],'translations':[],'time_reversals':[]}
	for i in range(len(rotations)):
		magmoms_rot = []
		for j in range(len(magmoms)):
			rotmag = np.dot(rotations[i],magmoms[j])
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
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
	rots = symmetry['rotations']
	trans = symmetry['translations']
	time = symmetry['time_reversals']
	msg_w_soc = get_magnetic_spacegroup_type_from_symmetry(rots, trans, time)
	#st.write(f"MSG with SOC: {msg_w_soc}")
	return msg_w_soc

def is_Centrosymmetric(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	symmetry = get_symmetry((lattice, positions, numbers))
	rots = symmetry['rotations']
	for r in rots:
		if np.all(r+[[-1,0,0],[0,-1,0],[0,0,-1]]==0): 
			return True
		else:
			continue
	return False
	

def has_ThetaI(mcif_file):
	lattice, positions, numbers, magmoms = read_mcif(mcif_file)
	symmetry = get_magnetic_symmetry((lattice, positions, numbers, magmoms))
	rots = symmetry['rotations']
	time = symmetry['time_reversals']
	for r,t in zip(rots,time):
		if np.all(r+[[1,0,0],[0,1,0],[0,0,1]]==0) and t:
			return True
		else:
			continue
	return False

def is_compensated_mag(magmoms):
	net_mag = [0,0,0]
	for magmom in magmoms:
		net_mag = net_mag + magmom
	if norm(net_mag) < 0.001: 
		return True
	else:
		return False

def main():
	st.title("FINDMAGSYM")
	multi ='''Version 1.0, Aug 2024  
Linding Yuan, James Rondinelli, Department of Materials Science and Engineering, Northwestern University, Evanston, Illinois  60208, USA
	'''
	st.markdown("Version 1.0, Aug 2024")
	st.markdown("**Description:** findmagsym is a tool to identify the magnetic space group (with and without spin-orbit coupling) of a magnetic crystal, given the positions and magnetic moments of the atoms in a unit cell.")
	st.markdown("**How to cite:** Lin-Ding Yuan, Alexandru B. Georgescu, and James M. Rondinelli. *Nonrelativistic Spin Splitting at the Brillouin Zone Center in Compensated Magnets*, Phys. Rev. Lett. 133, 216701 (2024).")
	st.page_link("pages/help.py",icon=":material/help:")

	st.header("INPUT")
	mcif_file = st.file_uploader("Import structure from a mcif file",type="mcif")
	df = pd.read_csv("msg_list.cvs",dtype=str)
	df = df.set_index(['UNI_NUM'])
	if mcif_file is not None:
		st.header("OUTPUT")
		lattice, positions, numbers, magmoms = read_mcif(mcif_file)
		sog, rotations, translations, spin_rotations = find_spinspacegroup(mcif_file)
		#st.write(magmoms)
		#st.write(f"Spin-only group: {sog}")
		#st.write(f"{sog.spin_only_group_type} spin")

		#MSG without SOC
		msg_wo_soc = find2_msg_wo_soc(mcif_file)
		bns_symbol = df.loc[df['BNS_NUM']==msg_wo_soc.bns_number].iloc[0]['BNS_SYM']
		st.markdown('''**MSG without SOC:**''')
		st.markdown(f"${bns_symbol}$ (BNS);   MSG Type {msg_wo_soc.type}")
		#st.write(msg_wo_soc)
		msg_symm = get_magnetic_symmetry_from_database(msg_wo_soc.uni_number)
		with st.container(height=300):
			st.write(msg_symm)

		#MSG with SOC
		msg_w_soc = find2_msg_w_soc(mcif_file)
		bns_symbol = df.loc[df['BNS_NUM']==msg_w_soc.bns_number].iloc[0]['BNS_SYM']
		st.markdown("**MSG with SOC**")
		st.markdown(f"${bns_symbol}$ (BNS);   MSG Type {msg_w_soc.type}")
		msg_symm = get_magnetic_symmetry_from_database(msg_w_soc.uni_number)
		with st.container(height=300):
			st.write(msg_symm)

		if msg_w_soc.type == 2 and is_Centerosymmetric(mcif_file):
			sst_key = "SST-6"  
		elif msg_w_soc.type == 2 and not is_Centerosymmetric(mcif_file):
			sst_key = "SST-7"  
		elif msg_wo_soc.type == 3 and has_ThetaI(mcif_file):  
			sst_key = "SST-1"
		elif msg_wo_soc.type == 4 and has_ThetaI(mcif_file):
			sst_key = "SST-2"
		elif msg_wo_soc.type == 4 and not has_ThetaI(mcif_file):
			sst_key = "SST-3"
		elif msg_wo_soc.type == 3 and not has_ThetaI(mcif_file):
			sst_key = "SST-4x"
		elif msg_wo_soc.type == 1 and is_compensated_mag(magmoms):
			sst_key = "SST-4y"
		elif msg_wo_soc.type == 1 and not is_compensated_mag(magmoms):
			sst_key = "SST-5"
		else:
			print("Unknown type!!!")
		
		sst = {'SST-1':'antiferromanget (no NRSS, no RSS)',\
			'SST-2':'antiferromagnet (no NRSS, no RSS)',\
			'SST-3':'antiferromagnet (no NRSS, yes RSS)',\
			'SST-4x':r"$\Gamma$-degenerate NRSS antiferromagnet",\
			'SST-4y':r"$\Gamma$-split NRSS antiferromagnet",\
			'SST-5':"ferromagnet or ferrimagnet",\
			'SST-6':"centerosymmetric nonmagent",\
			'SST-7':"non-centrosymmetric nonmagnet"}
		if str(sog.spin_only_group_type) == "COLLINEAR":
			if sst_key in ['SST-4x','SST-4y']:
				st.markdown(f"This material is a :blue[{sog.spin_only_group_type}] :green[{sst[sst_key]}] :red[(SST-4)]")
			else:
				st.markdown(f"This material is a :blue[{sog.spin_only_group_type}] :green[{sst[sst_key]}] :red[({sst_key})]")
		elif is_compensated_mag(magmoms):
			st.markdown(f"This material is a :blue[{sog.spin_only_group_type}] :green[antiferromagnet]")
		else:
     			st.markdown(f"This material is a :blue[{sog.spin_only_group_type}] :green[ferromagnet or ferrimagnet]")

if __name__ == "__main__":
	main()
