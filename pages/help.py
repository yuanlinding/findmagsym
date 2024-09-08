import streamlit as st

st.header("INPUT")
inp ='''
A MCIF format structure file. The MCIF file is an extended CIF format structure file in which magnetic space group and magnetic moment information are added to describe magnetic material structure. Practially, the MCIF file can be generated from CIF format structure file using [FINDSYM](https://iso.byu.edu/iso/findsym.php) or [STRCONVERT](https://cryst.ehu.es/cgi-bin/cryst/programs/mcif2vesta/index.php) program by supplying the magnetic moments at each individual sites mannually.
'''
st.markdown(inp)

st.header("OUTPUT")
outp ='''
Once a MICF format structure file in imported, the program will genearte the following magnetic symmetry information. 
1. Identified magnetic space group without spin-orbit coupling (MSG without SOC) described by: 
    - MSG symbols given in [BNS setting](https://stokes.byu.edu/iso/magneticspacegroupshelp.php)  
    - MSG type  
    - The set of space and time smmetries that preserves the crystal and its spin structure invariant without considering the spin orbit coupling (such that the spin and spatial space are decoupled). These symmetries consist of three parts (i) proper and improper crystal rotations represented by 3x3 matrixes; and (ii) translations represented by 1x3 matrixes; and (iii) time reversal symmetries
 represented by a list of 0 (absence) or 1 (presence).  

2. Identified magnetic space group with SOC described by: 
    - MSG symbols given in [BNS setting](https://stokes.byu.edu/iso/magneticspacegroupshelp.php)  
    - MSG type  
    - The set of space and time smmetries that preserves the crystal and its spin structure invariant considering the spin orbit coupling. These symmetries consist of three parts (i) proper and improper crystal rotations represented by 3x3 matrixes; and (ii) translations represented by 1x3 matrixes; and (iii) time reversal symmetries represented by a list of 0 (absence) or 1 (presence).  

3. Identified magnetic type, spin splitting behavior, and spin splitting type given by: 
    - COLLINEAR, COPLANNAR, NON-COPLANNAR 
    - YES/NO non-relativistic spin splitting (NRSS), YES/NO relativistic spin splitting (RSS), $\Gamma$-degenerate/$\Gamma$-split  
    - Spin splitting type (SST) following the definition of [Physcial Review Materials 5,014409 (2021)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.014409) 

'''
st.markdown(outp)
