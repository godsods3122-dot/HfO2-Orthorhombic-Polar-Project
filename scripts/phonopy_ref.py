import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
from phonopy.interface.calculator import read_crystal_structure

def get_ph():
    u,_=read_crystal_structure(filename='POSCAR',interface_mode='vasp')
    ph=Phonopy(u,supercell_matrix=np.diag([2,2,2]),primitive_matrix='P')
    ph.dataset=parse_FORCE_SETS(filename='FORCE_SETS')
    ph.produce_force_constants(calculate_full_force_constants=True)
    ph.symmetrize_force_constants_by_space_group(); ph.symmetrize_force_constants(level=3)
    P=parse_BORN(ph.primitive,filename='BORN'); P['factor']=14.399652; P['method']='gonze'
    ph.nac_params=P
    return ph
