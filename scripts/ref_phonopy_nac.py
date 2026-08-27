import numpy as np
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
from phonopy.interface.calculator import read_crystal_structure

unitcell, _ = read_crystal_structure(filename="POSCAR", interface_mode="vasp")
ph = Phonopy(unitcell, supercell_matrix=np.diag([2,2,2]), primitive_matrix='P')
fs = parse_FORCE_SETS(filename="FORCE_SETS")
ph.dataset = fs
ph.produce_force_constants(calculate_full_force_constants=True)
ph.symmetrize_force_constants_by_space_group()
ph.symmetrize_force_constants(level=3)

nac_params = parse_BORN(ph.primitive, filename="BORN")
nac_params["factor"] = 14.399652
nac_params["method"] = "gonze"
ph.nac_params = nac_params

# G -> X path, with q_direction at Gamma set toward X (1,0,0) for the LO-TO limit
qs = [[0.5*i/20.0, 0.0, 0.0] for i in range(21)]
ph.run_qpoints(qs, with_group_velocities=False)
freqs = ph.qpoints.frequencies

# Gamma with explicit q_direction toward X
from phonopy.phonon.qpoints import QpointsPhonon
qph = QpointsPhonon([[0.,0.,0.]], ph.dynamical_matrix, nac_q_direction=[1.,0.,0.], factor=ph.unit_conversion_factor)
gamma_dir_freq = qph.frequencies[0]

np.savetxt("phonopy_ref_GX_NAC.dat", freqs, fmt="%.6f")
print("Gamma (q_direction=None, phonopy default) freqs:", np.round(freqs[0],4))
print("Gamma (q_direction=X)                    freqs:", np.round(gamma_dir_freq,4))
print("X freqs:", np.round(freqs[-1],4))
print("point near Gamma (i=1) freqs:", np.round(freqs[1],4))
