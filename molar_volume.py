import math
from rdkit import Chem

def vdw_volume(smiles):

    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    atom_volumes = {
        'H': 7.24, 'C': 20.58, 'N': 15.6, 'O': 14.71, 'F': 13.31,
        'Cl': 22.45, 'Br': 26.52, 'I': 32.52, 'P': 24.43, 'S': 24.43,
        'As': 26.52, 'B': 40.48, 'Si': 38.79, 'Se': 28.73, 'Te': 36.62 # values from Zhao YH, Abraham MH, Zissimos AM. J Org Chem. 2003 68(19):7368-73. doi: 10.1021/jo034808o.
        }
    alvarez_radii = {
        'H': 1.2, 'He': 1.43, 'Li': 2.12, 'Be': 1.98, 'B': 1.91, 'C': 1.77, 'N': 1.66, 'O': 1.5, 'F': 1.46, 'Ne': 1.58, 'Na': 2.5, 'Mg': 2.51, 'Al': 2.25, 'Si': 2.19, 'P': 1.9, 'S': 1.89, 'Cl': 1.82, 'Ar': 1.83,
        'K': 2.73, 'Ca': 2.62, 'Sc': 2.58, 'Ti': 2.46, 'V': 2.42, 'Cr': 2.45, 'Mn': 2.45, 'Fe': 2.44, 'Co': 2.4, 'Ni': 2.4, 'Cu': 2.38, 'Zn': 2.39, 'Ga': 2.32, 'Ge': 2.29, 'As': 1.88, 'Se': 1.82, 'Br': 1.86, 'Kr': 2.25,
        'Rb': 3.21, 'Sr': 2.84, 'Y': 2.75, 'Zr': 2.52, 'Nb': 2.56, 'Mo': 2.45, 'Tc': 2.44, 'Ru': 2.46, 'Rh': 2.44, 'Pd': 2.15, 'Ag': 2.53, 'Cd': 2.49, 'In': 2.43, 'Sn': 2.42, 'Sb': 2.47, 'Te': 1.99, 'I': 2.04, 'Xe': 2.06,
        'Cs': 3.48, 'Ba': 3.03, 'La': 2.98, 'Ce': 2.88, 'Pr': 2.92, 'Nd': 2.95, 'Sm': 2.9, 'Eu': 2.87, 'Gd': 2.83, 'Tb': 2.79, 'Dy': 2.87, 'Ho': 2.81, 'Er': 2.83, 'Tm': 2.79, 'Yb': 2.8, 'Lu': 2.74, 'Hf': 2.63, 'Ta': 2.53,
        'W': 2.57, 'Re': 2.49, 'Os': 2.48, 'Ir': 2.41, 'Pt': 2.29, 'Au': 2.32, 'Hg': 2.45, 'Tl': 2.47, 'Pb': 2.6, 'Bi': 2.54, 'Ac': 2.8, 'Th': 2.93, 'Pa': 2.88, 'U': 2.71, 'Np': 2.82, 'Pu': 2.81, 'Am': 2.83, 'Cm': 3.05,
        'Bk': 3.4, 'Cf': 3.05, 'Es': 2.7 # values from Alvarez S. Dalton Trans. 2013 42:8617-8636. doi: 10.1039/C3DT50599E.
        }
    for element, radius in alvarez_radii.items():
        if element not in atom_volumes:
            atom_volumes[element] = (4 / 3) * math.pi * (radius ** 3)
    atom_sum = sum(atom_volumes.get(atom.GetSymbol(), 0) for atom in mol_h.GetAtoms()) # zero volume for unknown elements, be careful
    num_bonds = mol_h.GetNumBonds()
    ring_info = mol.GetRingInfo()
    num_aromatic_rings = 0
    num_nonaromatic_rings = 0
    
    for ring_atoms in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            num_aromatic_rings += 1
        else:
            num_nonaromatic_rings += 1
    
    vdw_volume = atom_sum - (5.92 * num_bonds) - (14.7 * num_aromatic_rings) - (3.8 * num_nonaromatic_rings)
    return vdw_volume

def il_molar_volume(smiles, num_ions = None): # Only for ionic liquids, do not expect good accuracy for molecular liquids
    if num_ions is None:
        num_ions = smiles.count('.')+1
    vdw_volume_cm3mol = vdw_volume(smiles)*6.02214/10
    il_molar_volume = 1.58718928*vdw_volume_cm3mol-4.87034696*num_ions
    return il_molar_volume
