from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

def calc(mol):
    weight = round(Descriptors.MolWt(mol),3)
    logp = round(Descriptors.MolLogP(mol),3)
    tpsa = round(Descriptors.TPSA(mol),3)
    hba = round(Lipinski.NumHAcceptors(mol),3)
    hbd = round(Lipinski.NumHDonors(mol),3)
    n = 0   # NN
    o = 0   # NO

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            n += 1
        elif atom.GetSymbol() == 'O':
            o += 1
    no = Lipinski.NOCount(mol)  # N+O

    return [weight,logp,tpsa,hba,hbd,n,o,no]