from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
import os
import pandas as pd
dir = os.getcwd()

# data for test
basePDB = dir + '/data/1pef.pdb'
baseMol = dir + '/data/mol_files/BBB+/mol_18.mol'
baseCSV = dir + '/data/brainPeps..csv'
baseFASTA = dir + '/data/pep.FASTA'

# descriptors to calculate
descriptors = ['logP', 'TPSA(Tot)', 'HBA', 'HBD', 'nN', 'nO', 'n(N+O)']

def calc(mol):
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

    return [logp,tpsa,hba,hbd,n,o,no]

def typeFilesCalc(data):
    # extension file
    ext = (os.path.basename(data)).split('.')[-1]

    if ext == 'pdb':
        molPDB = Chem.MolFromPDBFile(data)
        PDBresult = pd.DataFrame(calc(molPDB)).T
        PDBresult.columns = descriptors
        return PDBresult.to_csv('results/PDBresult.csv',index=False)

    elif ext == 'mol':
        molMol = Chem.MolFromMolFile(data)
        MOLresult = pd.DataFrame(calc(molMol)).T
        MOLresult.columns = descriptors
        return MOLresult.to_csv('results/MOLresult.csv',index=False)

    elif ext == 'csv':
        molCSV = pd.read_csv(data)
        peptides = []
        CSVresult = []

        for i in range(len(molCSV)):
            # removing break lines
            molCSV.iloc[i][0] = str(molCSV.iloc[i][0]).strip()
            peptides.append([molCSV.iloc[i][0],Chem.MolFromSequence(molCSV.iloc[i][0])])

        # removing invalids
        peptides = pd.DataFrame(peptides).dropna()

        for i in range(len(peptides)):
            CSVresult.append(calc(peptides.iloc[i][1]))

        CSVresult = pd.DataFrame(CSVresult)
        CSVresult = pd.concat([peptides[0],CSVresult], axis=1)
        CSVresult.columns = ['Peptide'] + descriptors

        return CSVresult.to_csv('results/CSVresult.csv', index=False)

    else:
        return 'Unsupported file type.'

typeFilesCalc(baseMol)
typeFilesCalc(basePDB)
typeFilesCalc(baseCSV)
#------------------------------------------------