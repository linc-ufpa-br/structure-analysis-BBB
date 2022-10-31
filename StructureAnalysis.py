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

def calc(data):
    logp = Descriptors.MolLogP(data)
    tpsa = Descriptors.TPSA(data)
    hba = Lipinski.NumHAcceptors(data)
    hbd = Lipinski.NumHDonors(data)
    n = 0   # NN
    o = 0   # NO
    for atom in data.GetAtoms():
        if atom.GetSymbol() == 'N':
            n += 1
        elif atom.GetSymbol() == 'O':
            o += 1
    no = Lipinski.NOCount(data)  # N+O

    return [logp,tpsa,hba,hbd,n,o,no]

def typeFilesCalc(data):
    # extension file
    ext = (os.path.basename(data)).split('.')[-1]

    if ext == 'pdb':
        dataPDB = Chem.MolFromPDBFile(data)
        PDBresult = pd.DataFrame(calc(dataPDB)).T
        PDBresult.columns = descriptors
        return PDBresult.to_csv('results/PDBresult.csv',index=False)

    elif ext == 'mol':
        dataMol = Chem.MolFromMolFile(data)
        MOLresult = pd.DataFrame(calc(dataMol)).T
        MOLresult.columns = descriptors
        return MOLresult.to_csv('results/MOLresult.csv',index=False)

    elif ext == 'csv':
        dataCSV = pd.read_csv(data)
        peptides = []
        CSVresult = []

        for i in range(len(dataCSV)):
            # removing line break
            dataCSV.iloc[i][0] = str(dataCSV.iloc[i][0]).strip()
            peptides.append(Chem.MolFromSequence(dataCSV.iloc[i][0]))

        # removing invalids
        peptides = pd.DataFrame(filter(lambda item: item is not None, peptides))

        for i in range(len(peptides)):
            CSVresult.append(calc(peptides.iloc[i, 0]))

        CSVresult = pd.DataFrame(CSVresult)
        CSVresult.columns = descriptors
        return CSVresult.to_csv('results/CSVresult.csv', index=False)

    else:
        return 'Unsupported file type.'

typeFilesCalc(baseMol)
typeFilesCalc(basePDB)
typeFilesCalc(baseCSV)
#------------------------------------------------