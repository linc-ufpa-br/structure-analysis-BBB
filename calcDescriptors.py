import os
import pandas as pd
from Bio import SeqIO
from rdkit import Chem
import filter
from calcRdkit import calc

# data for test
basePDB = os.getcwd() + '/data/1pef.pdb'
baseMol = os.getcwd() + '/data/mol_files/BBB+/mol_18.mol'
baseCSV = os.getcwd() + '/data/brainPeps.csv'

# descriptors calculated in dichiara dataset
descriptors = ['logP', 'TPSA(Tot)', 'HBA', 'HBD', 'nN', 'nO', 'n(N+O)']

# folder to save results
if not os.path.exists("results"):
    os.makedirs("results")

def typeFilesCalc(data, colPeptides = 'peptides', colKeep = None):
    ext = (os.path.basename(data)).split('.')[-1]

    # for pdb files
    if ext == 'pdb':
        chain = {record.id: record.seq for record in SeqIO.parse(data, 'pdb-seqres')}
        for k, v in chain.items():
            seqPDB = v
            naturalPep = filter.naturalPep(seqPDB)

        if naturalPep is not False:
            molPDB = Chem.MolFromPDBFile(data)
            PDBresult = pd.DataFrame(calc(molPDB)).T
            PDBresult.columns = descriptors
            return PDBresult.to_csv('results/PDBresult.csv', index=False)
        else:
            print("Error: Not natural peptide!")

    # for mol files (cannot extract sequence for filter)
    elif ext == 'mol':
        molMol = Chem.MolFromMolFile(data)
        MOLresult = pd.DataFrame(calc(molMol)).T
        MOLresult.columns = descriptors
        return MOLresult.to_csv('results/MOLresult.csv',index=False)

    # for csv files
    elif ext == 'csv':
        molCSV = pd.read_csv(data)
        results = []

        for i in range(len(molCSV)):
            # removing break lines (don't always need)
            molCSV.iloc[i][colPeptides] = str(molCSV.iloc[i][colPeptides]).strip()

        # only natural peptide
        naturalmolCSV = []
        for x,y in zip(molCSV[colPeptides],molCSV[colKeep]):
            if filter.naturalPep(x) is not None:
                naturalmolCSV.append([x,y])
        naturalmolCSV = pd.DataFrame(naturalmolCSV)

        for i in range(len(naturalmolCSV)):
            results.append(calc(Chem.MolFromSequence(naturalmolCSV[0][i])))

        if not os.path.exists("results"):
            os.makedirs("results")

        CSVresult = pd.concat([naturalmolCSV, pd.DataFrame(results)], axis=1)
        CSVresult.columns = [colPeptides] + [colKeep] + descriptors
        print("The results refer only to natural peptides.")
        return CSVresult.to_csv('results/CSVresult.csv', index=False)
    else:
        return 'Unsupported file type.'

if __name__ == '__main__':
    typeFilesCalc(baseMol)
    typeFilesCalc(basePDB)
    typeFilesCalc(baseCSV,'onelettersequence','label')
