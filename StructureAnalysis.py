import os
import pandas as pd
from Bio import SeqIO
from rdkit import Chem
from filter import naturalPepFilter
from calcDescriptors import calc

# data for test
basePDB = os.getcwd() + '/data/1pef.pdb'
baseMol = os.getcwd() + '/data/mol_files/BBB+/mol_18.mol'
baseCSV = os.getcwd() + '/data/brainPeps..csv'
baseFASTA = os.getcwd() + '/data/pep.FASTA'

# descriptors to calculate
descriptors = ['logP', 'TPSA(Tot)', 'HBA', 'HBD', 'nN', 'nO', 'n(N+O)']

def typeFilesCalc(data, colPeptides = 'peptides'):
    ext = (os.path.basename(data)).split('.')[-1]

    # for pdb files
    if ext == 'pdb':
        chain = {record.id: record.seq for record in SeqIO.parse(data, 'pdb-seqres')}
        for k, v in chain.items():
            seqPDB = v
            naturalPep = naturalPepFilter(seqPDB)

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
        data = baseCSV

        molCSV = pd.read_csv(data)
        results = []

        '''for i in range(len(molCSV)):
            # removing break lines
            molCSV.iloc[i][colPeptides] = str(molCSV.iloc[i][colPeptides]).strip()'''

        # only natural peptide
        molCSV = [x for x in molCSV[colPeptides] if naturalPepFilter(x) is not None]

        for i in range(len(molCSV)):
            results.append(calc(Chem.MolFromSequence(molCSV[i])))

        CSVresult = pd.concat([pd.DataFrame(molCSV), pd.DataFrame(results)], axis=1)
        CSVresult.columns = ['Peptide'] + descriptors
        return CSVresult.to_csv('results/CSVresult.csv', index=False)

    else:
        return 'Unsupported file type.'

#typeFilesCalc(baseMol)
#typeFilesCalc(basePDB)
#typeFilesCalc(baseCSV,'onelettersequence')
