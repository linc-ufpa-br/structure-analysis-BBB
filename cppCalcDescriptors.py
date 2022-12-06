'''
ja calculado (MW)
ja calculado (tPSA)
ja calculado (Fsp3)
ja calculado (cLogP)
ja calculado (HBA)
ja calculado (HBD)
número de anéis aromáticos (NAR)
numero de ligacoes rotativas (NRB)
carga líquida (NetC)
primary amine groups (NPA)
number of guanidine groups (NG)
number of negatively charged amino acid groups (NNCAA)
amino acid composition (AAC)
pseudo-amino acid composition (PseAAC)
dipeptide composition (DPC)
número de grupos guanidínio (NG)
'''

import os
import pandas as pd
from rdkit import Chem
from filter import naturalPepFilter
from calcDescriptors import calc
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors

def calc(mol):
    mw = round(Descriptors.MolWt(mol), 3)
    logp = round(Descriptors.MolLogP(mol),3)
    tpsa = round(Descriptors.TPSA(mol),3)
    hba = round(Lipinski.NumHAcceptors(mol),3)
    hbd = round(Lipinski.NumHDonors(mol),3)
    n = 0   # NN
    o = 0   # NO
    Fsp3 = Chem.rdMolDescriptors.CalcFractionCSP3(mol)
    nar = Chem.rdMolDescriptors.CalcNumAromaticCarbocycles(mol)

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            n += 1
        elif atom.GetSymbol() == 'O':
            o += 1
    no = Lipinski.NOCount(mol)  # N+O

    return [mw,logp,tpsa,hba,hbd,n,o,no,Fsp3,nar]

base = os.getcwd() + '/data/seq_train.csv'

# descriptors to calculate
descriptors = ['MW','logP', 'TPSA(Tot)', 'HBA', 'HBD', 'nN', 'nO', 'n(N+O)','Fsp3','NAR']

molCSV = pd.read_csv(base)
results = []

# only natural peptide
molCSV = [x for x in molCSV["Seq"] if naturalPepFilter(x) is not None]

for i in range(len(molCSV)):
    results.append(calc(Chem.MolFromSequence(molCSV[i])))

CSVresult = pd.concat([pd.DataFrame(molCSV), pd.DataFrame(results)], axis=1)
CSVresult.columns = ['Peptide'] + descriptors
CSVresult.to_csv('results/CSV_CPP.csv', index=False)