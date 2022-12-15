from mordred import Calculator, descriptors
from rdkit import Chem
import os

if __name__ == '__main__':
    basePDB = os.getcwd() + '/data/1pef.pdb'
    molecules = Chem.MolFromPDBFile(basePDB)

    calc = Calculator(descriptors,ignore_3D=True)
    df_mordred = calc.pandas([molecules])

    print(df_mordred)