from mordred import Calculator, descriptors ,TopoPSA,SLogP, AtomCount
from mordred.Lipinski import HBondDonor, HBondAcceptor
from rdkit import Chem

if __name__ == '__main__':
    molecules = Chem.MolFromSequence('TRSSRAGLQFPVGRVHRLLRK')

    # all descriptors
    #calc = Calculator(descriptors,ignore_3D=False)
    #df_mordred = calc.pandas([molecules])

    # dichiara descriptors
    logp = Calculator(SLogP)
    tpsa = Calculator(TopoPSA)
    AtomCount = Calculator(AtomCount)
    HBondDonor = Calculator(HBondDonor)
    HBondAcceptor = Calculator(HBondAcceptor)

    print('tpsa\n',tpsa.pandas([molecules])['TopoPSA'],
          '\n\nlogP\n',logp.pandas([molecules])['SLogP'],
          '\n\nAtom Count\n',AtomCount.pandas([molecules])['nAtom'],
          '\n\nNitrogen Count\n', AtomCount.pandas([molecules])['nN'],
          '\n\nOxygen Count\n', AtomCount.pandas([molecules])['nO'],
          '\n\nNitrogen + Oxygen Count\n', AtomCount.pandas([molecules])['nN'] + AtomCount.pandas([molecules])['nO'],
          '\n\nHydrogen Bond Donor\n',HBondDonor.pandas([molecules]),
          '\n\nHydrogen Bond Acceptor\n',HBondAcceptor.pandas([molecules])
          )