from mordred import Calculator, descriptors,TopoPSA,SLogP, Weight
import mordred #some functions do not work importing specifically or even importing whole mordred
from mordred import Lipinski

# all descriptors
def calcAllDescriptors(mol):
    calc = Calculator(descriptors,ignore_3D=False)
    result = calc.pandas([mol])
    return result

# dichiara descriptors + weight
def calcSpecific(mol):
    weight = Calculator(Weight)
    logp = Calculator(SLogP)
    tpsa = Calculator(TopoPSA)
    AtomCount = Calculator(mordred.AtomCount)
    HBondDonor = Calculator(Lipinski.HBondDonor)
    HBondAcceptor = Calculator(Lipinski.HBondAcceptor)

    resultWeight = round(float(weight.pandas([mol])['MW']),3)
    resultLogp = round(float(logp.pandas([mol])['SLogP']),3)
    resultTpsa = round(float(tpsa.pandas([mol])['TopoPSA']),3)
    resultAtomCount = round(float(AtomCount.pandas([mol])['nAtom']),3)
    resultNitrogenCount = round(float(AtomCount.pandas([mol])['nN']),3)
    resultOxygenCount = round(float(AtomCount.pandas([mol])['nO']),3)
    resultNOCount = round(float(AtomCount.pandas([mol])['nN'] + AtomCount.pandas([mol])['nO']),3)
    resultHBondDonor = round(float(HBondDonor.pandas([mol])['nHBDon']),3)
    resultHBondAcceptor = round(float(HBondAcceptor.pandas([mol])['nHBAcc']),3)

    return[resultWeight,resultLogp,resultTpsa,resultHBondAcceptor,resultHBondDonor,resultNitrogenCount,resultOxygenCount,resultNOCount]