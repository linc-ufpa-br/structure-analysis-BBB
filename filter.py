# natural peptides filter
def mn(peptide):
    oneLetterCode = ['A', 'R', 'N', 'D', 'B',
                     'C','E', 'Q', 'Z', 'G',
                     'H','I','L', 'K', 'M', 'F',
                     'P','S','T', 'W','Y', 'V']
    verify = True
    for letter in peptide:
        verify = letter in oneLetterCode
        if verify is False:
            return None
    if verify is True:
        return peptide