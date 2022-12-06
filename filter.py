oneLetterCode = ['A', 'R', 'N', 'D', 'B',
                     'C','E', 'Q', 'Z', 'G',
                     'H','I','L', 'K', 'M', 'F',
                     'P','S','T', 'W','Y', 'V']

# natural peptides filter (returns the sequence)
def naturalPepFilter(peptide):
    verify = True
    for letter in peptide:
        verify = letter in oneLetterCode
        if verify is False:
            return None
    if verify is True:
        return peptide

# natural and modified peptides filter (returns natural or unnatural)
def nmPepFilter(peptide):
    verify = True
    for letter in peptide:
        verify = letter in oneLetterCode
        if verify is False:
            return 'M'
    if verify is True:
        return 'N'
