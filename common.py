'''
    Common used functions
'''


def length(chain):
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    counter = 0
    for res in chain.child_list:
        if aas.count(res.resname):
            counter += 1
    return counter
