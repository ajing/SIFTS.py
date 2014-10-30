'''
    Get the relation of sequence distance and spatial distance
'''

def GetSpatialDistance(res1, res2):
    return res1['CA'] - res2['CA']

def GetSeqDistance():
    pass

def GetSeqSpaPair(protein):
    reslist  = protein.getResidues()
    pairlist = []
    for res1 in reslist:
        for res2 in reslist:
            dist = res1.getSeqDistance(res2)
            if dist:
                pairlist.append((GetSpatialDistance(res1, res2), dist))
    return pairlist

def PlotPairList(pairlist):




#given pdb find number of residues in chain
def length(chain):
    aas = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'SER', 'THR', 'ASN', 'GLN',
           'PHE', 'TYR', 'TRP', 'CYS', 'MET', 'PRO', 'ASP', 'GLU', 'LYS',
           'ARG', 'HIS']
    counter = 0
    for res in chain.child_list:
        if aas.count(res.resname):
            counter += 1
    return counter
