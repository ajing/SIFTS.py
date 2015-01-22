'''
 convert humsavar so MySQL database can read
'''

import collections

def Convert(infile):
    fileobj = open(infile)
    outobj  = open(infile + ".out_test_1", "w")
    for line in fileobj:
        content = line.strip().split()
        aachange =  String2Info(content[3])
        newline = content[:6]
        newline.append(" ".join(content[6:]))
        newline[3] = aachange
        newline    = list(flatten(newline))
        outobj.write("\t".join(newline) + "\n")
    outobj.close()

def String2Info(astring):
    content = astring.split(".")[1]
    res_before = content[:3].upper()
    res_after = content[-3:].upper()
    res_num   = content[3:-3]
    return res_before, res_num, res_after

def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


if __name__ == "__main__":
    Convert("humsavar_tru.txt")
