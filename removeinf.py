'''
 remove inf in line
'''

def RemoveInf(infile):
    outfile = infile + "_noinf"
    outobj  = open(outfile, "w")
    for line in open(infile):
        if not "inf" in line:
            outobj.write(line)
    outobj.close()

if __name__ == "__main__":
    RemoveInf("onlydist.txt")
