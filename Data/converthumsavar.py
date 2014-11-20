'''
 convert humsavar so MySQL database can read
'''

def Convert(infile):
    fileobj = open(infile)
    outobj  = open(infile + ".out", "w")
    for line in fileobj:
        content = line.strip().split()
        newline = content[:6]
        newline.append(" ".join(content[6:]))
        outobj.write("\t".join(newline) + "\n")
    outobj.close()


if __name__ == "__main__":
    Convert("humsavar_tru.txt")
