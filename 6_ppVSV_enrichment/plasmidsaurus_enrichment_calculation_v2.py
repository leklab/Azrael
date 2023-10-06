
import sys

# input the files

infile=open(sys.argv[1],'rt')
sta=int(sys.argv[2])
total=int(sys.argv[3])
rev=int(sys.argv[4])
cut=int(sys.argv[5])

# read the sequence

def revornot(x):
    if rev==1:
        return(total-x+1)
    if rev==0:
        return(x)

def output(x):
    if rev==1:
        if x==7:
            return(9)
        if x==9:
            return(7)
        if x==8:
            return(10)
        if x==10:
            return(8)
    if rev==0:
        return(x)

def score(x, site, wt, v, total_s, no):
    if x==str(revornot(site+sta)):
        print("c."+str(site)+ wt + ">" + v +"\t"+str(no/total_s))


inlines=[line.rstrip('\n') for line in infile]



for i in range(cut, len(inlines)):
    tl=inlines[i].split(",")
    total_s=float(tl[7])+float(tl[8])+float(tl[9])+float(tl[10])
    no7=float(tl[output(7)])
    no8=float(tl[output(8)])
    no9=float(tl[output(9)])
    no10=float(tl[output(10)])
    score(tl[0], 435, "C", "T", total_s, no9)
    score(tl[0], 992, "C", "T", total_s, no9)
    score(tl[0], 639, "T", "C", total_s, no8)
    score(tl[0], 642, "A", "G", total_s, no10)
    score(tl[0], 690, "T", "A", total_s, no7)
    score(tl[0], 693, "T", "G", total_s, no10)
    score(tl[0], 699, "C", "G", total_s, no10)
    score(tl[0], 26, "G", "C", total_s, no8)
    score(tl[0], 1050, "C", "G", total_s, no10)
    score(tl[0], 1720, "G", "T", total_s, no9)
    score(tl[0], 2257, "G", "T", total_s, no9)




infile.close()
