myseq = "potamo_transcripts.fasta"
myheaders = "contig_diff_n_s_2.txt"
myout = "diff_n_s.fasta"

SEQ = open(myseq,"r")
HEAD = open(myheaders,"r")
OUT = open(myout,"w")

mydict = {}

for line in SEQ:
    line = line.rstrip()
    if ">" in line:
        header = line.lstrip(">")
        mydict[header] = ""
    else:
        mydict[header] += line

SEQ.close()

for line in HEAD:
    header = line.rstrip()
    OUT.write(">"+header+"\n")
    OUT.write(mydict[header]+"\n")

HEAD.close()
OUT.close()
