#!/usr/bin/env python
import re

#dictionary 't' is made
list1 = ['TTT','TTC','TTA','TTG','TCT','TCC','TCA','TCG','TAT','TAC','TAA','TAG','TGT','TGC','TGA','TGG','CTT','CTC','CTA','CTG','CCT','CCC','CCA','CCG','CAT','CAC','CAA','CAG','CGT','CGC','CGA','CGG','ATT','ATC','ATA','ATG','ACT','ACC','ACA','ACG','AAT','AAC','AAA','AAG','AGT','AGC','AGA','AGG','GTT','GTC','GTA','GTG','GCT','GCC','GCA','GCG','GAT','GAC','GAA','GAG','GGT','GGC','GGA','GGG']
list2 = ['F','F','L','L','S','S','S','S','Y','Y','*','*','C','C','*','W','L','L','L','L','P','P','P','P','H','H','Q','Q','R','R','R','R','I','I','I','M','T','T','T','T','N','N','K','K','S','S','R','R','V','V','V','V','A','A','A','A','D','D','E','E','G','G','G','G']
triNuc = zip(list1,list2)
t = dict(triNuc)

#reverse dictionary is made
invTriNuc = zip(list2,list1)
iD = dict(invTriNuc)

#file is opened and is read into a lines list
with open("lab7.fa") as f:
	f.readline()
	line = f.readlines()
lines = []
for i in range(0,len(line)):
	lines.append(line[i].rstrip())
lines = "".join(lines)

#dna lists 1,2,3 are created
dna1 = []
dna2 = []
dna3 = []

#dna lists 1,2,3 are populated with codons
for pos in range(0,len(lines),3):
	codon = lines[pos:pos+3]
	dna1.append(codon)

for pos in range(1,len(lines),3):
        codon = lines[pos:pos+3]
        dna2.append(codon)

for pos in range(2,len(lines),3):
        codon = lines[pos:pos+3]
        dna3.append(codon)

#protein lists 1,2,3 are created
protein1 = []
protein2 = []
protein3 = []

#protein lists 1,2,3 are populated by changing codons from dna list to proteins
for x in range(0,len(dna1)):
	for y in range(0,len(t)):
		if list1[y] ==  dna1[x]:
			protein1.append(t.get(dna1[x]))
		else:
			continue		
protein1 = "".join(protein1)

for x in range(0,len(dna2)):
        for y in range(0,len(t)):
                if list1[y] ==  dna2[x]:
                        protein2.append(t.get(dna2[x]))
                else:
                        continue
protein2 = "".join(protein2)

for x in range(0,len(dna3)):
        for y in range(0,len(t)):
                if list1[y] ==  dna3[x]:
                        protein3.append(t.get(dna3[x]))
                else:
                        continue
protein3 = "".join(protein3)

#genes are found and are set to p1,p2,p3			
p1 = re.findall(r"M\w+\*",protein1)
p2 = re.findall(r"M\w+\*",protein2)
p3 = re.findall(r"M\w+\*",protein3)

#use iterator to find start and stop of codons for first 3 genes
reF1 = re.finditer(r"M\w+\*",protein1)
g1 = []

for i in reF1:
	start = i.start()
	end = i.end()
	g1.append(lines[start*3:end*3])

reF2 = re.finditer(r"M\w+\*",protein2)
g2 = []

for i in reF2:
        start = i.start()
        end = i.end()
        g2.append(lines[start*3:end*3])

reF3 = re.finditer(r"M\w+\*",protein3)
g3 = []

for i in reF3:
        start = i.start()
        end = i.end()
        g3.append(lines[start*3:end*3])

#end of first 3 frames

#reverse the dna sequence
temp = lines[::-1]
revLines = []
for i in range(0,len(temp)):
	if(temp[i] == "G"):
		revLines.append("C")
	if(temp[i] == "C"):
		revLines.append("G")
	if(temp[i] == "A"):
		revLines.append("T")
	if(temp[i] == "T"):
		revLines.append("A")
revLines = "".join(revLines)

#dna lists 4,5,6 are created
dna4 = []
dna5 = []
dna6 = []

#dna lists 4,5,6 are populated with codons
for pos in range(0,len(revLines),3):
        codon = revLines[pos:pos+3]
        dna4.append(codon)

for pos in range(1,len(revLines),3):
        codon = revLines[pos:pos+3]
        dna5.append(codon)

for pos in range(2,len(revLines),3):
        codon = revLines[pos:pos+3]
        dna6.append(codon)

#protein lists 4,5,6 are created
protein4 = []
protein5 = []
protein6 = []

#protein lists 4,5,6 are populated by changing codons from dna list to proteins
for x in range(0,len(dna4)):
        for y in range(0,len(t)):
                if list1[y] ==  dna4[x]:
                        protein4.append(t.get(dna4[x]))
                else:
                        continue
protein4 = "".join(protein4)

for x in range(0,len(dna5)):
        for y in range(0,len(t)):
                if list1[y] ==  dna5[x]:
                        protein5.append(t.get(dna5[x]))
                else:
                        continue
protein5 = "".join(protein5)

for x in range(0,len(dna6)):
        for y in range(0,len(t)):
                if list1[y] ==  dna6[x]:
                        protein6.append(t.get(dna6[x]))
                else:
                        continue
protein6 = "".join(protein6)

#genes are found and are set to p4,p5,p6                      
p4 = re.findall(r"M\w+\*",protein4)
p5 = re.findall(r"M\w+\*",protein5)
p6 = re.findall(r"M\w+\*",protein6)

#use iterator to find start and stop of codons for last 3 genes
reF4 = re.finditer(r"M\w+\*",protein4)
g4 = []

for i in reF4:
        start = i.start()
        end = i.end()
        g4.append(revLines[start*3:end*3])

reF5 = re.finditer(r"M\w+\*",protein5)
g5 = []

for i in reF5:
        start = i.start()
        end = i.end()
        g5.append(revLines[start*3:end*3])

reF6 = re.finditer(r"M\w+\*",protein6)
g6 = []

for i in reF6:
        start = i.start()
        end = i.end()
        g6.append(revLines[start*3:end*3])

#prints contents to output file
with open("output.txt","w") as f:
	f.write(">Gene1" + "\n" + g1[0] + "\n" + ">Protein1" + "\n" + p1[0] + "\n" + "\n")
	f.write(">Gene2" + "\n" + g1[1] + "\n" + ">Protein2" + "\n" + p1[1] + "\n" + "\n")
	f.write(">Gene3" + "\n" + g2[0] + "\n" + ">Protein3" + "\n" + p2[0] + "\n" + "\n")
	f.write(">Gene4" + "\n" + g2[1] + "\n" + ">Protein4" + "\n" + p2[1] + "\n" + "\n")
	f.write(">Gene5" + "\n" + g3[0] + "\n" + ">Protein5" + "\n" + p3[0] + "\n" + "\n")
	f.write(">Gene6" + "\n" + g3[1] + "\n" + ">Protein6" + "\n" + p3[1] + "\n" + "\n")
	f.write(">Gene7" + "\n" + g4[0] + "\n" + ">Protein7" + "\n" + p4[0] + "\n" + "\n")
	f.write(">Gene8" + "\n" + g4[1] + "\n" + ">Protein8" + "\n" + p4[1] + "\n" + "\n")
	f.write(">Gene9" + "\n" + g5[0] + "\n" + ">Protein9" + "\n" + p5[0] + "\n" + "\n")
        f.write(">Gene10" + "\n" + g5[1] + "\n" + ">Protein10" + "\n" + p5[1] + "\n" + "\n")
	f.write(">Gene11" + "\n" + g5[2] + "\n" + ">Protein11" + "\n" + p5[2] + "\n" + "\n")
        f.write(">Gene12" + "\n" + g5[3] + "\n" + ">Protein12" + "\n" + p5[3] + "\n" + "\n")
	f.write(">Gene13" + "\n" + g5[4] + "\n" + ">Protein13" + "\n" + p5[4] + "\n" + "\n")
        f.write(">Gene14" + "\n" + g5[5] + "\n" + ">Protein14" + "\n" + p5[5] + "\n" + "\n")
	f.write(">Gene15" + "\n" + g6[0] + "\n" + ">Protein15" + "\n" + p6[0] + "\n" + "\n")
        f.write(">Gene16" + "\n" + g6[1] + "\n" + ">Protein16" + "\n" + p6[1] + "\n" + "\n")
	f.close()


