import sys,re,fileinput,os

# This program calculates allelic ratio for list of positions 
# Input files are heterozygous SNPs, CNV regions, Pileup files
# Two pileup files are read, each for the reference genome and enhanced genome

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 5:
        print "Usage: Input_pileup_file_1 Input_pileup_file_2  SNP_file CNV_regions Output_file" 
        sys.exit()

File_Pileup_1 = Argument[0]
File_Pileup_2 = Argument[1]
File_SNP = Argument[2]
File_CNV = Argument[3]

output = open(str(Argument[4]),"w")

def numeric_compare(x, y):
        x1 = int(x)
        y1 = int(y)
        return x1 - y1

def CNV_span(chr,loc):
	chromosome = ""
	location = 0

	chromosome = chr
	location = int(loc)

	for cnv in CNVs[chromosome]:
		if location >= int(cnv[0]) and location <= int(cnv[1]):
			return True
		else:
			return False

CNVs = {}

for line in fileinput.input([File_CNV]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#"):
                continue
        else:
                if rowlist[0] not in CNVs:
			CNVs[rowlist[0]] = []
			CNVs[rowlist[0]].append(rowlist[1]+"\t"+rowlist[2])
		else:
			CNVs[rowlist[0]].append(rowlist[1]+"\t"+rowlist[2])

print "CNVs read"

SNPs = {}

for line in fileinput.input([File_SNP]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#") or rowlist[0] == "chrX" or rowlist[7].startswith("INDEL;"):
                continue
        else:
		if rowlist[0]+"\t"+rowlist[1] not in SNPs:
			#print rowlist[0],rowlist[1]
			if not CNV_span(rowlist[0],rowlist[1]):
				SNPs[rowlist[0]+"\t"+rowlist[1]] = []
				SNPs[rowlist[0]+"\t"+rowlist[1]].append(rowlist[3])
				SNPs[rowlist[0]+"\t"+rowlist[1]].append(rowlist[4])

print "SNPs read"

nucleotides = ["A","T","C","G","a","t","c","g"]
complement = {'a':'t','g':'c','t':'a','c':'g'}

Allelic_count_B6 = {}
Allelic_count_D2 = {}

B6_info = {}
D2_info = {}

for line in fileinput.input([File_Pileup_1]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#"):
                continue
        else:
		if rowlist[0]+"\t"+rowlist[1] in SNPs:
			if rowlist[2] == SNPs[rowlist[0]+"\t"+rowlist[1]][0]:
				count_1 = 0
				count_2 = 0

				for i in rowlist[8]:
					if i == "." or i == ",":
						count_1 = count_1 + 1
					else:
						if i in nucleotides:
							if i.isupper():
								if i.lower() == SNPs[rowlist[0]+"\t"+rowlist[1]][1].lower():
									count_2 = count_2 + 1
							else:
								if complement[i.lower()] == SNPs[rowlist[0]+"\t"+rowlist[1]][1].lower():
									count_2 = count_2 + 1

				if rowlist[0]+"\t"+rowlist[1] not in Allelic_count_B6:
					Allelic_count_B6[rowlist[0]+"\t"+rowlist[1]] = []
					Allelic_count_B6[rowlist[0]+"\t"+rowlist[1]].append(count_1)
					Allelic_count_B6[rowlist[0]+"\t"+rowlist[1]].append(count_2) 
					B6_info[rowlist[0]+"\t"+rowlist[1]] = line.rstrip("\n")
			
print "Pileup one done"
												
for line in fileinput.input([File_Pileup_2]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')

        if line.startswith("#"):
                continue
        else:
                if rowlist[0]+"\t"+rowlist[1] in SNPs:
                        if rowlist[2] == SNPs[rowlist[0]+"\t"+rowlist[1]][1]:
                        	count_1 = 0
                                count_2 = 0

                        	for i in rowlist[8]:
                                	if i == "." or i == ",":
                                        	count_1 = count_1 + 1
					else:
						if i in nucleotides:
							if i.isupper():
								if i.lower() == SNPs[rowlist[0]+"\t"+rowlist[1]][0].lower():
									count_2 = count_2 + 1
							else:
								if complement[i.lower()] == SNPs[rowlist[0]+"\t"+rowlist[1]][0].lower():
                                                                        count_2 = count_2 + 1

                               	if rowlist[0]+"\t"+rowlist[1] not in Allelic_count_D2:
                                	Allelic_count_D2[rowlist[0]+"\t"+rowlist[1]] = []
                                        Allelic_count_D2[rowlist[0]+"\t"+rowlist[1]].append(count_1)
					Allelic_count_D2[rowlist[0]+"\t"+rowlist[1]].append(count_2)
					D2_info[rowlist[0]+"\t"+rowlist[1]] = line.rstrip("\n")

print "Pileup two done"
 
for snp in SNPs:
	if snp in Allelic_count_B6 and snp in Allelic_count_D2:
		output.write(str(snp)+"\t"+str(SNPs[snp][0])+"\t"+str(SNPs[snp][1])+"\t"+str(Allelic_count_B6[snp][0])+"\t"+str(Allelic_count_B6[snp][1])+"\t"+str(Allelic_count_D2[snp][0])+"\t"+str(Allelic_count_D2[snp][1])+"\t"+str(B6_info[snp])+"\t"+str(D2_info[snp])+"\n")
		continue

	if snp not in Allelic_count_B6 and snp not in Allelic_count_D2:
                continue

	if snp in Allelic_count_B6 and snp not in Allelic_count_D2:
                output.write(str(snp)+"\t"+str(SNPs[snp][0])+"\t"+str(SNPs[snp][1])+"\t"+str(Allelic_count_B6[snp][0])+"\t"+str(Allelic_count_B6[snp][1])+"\t0\t0\t"+str(B6_info[snp])+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                continue

	if snp not in Allelic_count_B6 and snp in Allelic_count_D2:
        	output.write(str(snp)+"\t"+str(SNPs[snp][0])+"\t"+str(SNPs[snp][1])+"\t0\t0\t"+str(Allelic_count_D2[snp][0])+"\t"+str(Allelic_count_D2[snp][1])+"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+str(D2_info[snp])+"\n")
                continue
				
			
output.close()
