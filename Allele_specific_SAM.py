import sys,re,itertools,os

#This script is used to detect allele of origin for a RNA-seq read
# RNA-seq data was aligned to two diffferent strains comprising of parental genome of a cross 
#This program assigns aligned paired-end reads to a particular based on mapping quality and mismatches
#This program compares mapped Two sam files for the two strains and then produce new sam files for each of the strain
#It only seleects reads with high mapping quality that are aligned in conrcondance with the library design

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 3:
        print "Usage: Input_sam_file_1 Input_sam_file_2 Sam_header_info_1 Output_sam_file_1 Output_sam_file_2" 
        sys.exit()

Filepath_1 = Argument[0]
Filepath_2 = Argument[1]

Header = Argument[2]

Bam_Header = ""

for a in Header.split(","):
	Bam_Header = Bam_Header+"\t"+str(a)

print Bam_Header

output_1 = open(Argument[3],"w")
output_2 = open(Argument[4],"w")

#discard = open(str(Argument[3])+"_discarded","w")

def numeric_compare(x, y):
        x1 = int(x)
        y1 = int(y)
        return x1 - y1

def filter_check(list1,list2):
	read1 = []
	read2 = []

	Chromosome = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX"]
	Flags = ["99","147","83","163","73","133","69","137","89","165","101","153"] #Samtools flags

	read1 = list1
	read2 = list2

	print read1[0], read2[0]

	if read1[0] != read2[0]:
		print read1[0],read2[0]
		print "Reads not in synchrony within sam"
		sys.exit()

	if read1[2] in Chromosome or read2[2] in Chromosome:
		if (read1[1] in Flags or read2[1] in Flags) and int(read1[8]) <= 300000:
			NH1=0
			nM1=0
			AS1=0
	
			NH1 = read1[11].lstrip("NH:i:")	
			if int(NH1) == 1:
				nM1 = int(read1[13].lstrip("nM:i:"))
				AS1 = int(read1[14].lstrip("AS:i:"))
				return nM1,AS1
							
	return 100,0
			
def loc_comparison(pos1,pos2):
	loc1 = 0
	loc2 = 0

	loc1 = int(pos1)
	loc2 = int(pos2)

	if loc1 == loc2:
		return "True"

	if loc1 < loc2 and loc1+50 >= loc2:
		return "True"

	if loc1 > loc2 and loc2+50 >= loc1:
		return "True"

	return "False"

with open(Filepath_1) as A, open(Filepath_2) as B:
	while True:
        	try:
			#output_1.flush()
			#output_2.flush()

            		line1, line2, line_1, line_2 = next(A), next(A), next(B), next(B)

			if line1.startswith("@") and line2.startswith("@"):
				output_1.write(str(line1))
                		output_1.write(str(line2))
				
        		if line_1.startswith("@") and line_2.startswith("@"):
        			output_2.write(str(line_1))
                		output_2.write(str(line_2))

				if line_2.startswith("@PG"):
                                	output_1.write("@RG"+str(Bam_Header)+"\n")
                                	output_2.write("@RG"+str(Bam_Header)+"\n")
                		
				continue

			
        		rowlist1 = []
			rowlist2 = []

			rowlist_1 = []
                        rowlist_2 = []

			mismatch1 = 100
			mismatch2 = 100

			alignscore1 = 0
			alignscore2 = 0

        		rowlist1 = (line1.rstrip("\n")).split('\t')
			rowlist2 = (line2.rstrip("\n")).split('\t')

			rowlist_1 = (line_1.rstrip("\n")).split('\t')
                        rowlist_2 = (line_2.rstrip("\n")).split('\t')

			if rowlist1[0] != rowlist_1[0]:
				print "Reads not in synchrony between sam files"
  				sys.exit()

			if rowlist1[2] != rowlist_1[2] and rowlist1[2] != rowlist_2[2] and rowlist2[2] != rowlist_1[2] and rowlist2[2] != rowlist_2[2]:
				continue

			mismatch1,alignscore1 = filter_check(rowlist1,rowlist2)
			mismatch2,alignscore2 = filter_check(rowlist_1,rowlist_2)
			
			#print rowlist1[0], rowlist2[0], rowlist_1[0], rowlist_2[0]
			
			if mismatch1 == mismatch2:
				if alignscore1 == alignscore2:
					output_1.write(str(line1)+str(line2))
					output_2.write(str(line_1)+str(line_2))
					continue
				elif alignscore1 > alignscore2:
 	                        	output_1.write(str(line1)+str(line2))
        	                        continue
				else:
					output_2.write(str(line_1)+str(line_2))
                                	continue

			if mismatch1 < mismatch2:
				output_1.write(str(line1)+str(line2))
				continue

			if mismatch1 > mismatch2:
				output_2.write(str(line_1)+str(line_2))
				continue

			if alignscore1 > alignscore2:
				output_1.write(str(line1)+str(line2))
                                continue
						
			if  alignscore1 <  alignscore2:
                                output_2.write(str(line_1)+str(line_2))
                                continue
		
		except StopIteration:
        		break

output_1.close()
output_2.close()


