import sys,re,fileinput,os

Argument = []
Argument = sys.argv[1:] 

if (len(Argument)) < 3:
        print "Usage: Allelic_ratio_file gtf_file  Output_file" 
        sys.exit()

File = Argument[0]
refGene = Argument[1]
Output = Argument[2]

chr2gene = {}

ref3UTR = {}
ref5UTR = {}
refcodExons = {}
refIntron = {}
refInfo = {}

for line in fileinput.input([refGene]):
	if line.startswith("#"):
		continue

        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	strand = ""
	strand = rowlist[3]	

	if rowlist[2] not in chr2gene:
		chr2gene[rowlist[2]] = {}
		if rowlist[12] not in chr2gene[rowlist[2]]:
			chr2gene[rowlist[2]][rowlist[12]] = ""
			
	else:
		if rowlist[12] not in chr2gene[rowlist[2]]:
                        chr2gene[rowlist[2]][rowlist[12]] = ""
                    
        if rowlist[12] not in refInfo:
                refInfo[rowlist[12]] = {}
		ref5UTR[rowlist[12]] = {}
		refcodExons[rowlist[12]] = {}
		ref3UTR[rowlist[12]] = {}
		refIntron[rowlist[12]] = {}		

		refInfo[rowlist[12]][rowlist[1]] = []
                refInfo[rowlist[12]][rowlist[1]].append(int(rowlist[8])+2)
               
        
                ref5UTR[rowlist[12]][rowlist[1]] = []

		if strand == "+":
			ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[4])
			ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[6])
		else:
			ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[7])
                        ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[5])

		
                ref3UTR[rowlist[12]][rowlist[1]] = []

		if strand == "+":
               		ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[7])
                	ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[5])
		else:
			ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[4])
                        ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[6])


                refcodExons[rowlist[12]][rowlist[1]] = {}
		refIntron[rowlist[12]][rowlist[1]] = {}

		exons_start = []
		exons_end = []

                exons_start = map(int,rowlist[9].split(",")[:-1])
                exons_end = map(int,rowlist[10].split(",")[:-1])
		
		if strand == "+":
			exon_count = 0
                	for i in range(0,len(exons_start)):
				exon_count = exon_count+1
                        
				start = 0
				end = 0
				start = exons_start[i]
				end = exons_end[i]
                                refcodExons[rowlist[12]][rowlist[1]][exon_count] = []   
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(start)
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(end)
			
				if i != len(exons_start)-1:
					istart = 0
					iend = 0
					istart = int(exons_end[i])+1
					iend = int(exons_start[i+1])-1
					refIntron[rowlist[12]][rowlist[1]][exon_count] = []
					refIntron[rowlist[12]][rowlist[1]][exon_count].append(istart) 
					refIntron[rowlist[12]][rowlist[1]][exon_count].append(iend)
		else:
			exon_count = 0
			for i in range(0,len(exons_start)):
                                exon_count = exon_count+1

				start = 0
                                end = 0
                                start = exons_end[::-1][i]
                                end = exons_start[::-1][i]
                                refcodExons[rowlist[12]][rowlist[1]][exon_count] = []
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(end)
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(start)

				if i != len(exons_start)-1:
                                	istart = 0
                                	iend = 0
                                	istart = int(exons_start[::-1][i])+1
                                	iend = int(exons_end[::-1][i+1])-1
                                	refIntron[rowlist[12]][rowlist[1]][exon_count] = []
                                	refIntron[rowlist[12]][rowlist[1]][exon_count].append(iend)
                                	refIntron[rowlist[12]][rowlist[1]][exon_count].append(istart)


	else:

                refInfo[rowlist[12]][rowlist[1]] = []
                refInfo[rowlist[12]][rowlist[1]].append(int(rowlist[8])+2)
               
         
                ref5UTR[rowlist[12]][rowlist[1]] = []

                if strand == "+":
                        ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[4])
                        ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[6])
                else:
                        ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[7])
                        ref5UTR[rowlist[12]][rowlist[1]].append(rowlist[5])

                ref3UTR[rowlist[12]][rowlist[1]] = []

                if strand == "+":
                        ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[7])
                        ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[5])
                else:
                        ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[4])
                        ref3UTR[rowlist[12]][rowlist[1]].append(rowlist[6])

                refcodExons[rowlist[12]][rowlist[1]] = {}

                exons_start = []
                exons_end = []
		
		refIntron[rowlist[12]][rowlist[1]] = {}		

                exons_start = map(int,rowlist[9].split(",")[:-1])
                exons_end = map(int,rowlist[10].split(",")[:-1])

		if strand == "+":
                	exon_count = 0
                	for i in range(0,len(exons_start)):
                        	exon_count = exon_count+1
                                
                                start = 0
                                end = 0 
                                start = exons_start[i]
                                end = exons_end[i]
                                refcodExons[rowlist[12]][rowlist[1]][exon_count] = []    
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(start)
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(end)
				
				if i != len(exons_start)-1:
                                        istart = 0
                                        iend = 0
                                        istart = int(exons_end[i])+1
                                        iend = int(exons_start[i+1])-1
                                        refIntron[rowlist[12]][rowlist[1]][exon_count] = []
                                        refIntron[rowlist[12]][rowlist[1]][exon_count].append(istart)
                                        refIntron[rowlist[12]][rowlist[1]][exon_count].append(iend)

                else:
			exon_count = 0
                	for i in range(0,len(exons_start)):
                                exon_count = exon_count+1

                                start = 0
                                end = 0
                                start = exons_end[::-1][i]
                                end = exons_start[::-1][i]
                                refcodExons[rowlist[12]][rowlist[1]][exon_count] = []
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(end)
                                refcodExons[rowlist[12]][rowlist[1]][exon_count].append(start)

				if i != len(exons_start)-1:
                                        istart = 0
                                        iend = 0
                                        istart = int(exons_start[::-1][i])+1
                                        iend = int(exons_end[::-1][i+1])-1
                                        refIntron[rowlist[12]][rowlist[1]][exon_count] = []
                                        refIntron[rowlist[12]][rowlist[1]][exon_count].append(iend)
                                        refIntron[rowlist[12]][rowlist[1]][exon_count].append(istart) 


#print ref3UTR["Lmx1a"]
#print ref5UTR["Lmx1a"]
#print refcodExons["Lmx1a"]
#print refIntron["Lmx1a"]
#print refInfo["Lmx1a"]

#print ref3UTR["Adora1"]
#print ref5UTR["Adora1"]
#print refcodExons["Adora1"]
#print refIntron["Adora1"]
#print refInfo["Adora1"]

def annotate(gene,location):
	#flag = ""
	location = int(location)

	gene_info = {}
	gene_info[gene] = {}
	
	if gene not in ref5UTR.keys():
		sys.exit()

        for transcript in ref5UTR[gene]:
                if int(ref5UTR[gene][transcript][0]) <= location and location <= int(ref5UTR[gene][transcript][1]):
                        if not transcript in gene_info[gene]:
                                gene_info[gene][transcript] = ""
                                gene_info[gene][transcript] = "UTR5;"
                        else:
                                gene_info[gene][transcript] = gene_info[gene][transcript]+"UTR5;"

                        #flag = flag+"UTR5("+str(transcript)+");"

        for transcript in ref3UTR[gene]:
                if int(ref3UTR[gene][transcript][0]) <= location and location <= int(ref3UTR[gene][transcript][1]):
                        if not transcript in gene_info[gene]:
                                gene_info[gene][transcript] = ""
                                gene_info[gene][transcript] = "UTR3;"
                        else:
                                gene_info[gene][transcript] = gene_info[gene][transcript]+"UTR3;"

                        #flag = flag+"UTR3("+str(transcript)+");"

	
	for transcript in refcodExons[gene]:
		if transcript in gene_info[gene]:
			continue
	
		for exon in refcodExons[gene][transcript]:
			if refcodExons[gene][transcript][exon][0] <= location and location <=  refcodExons[gene][transcript][exon][1]:	
				if not transcript in gene_info[gene]:
                                	gene_info[gene][transcript] = ""
                                	gene_info[gene][transcript] = "EXON:"+str(exon)+";"
                        	else:
                                	gene_info[gene][transcript] = gene_info[gene][transcript]+"EXON:"+str(exon)+";"

				#flag = flag+"EXON:"+str(exon)+"("+str(transcript)+");"

	for transcript in refIntron[gene]:
		if transcript in gene_info[gene]:
			continue

		for intron in refIntron[gene][transcript]:
			if refIntron[gene][transcript][intron][0]<=location and location <= refIntron[gene][transcript][intron][1]:
				if not transcript in gene_info[gene]:
                                        gene_info[gene][transcript] = ""
                                        gene_info[gene][transcript] = "INTRON:"+str(intron)+";"
                                else:
                                        gene_info[gene][transcript] = gene_info[gene][transcript]+"INTRON:"+str(intron)+";"

				#flag =  flag+"INTRON:"+str(intron)+"("+str(transcript)+");"	
	
	return gene_info


import re 
outfile = open(Output,"w")

for line in fileinput.input([File]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	#print rowlist[0]
        if line.startswith("#") or rowlist[0] == "chrX":
                continue
        else:
		for gene in chr2gene[rowlist[0]]:
			result = {}
			result = annotate(gene,rowlist[1])
			for gene_name in result:
				if len(result[gene_name]) > 0:
					for trid in result[gene_name]:
						outfile.write(str(line.rstrip("\n"))+"\t"+str(gene_name)+"\t"+str(trid)+"\t"+str(result[gene_name][trid])+"\n")

outfile.close()			

