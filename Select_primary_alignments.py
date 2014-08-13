import os, sys

Argument = []
Argument = sys.argv[1:]

if (len(Argument)) < 2:
        print "Usage:Input_Sam Output_Sam"
        sys.exit()


output = open(Argument[1],"w")
input = open(Argument[0])

output.write("@HD"+"\t"+"VN:1.0"+"\t"+"SO:unsorted\n")

Counter = {}

read = ""
lastread = ""

def smallest_insert(IS_dict):
	IS = {}
	IS = IS_dict

	IS_sort = []
	IS_sort = sorted(IS.keys(),key=int)
		
	return IS_sort[0],IS_sort[0].lstrip("-")


for line in input:
        if line.startswith("@"):
        	output.write(str(line))
                continue
		
        v = []
        v = line.split("\t")
	
	read = v[0]

	if read == lastread:
		lastread = read
		if v[0] not in Counter:
			Counter[v[0]] = {}
			Counter[v[0]][v[1]] = line 
		else:
			Counter[v[0]][v[1]] = line

	else:
		if lastread == "":
			pass
		else:
			primary = []
			primary = sorted(Counter[lastread].keys(), key=int)[0:2]
	
			for i in primary:
				output.write(str(Counter[lastread][i]))

		lastread = read
		Counter = {}
	
		if v[0] not in Counter:
                        Counter[v[0]] = {}
                        Counter[v[0]][v[1]] = line
                else:
                        Counter[v[0]][v[1]] = line

end = []
end = sorted(Counter[lastread].keys(), key=int)[0:2]

for i in end:
	output.write(str(Counter[lastread][i]))

output.close()
