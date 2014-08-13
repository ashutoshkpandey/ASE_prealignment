import sys,re,fileinput,os
from decimal import Decimal, getcontext

# This program calculates the p-value for allelic ratio
# Input files are heterozygous SNPs, Allelic_ratio file

Argument = []
Argument = sys.argv[1:]

if (len(Argument)) < 3:
        print "Usage: Input_allelic_ratio_file SNP_file Output_file"
        sys.exit()

A = Decimal('15')

def mygamma(z):
    """
        The constant SQRT2PI is defined as sqrt(2.0 * PI);
        For speed the constant is already defined in decimal
        form.  However, if you wish to ensure that you achieve
        maximum precision on your own machine, you can calculate
        it yourself using (sqrt(atan(1.0) * 8.0))
    """
    #const long double SQRT2PI = sqrtl(atanl(1.0) * 8.0);
    SQRT2PI = Decimal('2.5066282746310005024157652848110452530069867406099383')
    f = Decimal('1')
    sum_v = SQRT2PI
    sc = getcontext().power(z+A,z+Decimal('0.5'))
    sc *= Decimal(Decimal('-1') * (z+A)).exp()
    sc /= z
    for k in range(1,15):
        z+=Decimal('1')
        ck = getcontext().power(A - Decimal(k) , Decimal(k) - Decimal('0.5'))
        ck *= Decimal(A -Decimal(k)).exp()
        ck /= f
        sum_v += (ck / z)
        f *= (Decimal('-1') * k)
    return sum_v * sc

def approx_gamma(z):
    RECIP_E = Decimal('0.36787944117144232159552377016147')
    TWOPI = Decimal('6.283185307179586476925286766559')
    d = Decimal('1.0') / (Decimal('10.0') * z)
    d = Decimal('1.0') / ((Decimal('12') * z) - d)
    d = (d + z) * RECIP_E
    d = getcontext().power(d,z)
    d *= Decimal.sqrt(TWOPI / z)
    return d

def igf(s, z):
    if z < Decimal('0'):
        return Decimal('0')
    sc = Decimal('1') / s
    sc *= getcontext().power(z,s)
    sc *= Decimal(-z).exp()
    sum_v = Decimal('1')
    nom = Decimal('1')
    denom = Decimal('1')
    for i in range(0,200):
        nom *= z
        s+=Decimal('1')
        denom *= s
        sum_v += (nom / denom)
    return sum_v * sc

def chisqr(dof, cv):
    if cv < Decimal('0') or dof < Decimal('1'):
        return Decimal('0')
    k = dof * Decimal('0.5')
    x = cv * Decimal('0.5')
    if dof == Decimal('2'):
        print Decimal(Decimal('-1') * x).exp()
        return

    pvalue = igf(k,x)
   
    if pvalue.is_nan() or pvalue.is_infinite() or pvalue <= 1e-8:
        return 1e-14
    pvalue /= mygamma(k)   
    #print dof,cv,Decimal('1') - pvalue
    return Decimal('1') - pvalue


File_Allelic_ratio = Argument[0]
File_SNP = Argument[1]
output = open(Argument[2],"w")
output_2 = open(str(Argument[2])+"_20reads","w")

Allelic_ratio = {}
Chisqr_saved = {}
Chisqr_pvalue = {}

for line in fileinput.input([File_Allelic_ratio]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	
	if rowlist[0]+"\t"+rowlist[1] not in Allelic_ratio:
		Allelic_ratio[rowlist[0]+"\t"+rowlist[1]] = ""
		Allelic_ratio[rowlist[0]+"\t"+rowlist[1]] = line

	if rowlist[0]+"\t"+rowlist[1] not in Chisqr_pvalue:
		cv = 0.0
		sum = 0
		o1 = 0
		o2 = 0
		e1 = 0
		e2 = 0

		o1 = float(rowlist[4])+float(rowlist[7])
		o2 = float(rowlist[6])+float(rowlist[5])
		sum = o1+o2

		if sum < 1.00:
			continue

		e1 = float(sum)/2.0
		e2 = float(sum)/2.0

		chisqr_value = 0.0
		chisqr_value  = (((o1-e1)*(o1-e1))/e1)+(((o2-e2)*(o2-e2))/e2)
		
		if chisqr_value <= 0.00001:
			chisqr_value = .00001
			
		Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]] = 1.00,0.5,0.0

		if chisqr_value not in Chisqr_saved:
			p_value = 1.00
			p_value = chisqr(Decimal('1'), Decimal(chisqr_value))
			Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]] = (p_value,o1/sum,sum)
			Chisqr_saved[chisqr_value] = p_value

		else:
			Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]] = (Chisqr_saved[chisqr_value],o1/sum,sum)

#print Chisqr_pvalue['chr1'+"\t"+'3228500'] 

for line in fileinput.input([File_SNP]):
        rowlist = []
        rowlist = (line.rstrip("\n")).split('\t')
	if line.startswith("#") or rowlist[0] == "chrX" or rowlist[7].startswith("INDEL;"):
		continue
	
	if rowlist[0]+"\t"+rowlist[1] in Chisqr_pvalue:
		if int(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][2]) >= 20:
			output_2.write(str((Allelic_ratio[rowlist[0]+"\t"+rowlist[1]]).rstrip("\n"))+"\t"+str(line.rstrip("\n"))+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][0])+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][1])+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][2])+"\n")
		output.write(str((Allelic_ratio[rowlist[0]+"\t"+rowlist[1]]).rstrip("\n"))+"\t"+str(line.rstrip("\n"))+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][0])+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][1])+"\t"+str(Chisqr_pvalue[rowlist[0]+"\t"+rowlist[1]][2])+"\n")
			
output.close()
output_2.close()

