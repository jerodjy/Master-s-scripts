#This script merges the probabilities of a pair of ancestries from the vit file. 
#This approach is useful if two different ancestries should be treated as the same source after being called reliably as two separate ancestries by RFMix.

#Oct 1st 2018
#Juan Esteban Rodriguez

import sys,os
from itertools import izip
import numpy as np

vit_fwb=sys.argv[1] ## Prefix from input Viterbi and ForwardBackward files. (The script assumes this file syntax: MEX_individuals_chr22_shapeout.0.ForwardBackward.txt. MEX_individuals string must be provided to the script.)
merged_anc=sys.argv[2] ## Two ancestry numbers to merge (in increasing order)
anc_string=sys.argv[3] ## Ancestry strings. Only required to avoid mistakes in the final ancestry order.

#The output of the file can be easily provided to scripts plotting karyograms

#Example of how to run this script
#Here, European and North African ancestry are merged into one ancestry:
#python merge_ancestries_vit.py MEX_individuals '2,5' 'AFR,EUR,NAT,ASN,NAFR'
#Final ancestry order (merged ancestry is always put at the end):
#AFR,NAT,ASN,EUR-NAFR

#Read other parameters
ancs=anc_string.split(",")
merges=map(int,merged_anc.split(","))	#Str to int

#Read both Viterbi and ForwardBackward files (a and b) and open output files (c and d)
for chr in range(1,23):		#Read per chromosome
	chr=str(chr)
	with open(vit_fwb+"_chr"+chr+"_shapeout.0.ForwardBackward.txt", "r") as a, open(vit_fwb+"_chr"+chr+"_shapeout.0.Viterbi.txt","r") as b, open (vit_fwb+"_Merged_chr"+chr+"_shapeout.0.Viterbi.txt","w") as c, open(vit_fwb+"_Merged_chr"+chr+"_shapeout.0.ForwardBackward.txt","w") as d:
		#Read both input files simultaneously
		for line1,line2 in izip(a,b):
			line1=line1.rstrip()
			line2=line2.rstrip()
			fwb_line=map(float,line1.split(" "))
			vit_line=map(int,line2.split(" "))
			fwb_out_line=[]
			vit_out_line=[]

			#Read ancestry probabilities for each individual
			for i in range(0,len(fwb_line),len(ancs)):
				temp=fwb_line[i:(i+len(ancs))]	#Get all ancestry probabilities from a specific SNP and individual
	
				new_anc_prob = temp.pop(merges[1]-1) + temp.pop(merges[0]-1)	#First and second number in merged_anc: merges[]
				temp.append(new_anc_prob)
		
				vit_out_line.append(np.argmax(temp)+1)	#To get the index from the largest number. Will be printed in the new Viterbi file.
	
                	        temp=map(str,temp)                      		
				fwb_out_line.extend(temp)
	
			c.write(' '.join(map(str,vit_out_line))+'\n')	#Print modified line for the Viterbi file
			d.write(' '.join(fwb_out_line)+'\n')		#Print modified line for the ForwardBackward file

#Confirm the new ancestry order to the user
ancs.append(ancs.pop(merges[1]-1)+'-'+ancs.pop(merges[0]-1))
temp=','.join(ancs)
print "The new ancestry order is: "+temp+"\n"
