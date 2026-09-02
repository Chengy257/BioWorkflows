#!/usr/bin/python
#########################################################################
# File Name: getGroups.py
# Author: ChengYu
# Description: 
# Created Time: Thu 20 Jul 2023 08:51:33 PM CST
#########################################################################
import sys

infile = sys.argv[1]  ## input file sample_info.csv, list of sample IDs
DEG_dir = sys.argv[2]  ## DEGs directory
item = [] 
with open(infile,'r') as it:
	next(it)
	for line in it:
		line = line.strip().split(",")
		if line[-1] not in item:
			if line[-1] != "control":  ## 
				item.append(line[-1])
# print(item)
# item = ["A", "B", "C", "D", "E"]  ## test

for i in range(0,len(item),1):	
	for j in range(i+1,len(item),1):
		#print(i,j)
		for regu1 in ["UP", "DOWN", "ALL"]:
			for regu2 in ["UP", "DOWN", "ALL"]:
				# print(item[i] , regu1, item[j] , regu2, sep = " ")	
				out = "6.DEGcompare/combination/" + item[i] + "_" + regu1 + "_" + item[j] + "_" + regu2
				with open(out,'w') as f:					
					a = "".join([DEG_dir,item[i],"_vs_control_",regu1,".DEGs.txt"," ",item[i],"_",regu1])
					b = "".join([DEG_dir,item[j],"_vs_control_",regu2,".DEGs.txt"," ",item[j],"_",regu2])
					print(a,b,sep="\n",file = f)
					# print(item[i],item[j])
				f.close()
				