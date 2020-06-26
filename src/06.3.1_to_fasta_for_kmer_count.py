#! /bin/env python3
import os
import sys
#out_pos = open("")
#CB = os.getenv("CODEBASE")
data_dir = <your path to data>
subdir = os.path.join(data_dir, "results_analysis", "/02_predict_meth/02.2_test_on_other_species/screen/")


def save_to_fa(path_table, subdir, SP):
	print(path_table)
	if os.path.exists((path_table)):
		out_pos = open(subdir+SP+"/sequences/"+SP+"_high.fa", "w")
		out_neg = open(subdir+SP+"/sequences/"+SP+"_low.fa", "w")
		with open(path_table) as f:
			for line in f.readlines():
				line_array = line[:-1].split("\t")
				if line_array[-1]=="1":
					out_pos.write(">"+line_array[0]+"\n")
					out_pos.write(line_array[1]+"\n")
				elif line_array[-1]=="-1":
					out_neg.write(">"+line_array[0]+"\n")
					out_neg.write(line_array[1]+"\n")
 
		out_pos.close()
		out_neg.close()
		return 1
	else:
		return 0

with open(os.getenv("CODEBASE")+"/compEpi/meta/species_list.txt", "r") as sp_list:
	species = sp_list.readlines()

for sp in species:

	print(save_to_fa(subdir+sp[:-1]+"/sequences/sequence_table_"+sp[:-1]+".tsv", subdir, sp[:-1]))
