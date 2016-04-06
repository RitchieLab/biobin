#!/usr/bin/env python

import sys
import argparse
import glob

if __name__ == "__main__":
	
	parser=argparse.ArgumentParser()
	
	parser.add_argument("-p", "--prefix", help="Prefix of BioBin files (do not include the '-')", required=True)
	parser.add_argument("-c", "--control-value", help="Phenotype control value", default="0")
	parser.add_argument("-a", "--all-control", help="Use all non-missing subjects as controls", action="store_true", default=False)
	parser.add_argument("-s", "--sep", help="Separator character in bins file(s)", default=",")
	parser.add_argument("-o", "--output-sep", help="Separator character for output", default="\t")
	
	args=parser.parse_args()
	
	pv_names = set([])
	header_printed = False
	
	for fn in glob.glob(args.prefix + "-*bins.csv"):
	
		pheno_name = fn[len(args.prefix)+1:-9]
		f = file(fn, 'r')
		
		n_case=0
		n_control=0
		
		bin_names=f.readline().strip().split(args.sep)[2:]		
		
		bin_vars=f.readline().strip().split(args.sep)[2:]
		bin_loci=f.readline().strip().split(args.sep)[2:]
		
		ctl_vars=f.readline().strip().split(args.sep)[2:]
		case_vars=f.readline().strip().split(args.sep)[2:]
		
		ctl_cap=f.readline().strip().split(args.sep)[2:]
		case_cap=f.readline().strip().split(args.sep)[2:]
		
		# OK, now the header is read, let's read the p-values
		
		pv_dict = {}
		for l in f:
			pv_arr = l.strip().split(args.sep)
			if pv_arr[0].endswith("p-value"):
				pv_test = pv_arr[0].rsplit(" ", 1)
				if not header_printed:
					pv_names.add(pv_test[0])
				pv_dict[pv_test[0]] = pv_arr[2:]
			else:
				n_case += (pv_arr[1] != "nan" and (not args.all_control) and pv_arr[1] != args.control_value)
				n_control += (pv_arr[1] != "nan" and (args.all_control or pv_arr[1] == args.control_value))
	
		for l in f:
			ptid, status, junk = l.split(args.sep, 2)
		
			n_case += (status != "nan" and (not args.all_control) and status != args.control_value)
			n_control += (status != "nan" and (args.all_control or status == args.control_value))
		

		f.close()
		
		n_case=str(n_case)
		n_control=str(n_control)
		
		# OK, now print everything for this file
		# First print the header
		if not header_printed:
			header_printed = True
			print args.output_sep.join(
				[
					"Outcome",
					"Bin",
					"Num_Cases",
					"Num_Controls",
					"Num_Loci",
					"Num_Vars",
					"Case_Vars",
					"Control_Vars",
					"Case_Capacity",
					"Control_Capacity"
				] + [n.replace('-','_') for n in pv_names]
			)
		
		# Now, enumerate all the lists we have - they should have the same # of elements
		for i in xrange(len(bin_names)):
			print args.output_sep.join(
				[
					pheno_name,
					bin_names[i],
					n_case,
					n_control,
					bin_loci[i],
					bin_vars[i],
					case_vars[i],
					ctl_vars[i],
					case_cap[i],
					ctl_cap[i]
				] + [pv_dict[n][i] for n in pv_names]
			)
		
	
