#!/usr/bin/env python

import os
import sys

def overlaps(t1, t2):
	return (t1[0] <= t2[1] and t1[1] >= t2[0])

def getNonOverlap(l):
	ret = [l[0]]
	
	for p in l:
		if overlaps(p,ret[-1]):
			ret[-1] = (min(p[0],ret[-1][0]), max(ret[-1][1], p[1]))
		else:
			ret.append(p)
	
	return ret
		

if __name__ == "__main__":
	
	chr_sizes={}
	# Open the chromosome info
	#chr_f = file("chromInfo.txt",'r')
	#for l in chr_f:
	#	(ch, sz, loc) = l.split('\t',2)
	#	chr_sizes[ch] = int(sz)
		
	gene_f = file("knownGene.txt", 'r')
	
	exon_list = {}
	
	for l in gene_f:
		data = l.split('\t')
		is_fwd = (data[2] == '+')
		ch = data[1]
		if ch not in exon_list:
			exon_list[ch] = set()
			
		#if is_fwd:
		st_idx = 8
		en_idx = 9
		#else:
		#	st_idx = 9
		#	en_idx = 8
		
		st_list = [int(d) for d in data[st_idx].strip(',').split(',')]
		en_list = [int(d) for d in data[en_idx].strip(',').split(',')]
		
		#if not is_fwd:
		#	st_list = [chr_sizes[ch] - i for i in st_list]
		#	en_list = [chr_sizes[ch] - i for i in en_list]
			
		exon_list[ch].update(set(zip(st_list, en_list)))
		
	k = exon_list.keys()
	k.sort()
		
	for c in k:
		l = list(exon_list[c])
		l.sort()
		l = getNonOverlap(l)
		for (t, e) in l:
			print "%s\t%d\t%d" % (c, t, e)
			
	
