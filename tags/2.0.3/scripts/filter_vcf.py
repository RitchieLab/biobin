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
	
	vcf_fn = sys.argv[1]
	exon_fn = sys.argv[2]
	
	e_f = file(exon_fn, 'r')
	chr_exon = {}
	for l in e_f:
		(ch, s, e) = l.split('\t')
		if ch.lower().startswith('chr'):
			ch = ch[3:]
		
		if ch not in chr_exon:		
			chr_exon[ch] = []
			
		chr_exon[ch].append((int(s), int(e)))
		
	for (k, v) in chr_exon.iteritems():
		v.sort()
		chr_exon[k] = getNonOverlap(v)
		
	v_f = file(vcf_fn, 'r')
	
	curr_ch=''
	curr_loc=0
	curr_list=[]
	for l in v_f:
		if l[0] == '#':
			print l,
		else:
			(ch,pos,o) = l.split('\t',2)
			if ch.lower().startswith('chr'):
				ch = ch[3:]
			
			if ch != curr_ch:
				curr_ch = ch
				curr_loc = 0
				curr_list = chr_exon.get(ch,[])
				
			pos = int(pos)
			
			while curr_loc < len(curr_list) and curr_list[curr_loc][1] < pos:
				curr_loc = curr_loc + 1
					
			if curr_loc < len(curr_list) and curr_list[curr_loc][0] <= pos and curr_list[curr_loc][1] >= pos:
				print l,
