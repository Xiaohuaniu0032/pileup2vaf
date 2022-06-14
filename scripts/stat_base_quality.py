from __future__ import division
import os
import sys
import pysam
from collections import defaultdict


bam = sys.argv[1]
print(bam)

samfile = pysam.AlignmentFile(bam, "rb")

total_base = 0
seq_num = 0
unmap_read_num = 0
mapq_20_read_num = 0

quality_dict = defaultdict(lambda:0)

for read in samfile.fetch():
	# get all aligned reads
	if read.is_secondary or read.is_supplementary:
		# skip second or sup reads, but include dup reads
		next

	seq_num += 1

	mapq = read.mapping_quality
	if mapq >= 20:
		mapq_20_read_num += 1

	quality = list(read.query_qualities)
	
	for q in quality:
		quality_dict[q] += 1

	len_seq = len(quality)
	total_base += len_seq

print("total reads num:", end=" ")
print(seq_num)

print("total base:", end=" ")
print(total_base)


q20_base_num = 0
for q in sorted(quality_dict):
	if q >= 20:
		q20_base_num += quality_dict[q]

print("Q20 base num:", end=" ")
print(q20_base_num)

q20_pct = q20_base_num/total_base*100
print("Q20 pct:", end=" ")
print(q20_pct)

print("mapq >= 20 reads num:", end= " ")
print(mapq_20_read_num)
mapq_20_pct = mapq_20_read_num/seq_num*100
print("mapq >= 20 reads pct:", end= " ")
print(mapq_20_pct)

print("base.quality","base.num")
for q in sorted(quality_dict):
	print(q,quality_dict[q])
