import pysam
import sys
import argparse
import os
from collections import defaultdict

def parse_args():
        AP =  argparse.ArgumentParser("Detect snv and small indel using pysam pileup method")
        AP.add_argument('-bam',help="bam file",dest="bam")
        AP.add_argument('-bed',help="bed file",dest="bed")
        AP.add_argument('-d',help='max depth',dest='depth',default=8000)
        AP.add_argument('-n',help='samplel name',dest='name')
        AP.add_argument('-mq',help='map quality',dest='mapq',default=0)
        AP.add_argument('-bq',help='base quality',dest='baseq',default=0)
        AP.add_argument('-fa',help="fasta file",dest="fasta",default="/data1/database/b37/human_g1k_v37.fasta")
        AP.add_argument('-od',help="output dir",dest="outdir")
        return AP.parse_args()

def main():
    options = parse_args()
    pysam_sam = pysam.AlignmentFile(options.bam,"rb")
    pysam_fa = pysam.FastaFile(options.fasta)

    # get target
    target = []
    with open(options.bed,"r") as b_fh:
        for line in b_fh:
            arr = line.rstrip().split('\t')
            t = arr[0] + ':' + str(int(arr[1]) + 1) + '-' + arr[2] # chr:start-end, 1-based
            target.append(t)

    #print(target)

    # output file
    outfile = "%s/%s.variants.xls" % (options.outdir,options.name)
    of = open(outfile,'w')
    h = ['Chr','Pos','Ref','Alt','RefNum','AltNum','Depth','AlleleFrequency']
    hh = '\t'.join(h)
    of.write(hh)
    of.write('\n')

    for t in target:
        print("process %s target" % (t))
        for pileupcol in pysam_sam.pileup(region=t,reference=options.fasta,truncate=True,stepper="samtools",fastafile=pysam_fa,max_depth=options.depth,min_base_quality=options.baseq,min_mapping_quality=options.mapq): # return PileupColumn Obj
            '''
            Note: 'all' reads which overlap the region are returned. The first base returned will be the first base of the first read ‘not’ necessarily the first base of the region used in the query. you need to specifiy the 'truncate' param
            '''
            # max depth = 8000X, mapQ = 30, baseQ = 10
            querybase = pileupcol.get_query_sequences(mark_matches=False,mark_ends=False,add_indels=False)
            #print("querybase is: %s" % (querybase))
            col = pileupcol.reference_pos + 1 # 1-based
            t_chr = t.split(":")[0]
            reg = t_chr + ':' + str(col) + '-' + str(col) # chr1:2-2
            ref_base = pysam_fa.fetch(region=reg).upper() # ref base
            #print(ref_base)
            ref_n, alt_n = 0,0

            all_n = 0
            dict_all = defaultdict(lambda:0)
            for base in querybase:
                dict_all[base.upper()] += 1
                all_n += 1
            
            if all_n == 0:
                # skip zero depth pos
                continue
            
            #print("base info is: %s" % (dict_all))

            # check if exists alt base
            alt_base = []
            ref_n = 0
            for k in dict_all:
                #print("base is: %s" % (k))
                if k != ref_base: # has translated into upper
                    alt_base.append(k)
                if k == ref_base:
                    ref_n += dict_all[k]

            if len(alt_base) != 0:
                for alt in alt_base:
                    alt_n = dict_all[alt]
                    alt_freq = round(alt_n/all_n,3)
                    alt_freq = round(alt_n/all_n,3)
                    if alt == '':
                        # a del
                        val = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t_chr,col,ref_base,'-',ref_n,alt_n,all_n,alt_freq)
                    else:
                        # a mismatch
                        val = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t_chr,col,ref_base,alt,ref_n,alt_n,all_n,alt_freq)
                    of.write(val)

            else:
                # homo ref
                alt_n = 0
                alt_freq = 0
                alt = 'homo-ref'
                val = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t_chr,col,ref_base,alt,ref_n,alt_n,all_n,alt_freq)
                of.write(val)
    
    of.close()

if __name__ == "__main__":
    main()
