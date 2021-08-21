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

def parse_query(chrom,pos,ref,query):
    '''
    1.count depth
    2.mismatch
    3.ins
    4.del
    5.cal freq
    '''
    BASE = ['A','T','C','G']
    #print(chrom,pos,ref)

    dict_all_temp = defaultdict(lambda:0)
    for item in query:
        dict_all_temp[item.upper()] += 1

    dict_all = {}
    for k,v in dict_all_temp.items():
        dict_all[k] = v

    #print(dict_all)

    dict_alt = defaultdict(lambda:0) # store alt info
    #print(dict_all.keys())
    # check if homo-ref or alt
    if len(dict_all.keys()) == 1:
        for k,v in dict_all.items():
            if k == ref:
                # ref-homo
                alt = 'ref-homo'
                var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                dict_alt[var] = v
            elif k != ref:
                # alt-homo
                alt = k
                var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                dict_alt[var] = v
            else:
                pass
    else:
        # has >= two allels
        for k in dict_all:
        #print(k,ref)
            if k.upper() == ref:
                # skip ref
                pass
            else:
                # alt base
                if len(k) == 1:
                    if k != '*':
                        # SNV
                        alt = k.upper()
                        var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                        dict_alt[var] += dict_all[k]
                elif k[1] == '+':
                    # Ins
                    # 'C+4TACT' => VCF record is: chr pos C CTACT
                    # get alt allele
                    alt_allele = []
                    for i in k:
                        if i in BASE:
                            alt_allele.append(i.upper())

                    alt = "".join(alt_allele)
                    var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                    dict_alt[var] += dict_all[k]
                elif k[1] == '-':
                    # Del
                    # 'G-1C' => VCF record is: chr pos GC G
                    alt_allele = []
                    for i in k:
                        if i in BASE:
                            alt_allele.append(i.upper())

                    del_ref = "".join(alt_allele)
                    del_alt = del_ref[0]
                    var = "%s\t%s\t%s\t%s" % (chrom,pos,del_ref,del_alt)
                    dict_alt[var] += dict_all[k]
                else:
                    pass

    #print(dict_alt)

    return(dict_alt)

def cal_depth(query):
    depth = 0
    for item in query:
        depth += 1

    return(depth)

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
    h = ['Chr','Pos','Ref','Alt','AltNum','Depth','AlleleFrequency']
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
            querybase = pileupcol.get_query_sequences(mark_matches=False,mark_ends=False,add_indels=True)

            # cal depth

            col = pileupcol.reference_pos + 1 # 1-based
            t_chr = t.split(":")[0]
            #print(t_chr,col,querybase)
            reg = t_chr + ':' + str(col) + '-' + str(col) # chr1:2-2
            ref_base = pysam_fa.fetch(region=reg).upper() # ref base
            alt_dict = parse_query(t_chr,col,ref_base,querybase)
            
            depth = cal_depth(querybase)

            # output all alt allele
            for k in alt_dict:
                alt_n = alt_dict[k]
                freq = round(alt_n/depth,3)
                v = "%s\t%s\t%s\t%s" % (k,alt_n,depth,freq)
                of.write(v+'\n')
    of.close()

if __name__ == "__main__":
    main()
