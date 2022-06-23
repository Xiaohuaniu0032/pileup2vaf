import pysam
import sys
import argparse
import os
from collections import defaultdict

def parse_args():
        AP =  argparse.ArgumentParser("Detect snv and small indel using pysam pileup method")
        AP.add_argument('-bam',help="bam file",dest="bam")
        AP.add_argument('-bed',help="bed file",dest="bed")
        AP.add_argument('-d',help='max depth',dest='depth',default=100000)
        AP.add_argument('-n',help='samplel name',dest='name')
        AP.add_argument('-mq',help='map quality',dest='mapq',default=0)
        AP.add_argument('-bq',help='base quality',dest='baseq',default=0)
        AP.add_argument('-fa',help="fasta file",dest="fasta",default="/data1/database/b37/human_g1k_v37.fasta")
        AP.add_argument('-od',help="output dir",dest="outdir")
        return AP.parse_args()

BASE = ['A','T','C','G']

def parse_indel(s):
    '''
    s: string
    C+4TACT
    A-1G
    
    return value(ref,alt)
    '''
    val = []
    if '+' in s:
        # Ins
        # 'C+4TACT' => VCF record is: chr pos C CTACT
        alt_allele = []
        for i in s:
            if i in BASE:
                alt_allele.append(i.upper())

        ref = alt_allele[0] 
        alt = "".join(alt_allele)
    else:
        # Del
        # 'G-1C' => VCF record is: chr pos GC G
        alt_allele = []
        for i in s:
            if i in BASE:
                alt_allele.append(i.upper())
        
        ref = "".join(alt_allele) # GC
        alt = ref[0]              # G

    val.append(ref)
    val.append(alt)
    return(val)

def parse_query(chrom,pos,ref,query):
    '''
    1.count depth
    2.mismatch
    3.ins
    4.del
    5.cal freq
  
    query format
    HIV-1 1641 A 16 ['A', 'A', 'A', 'a', 'a', 'g', 'a', 'A', 'A', 'A', 'a', 'g', 'A', 'A', 'a', 'A']
    {'A': 9, 'a': 5, 'g': 2}
    HIV-1 1642 A 16 ['G', 'G', 'G', 'g', 'g', 'a', 'g', 'A-1G', 'G', 'A', 'a', 'a', 'G', 'G', 'g', 'G']
    {'G': 7, 'g': 4, 'a': 3, 'A-1G': 1, 'A': 1}
    HIV-1 1643 G 15 ['A', 'A', 'A', 'a', 'a', 'g+1a', 'a', '*', 'A', 'a', 'g', 'A', 'A', 'a', 'A']
    {'A': 7, 'a': 5, 'g+1a': 1, '*': 1, 'g': 1}
    HIV-1 1644 A 15 ['A', 'A', 'A', 'a', 'a', 'a', 'a', 'A', 'A', 'a', 'a+2at', 'A', 'A', 'a', 'A']
    {'A': 8, 'a': 6, 'a+2at': 1}
    HIV-1 1645 T 15 ['T', 'T', 'T', 't', 't', 't', 't', 'T', 'T', 't', 't', 'T', 'T', 't', 'T']
    {'T': 8, 't': 7}

    '''
    BASE = ['A','T','C','G']

    dict_all = defaultdict(lambda:0)
    for item in query:
        '''
        {'A': 9, 'a': 5, 'g': 2}     => {'A':14, 'G':2}
        {'A': 8, 'a': 6, 'a+2at': 1} => {'A':14, 'A+2AT':1}
        '''
        dict_all[item.upper()] += 1
    #print(dict_all)

    dict_alt = defaultdict(lambda:0) # store alt info
    # check if homo-ref or alt
    if len(dict_all.keys()) == 1:
        # this pos has only one type base
        for k,v in dict_all.items():
            if k == ref:
                # ref-homo
                alt = 'ref-homo'
                var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                dict_alt[var] = 0
            elif k != ref and k in BASE:
                # alt-homo
                alt = k
                var = "%s\t%s\t%s\t%s" % (chrom,pos,ref,alt)
                dict_alt[var] = v
            elif '+' in k or '-' in k:
                # Ins or Del
                indel_info = parse_indel(k)
                ref_tmp = indel_info[0]
                alt_tmp = indel_info[1]
                var = "%s\t%s\t%s\t%s" % (chrom,pos,ref_tmp,alt_tmp)
                dict_alt[var] = v
            else:
                pass
    else:
        # has >= two alleles
        for k in dict_all:
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
                elif '+' in k or '-' in k:
                    # Ins/Del
                    # 'C+4TACT' => VCF record is: chr pos C CTACT
                    # 'G-1C' => VCF record is: chr pos GC G
                    indel_info = parse_indel(k)
                    ref_tmp = indel_info[0]
                    alt_tmp = indel_info[1]
                    var = "%s\t%s\t%s\t%s" % (chrom,pos,ref_tmp,alt_tmp)
                    dict_alt[var] += dict_all[k]
                else:
                    pass

    return(dict_alt)

def depth(query):
    depth = 0
    for i in query:
        if i != '*':
            depth += 1
    return(depth)

def del_depth(query):
    del_depth = 0
    for i in query:
        if i == '*':
            del_depth += 1
    return(del_depth)

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
    h = ['Chr','Pos','Ref','Alt','AltNum','Depth','AltAlleleFrequency']
    hh = '\t'.join(h)
    of.write(hh)
    of.write('\n')

    for t in target:
        print("process %s target" % (t))
        for pileupcol in pysam_sam.pileup(region=t,reference=options.fasta,truncate=True,stepper="samtools",fastafile=pysam_fa,max_depth=options.depth,min_base_quality=options.baseq,min_mapping_quality=options.mapq): # return PileupColumn Obj
            '''
            Note: 'all' reads which overlap the region are returned. The first base returned will be the first base of the first read ‘not’ necessarily the first base of the region used in the query. you need to specifiy the 'truncate' param
            '''
            querybase = pileupcol.get_query_sequences(mark_matches=False,mark_ends=False,add_indels=True)

            # detail depth info
            my_dict = {}
            for i in querybase:
                if i not in my_dict:
                    my_dict[i] = 1
                else:
                    my_dict[i] += 1

            # cal depth
            col = pileupcol.reference_pos + 1 # 1-based
            t_chr = t.split(":")[0]
            reg = t_chr + ':' + str(col) + '-' + str(col) # chr1:2-2
            ref_base = pysam_fa.fetch(region=reg).upper() # ref base
            real_depth = depth(querybase)
            
            total_base = len(querybase)
            del_base = del_depth(querybase)
            if total_base == del_base:
                pass

            #print(t_chr,col,ref_base,my_dict)
            #print(t_chr,col,ref_base,real_depth,querybase)
            
            alt_dict = parse_query(t_chr,col,ref_base,querybase)
            #print(alt_dict)

            for k in alt_dict:
                #print(k)
                alt_n = alt_dict[k] # alt reads number
                freq = round(alt_n/real_depth,3)
                v = "%s\t%s\t%s\t%s" % (k,alt_n,real_depth,freq)
                of.write(v+'\n')
    of.close()

if __name__ == "__main__":
    main()
