import pysam
import sys
import argparse
import os

def parse_args():
        AP =  argparse.ArgumentParser("Detect snv and small indel using pysam pileup method")
        AP.add_argument('-bam',help="bam file",dest="bam")
        AP.add_argument('-bed',help="bed file",dest="bed")
        AP.add_argument('-d',help='max depth',dest='depth',default=8000)
        AP.add_argument('-n',help='samplel name',dest='name')
        AP.add_argument('-mq',help='map quality',dest='mapq',default=30)
        AP.add_argument('-bq',help='base quality',dest='baseq',default=10)
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
    h = ['Chr','Pos','Ref','RefNum','AltNum','Depth','AlleleFrequency']
    hh = '\t'.join(h)
    of.write(hh)
    of.write('\n')

    for t in target:
        for pileupcol in pysam_sam.pileup(region=t,reference=options.fasta,truncate=True,stepper="samtools",fastafile=pysam_fa,max_depth=options.depth,min_base_quality=options.baseq,min_mapping_quality=options.mapq): # return PileupColumn Obj
            '''
            Note: 'all' reads which overlap the region are returned. The first base returned will be the first base of the first read ‘not’ necessarily the first base of the region used in the query. you need to specifiy the 'truncate' param
            '''
            # max depth = 8000X, mapQ = 30, baseQ = 10
            querybase = pileupcol.get_query_sequences(mark_matches=False,mark_ends=False,add_indels=False)
            #print(querybase)
            col = pileupcol.reference_pos + 1 # 1-based
            t_chr = t.split(":")[0]
            reg = t_chr + ':' + str(col) + '-' + str(col) # chr1:2-2
            ref_base = pysam_fa.fetch(region=reg).upper() # ref base

            ref_n, alt_n = 0,0

            for base in querybase:
                if base.upper() == ref_base:
                    ref_n += 1
                else:
                    alt_n += 1

            all_n = ref_n + alt_n

            if all_n > 0:
                alt_freq = round(alt_n/all_n,3)
            else:
                # skip zero depth pos
                continue
               
            val = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (t_chr,col,ref_base,ref_n,alt_n,all_n,alt_freq)
            of.write(val)
    
    of.close()

    # output SNP info
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    snp_file = "%s/%s.BAF.HetPos.txt" % (options.outdir,options.name)
    cmd = "perl %s/bin/get_het_snp.pl %s %s %s" % (bin_dir,options.snp,outfile,snp_file)
    print("CMD is: %s" % (cmd))
    os.system(cmd)


if __name__ == "__main__":
    main()
