import argparse
import os
from Bio import SeqIO
import gzip


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='convert anchor.txt to VCF for merge')
parser.add_argument('-txt', '--anchor_txt', type=str, help='path to the anchor.txt')
parser.add_argument('-assembly', '--assembly_fasta', type=str, help='assembly fasta by NGS or self polished')
parser.add_argument('-ref', '--reference', type=str, help='reference fasta of hg38', 
                    default='/home/litong/database/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
parser.add_argument('-tr_region', '--TR_region', type=str, help='TR region for annotation',
                    default='/home/litong/software/tandem-genotypes-1.8.2/hg38/simpleRepeat_20210615.sorted.bed')
parser.add_argument('-sample_name', '--sample_name', type=str, help='sample name for prefix')


def main():
    args = parser.parse_args()
    anchor_txt = args.anchor_txt
    assembly_fasta = args.assembly_fasta
    reference = args.reference
    # tr_region = args.TR_region
    sample_name = args.sample_name
    # annotation_tr_dict = tr_annotation(tr_region, anchor_txt, sample_name)
    get_vcf(anchor_txt, assembly_fasta, reference, sample_name)


def get_vcf(anchor_txt, assembly_fasta, reference, sample_name):
    if assembly_fasta.endswith('gz'):
        assembly_fasta = gzip.open(assembly_fasta, 'rt')
    seq_dict = SeqIO.to_dict(SeqIO.parse(assembly_fasta, 'fasta'))
    reference_dict = SeqIO.to_dict(SeqIO.parse(reference, 'fasta'))
    f_shuchu = open('%s.vcf' % sample_name, 'w')
    f_shuchu.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample_name)
    wait_shuchu = []
    with open(anchor_txt, 'r') as f:
        for line in f:
            if not line.startswith('#NRS_ID'):
                line = line.strip('\n').split('\t')
                ref_base = str(reference_dict[line[2].split(' ')[0]].seq)[int(line[3]) - 1]
                if line[5] == '+':
                    alt_seq = str(seq_dict[line[6].split(' ')[0]][int(line[7])-1:int(line[8])].seq)
                else:
                    alt_seq = seq_dict[line[6].split(' ')[0]][int(line[7])-1:int(line[8])].reverse_complement()
                    alt_seq = str(alt_seq.seq)
                # if line[0] in annotation_tr_dict:
                #    newinfo = ['CHR2=%s' % line[2], 'END=%s' % line[4], 'SVTYPE=INS', 'SVLEN=%s' % line[1],
                #              'STRANDS=+-', 'REF_TRF=%s' % annotation_tr_dict[line[0]][2],
                #               'TRRBEGIN=%s' % annotation_tr_dict[line[0]][0],
                #              'TRREND=%s' % annotation_tr_dict[line[0]][1]]
                #else:
                newinfo = ['CHR2=%s' % line[2].split(' ')[0], 'END=%s' % line[4], 'SVTYPE=INS', 'SVLEN=%s' % line[1],
                               'STRANDS=+-']
                newinfo = ';'.join(newinfo)
                newline = [line[2].split(' ')[0], line[3], line[0], ref_base, alt_seq, '.', '.', newinfo, 'GT', '1/1']
                wait_shuchu.append(newline)
    wait_shuchu = sorted(wait_shuchu, key=lambda x: (x[0], int(x[1])))
    for per in wait_shuchu:
        f_shuchu.write('\t'.join(per) + '\n')
    f_shuchu.close()
    os.system('cat header.vcf %s.vcf > %s.final.vcf' % (sample_name, sample_name))


def tr_annotation(tr_region, anchor_txt, sample_name):
    annotation_tr_dict = {}
    f_shuchu = open('%s_anchor.bed' % sample_name, 'w')
    with open(anchor_txt, 'r') as f:
        for line in f:
            if not line.startswith('#NRS_ID'):
                line = line.strip('\n').split('\t')
                f_shuchu.write('\t'.join([line[2], line[3], line[4], line[0]]) + '\n')
    f_shuchu.close()
    os.system('sort -k 1,1 -k2,2n %s_anchor.bed > %s_anchor.sorted.bed' % (sample_name, sample_name))
    os.system('/home/litong/software/bedtools intersect -wa -wb -a %s_anchor.sorted.bed -b %s > %s_anchor.sorted.tr.bed'
              % (sample_name, tr_region, sample_name))
    with open('%s_anchor.sorted.tr.bed' % sample_name, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if line[3] not in annotation_tr_dict:
                annotation_tr_dict[line[3]] = [line[5], line[6], line[7]]
    return annotation_tr_dict


if __name__ == '__main__':
    main()

