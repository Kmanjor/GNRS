import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2021.06.18'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='select representive NRS based on jasmine VCF by kalign')
parser.add_argument('-vcf', '--merge_vcf', type=str, help='merge VCF by jasmine')
parser.add_argument('-dir', '--all_vcf_dir', type=str, help='path to the folder containing final.vcf file')
parser.add_argument('-temp', '--temporary_file_name', type=str)


def main():
    args = parser.parse_args()
    merge_vcf = args.merge_vcf
    all_vcf_dir = args.all_vcf_dir
    temporary_file_name = args.temporary_file_name
    clique_dict = read_merge_vcf(merge_vcf)
    all_nrs_dict = read_all_vcf(all_vcf_dir)
    get_cluster_fasta(clique_dict, all_nrs_dict, temporary_file_name)


def get_cluster_fasta(clique_dict, all_nrs_dict, temporary_file_name):
    for i, per in clique_dict.items():
        shuchu = []
        for x in per:
            shuchu.append(all_nrs_dict[x])
        SeqIO.write(shuchu, '%s_%s.fa' % (temporary_file_name, i), 'fasta')


def read_all_vcf(all_vcf_dir):
    all_nrs_dict = {}
    vcf_file_list = os.listdir(all_vcf_dir)
    vcf_file_list = [i for i in vcf_file_list if i.endswith('final.vcf')]
    for file in vcf_file_list:
        file_path = os.path.join(os.path.abspath(all_vcf_dir), file)
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip('\n').split('\t')
                    temp = SeqRecord(Seq(line[4]), id=line[2], name='', description='')
                    all_nrs_dict[line[2]] = temp
    return all_nrs_dict


def read_merge_vcf(merge_vcf):
    clique_dict = {}
    n = 0
    with open(merge_vcf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip('\n').split('\t')
                info = line[7].split(';')
                info_dict = {}
                for x in info:
                    x = x.split('=')
                    info_dict[x[0]] = x[1]
                if info_dict['SUPP'] != '1':  # 去除singleton，否则kalign会报错
                    id_list = info_dict['IDLIST'].split(',')
                    n += 1
                    clique_dict['clique_%s' % n] = id_list
    return clique_dict


if __name__ == '__main__':
    main()

