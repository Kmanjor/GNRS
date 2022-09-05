import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='create command of AGE')
parser.add_argument('-raw_txt', '--raw_NRS_txt', type=str, help='NRS_raw.txt from extract_unmap')
parser.add_argument('-id_map', '--id_map_file', type=str, help='GCF_000001405.39_GRCh38.p13.full.report.txt',
                    default='/home/litong/software/tandem-genotypes-1.7.1/hg38/GCF_000001405.39_GRCh38.p13.full.report.txt')
parser.add_argument('-sample_name', '--sample_name', type=str, help='sample name')
parser.add_argument('-per_fasta_dir', '--per_fasta_dir_age', type=str)


def main():
    args = parser.parse_args()
    raw_nrs_txt = args.raw_NRS_txt
    # id_map_file = args.id_map_file
    sample_name = args.sample_name
    per_fasta_dir_age = args.per_fasta_dir_age
    # chr_dict = read_id_map(id_map_file)
    read_txt(raw_nrs_txt, sample_name, per_fasta_dir_age)


def read_id_map(id_map_file):
    chr_dict = {}
    with open(id_map_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip('\n').split('\t')
                if 'NC_' in line[6]:  # 只保留了1-24号染色体，与M
                    chr_dict[line[6]] = line[-1]
    return chr_dict


def read_txt(raw_nrs_txt, sample_name, per_fasta_dir_age):
    with open(raw_nrs_txt, 'r') as f:
        for line in f:
            if not line.startswith('#NRS_ID'):
                line = line.strip('\n').split('\t')
                if line[1] == 'partial':
                    if line[3] == line[5]:
                        if abs(int(line[6]) - int(line[4])) <= 2000:
                            f_shuchu = open('%s.AGE.sh' % line[0], 'w')
                            if int(line[4]) <= int(line[6]):
                                f_shuchu.write('/home/litong/software/age_v0.4/src/age_align -both '
                                               '/home/litong/database/S288C/per_fasta_for_AGE/R64_%s.fa '
                                               '%s%s_%s.fa -coor1=%s-%s -coor2=%s-%s > %s.age'
                                               % (line[3], per_fasta_dir_age, sample_name, line[7],
                                                  int(line[4]) - 1000, int(line[6]) + 1000, int(line[8]) - 1000,
                                                  int(line[9]) + 1000, line[0]))
                                print('%s.AGE.sh' % line[0])
                            else:
                                f_shuchu.write('/home/litong/software/age_v0.4/src/age_align -both '
                                               '/home/litong/database/S288C/per_fasta_for_AGE/R64_%s.fa '
                                               '%s%s_%s.fa -coor1=%s-%s -coor2=%s-%s > %s.age'
                                               % (line[3], per_fasta_dir_age, sample_name, line[7],
                                                  int(line[6]) - 1000, int(line[4]) + 1000, int(line[8]) - 1000,
                                                  int(line[9]) + 1000, line[0]))
                                print('%s.AGE.sh' % line[0])
                            f_shuchu.close()


if __name__ == '__main__':
    main()

