import argparse
import os
from collections import Counter


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2021.06.18'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='remove redundancy and select representative NRS by cluster and kalign')
parser.add_argument('-dir_score', '--dir_kalign_score', type=str, help='path to the folder containing kalign score')
parser.add_argument('-vcf', '--merge_vcf', type=str, help='merge VCF by jasmine')
parser.add_argument('-dir', '--all_vcf_dir', type=str, help='path to the folder containing final.vcf file')
parser.add_argument('-out_name', '--output_name', type=str)


def main():
    args = parser.parse_args()
    dir_kalign_score = args.dir_kalign_score
    merge_vcf = args.merge_vcf
    all_vcf_dir = args.all_vcf_dir
    output_name = args.output_name
    singleton_set = read_merge_vcf(merge_vcf)
    all_nrs_dict = read_all_vcf(all_vcf_dir)
    total_representative_set = read_kalign_score(dir_kalign_score, all_nrs_dict)
    shuchu_vcf(total_representative_set, singleton_set, all_nrs_dict, output_name)


def shuchu_vcf(total_representative_set, singleton_set, all_nrs_dict, output_name):
    temp_shuchu = []
    for i, per in all_nrs_dict.items():
        if i in total_representative_set or i in singleton_set:
            temp_shuchu.append([per[0], per[1], per[2]])
    temp_shuchu = sorted(temp_shuchu, key=lambda x: (x[0], int(x[1])))
    f_shuchu = open('%s.vcf' % output_name, 'w')
    f_shuchu.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % output_name)
    for z in temp_shuchu:
        f_shuchu.write('\t'.join(all_nrs_dict[z[2]]) + '\n')
    f_shuchu.close()
    os.system('cat header.vcf %s.vcf > %s.final.vcf' % (output_name, output_name))


def read_all_vcf(all_vcf_dir):
    all_nrs_dict = {}
    vcf_file_list = os.listdir(all_vcf_dir)
    vcf_file_list = [i for i in vcf_file_list if i.endswith('final.vcf')]
    for file in vcf_file_list:
        file_path = os.path.join(os.path.abspath(all_vcf_dir), file)
        with open(file_path, 'r') as f:
            for raw_line in f:
                if not raw_line.startswith('#'):
                    line = raw_line.strip('\n').split('\t')
                    all_nrs_dict[line[2]] = line
    return all_nrs_dict


def read_merge_vcf(merge_vcf):
    singleton_set = set()
    with open(merge_vcf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip('\n').split('\t')
                info = line[7].split(';')
                info_dict = {}
                for x in info:
                    x = x.split('=')
                    info_dict[x[0]] = x[1]
                if info_dict['SUPP'] == '1':
                    singleton_set.add(line[2])
    return singleton_set


def read_kalign_score(dir_kalign_score, all_nrs_dict):
    total_representative_set = set()
    score_file_list = os.listdir(dir_kalign_score)
    score_file_list = [i for i in score_file_list if i.endswith('score')]
    for file in score_file_list:
        file_path = os.path.join(os.path.abspath(dir_kalign_score), file)
        # kalign score 已使用 parallel 并行化计算
        with open(file_path, 'r') as f:
            first_line = f.readline()
        # 获取 score 返回值，并优先选择 M（华西数据）
        temp_raw = first_line.strip('\n').split('\t')
        temp_length = [len(all_nrs_dict[i.split(',')[0]][4]) for i in temp_raw]
        most_frequent_number = Counter(temp_length).most_common(1)[0][1]
        most_frequent_length = [a for a, b in Counter(temp_length).items() if b == most_frequent_number]
        most_frequent_length = set(most_frequent_length)
        print(temp_raw)
        print(Counter(temp_length), most_frequent_length, most_frequent_number)
        temp_line = [(i.split(',')[0], float(i.split(',')[1])) for i in temp_raw
                     if len(all_nrs_dict[i.split(',')[0]][4]) in most_frequent_length]
        new_line = sort_again(temp_line)
        new_line = sorted(new_line, key=lambda x: (x[2], x[1]), reverse=True)
        print(new_line)
        print('-')
        total_representative_set.add(new_line[0][0])
    return total_representative_set


def sort_again(temp_line):
    new_line = []
    for i in temp_line:
        if 'M' in i[0]:
            new_line.append((i[0], i[1], 'M'))
        else:
            new_line.append((i[0], i[1], 'CN'))
    return new_line


if __name__ == '__main__':
    main()

