import argparse


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2021.06.09'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='extract unmap using quast')
parser.add_argument('-quast_unalign_info', '--quast_unalign_info', type=str, help='quast unalign info')
parser.add_argument('-quast_coords', '--quast_coords_filtered', type=str, help='quast coords.filtered file')
parser.add_argument('-out_txt_NRS', '--out_txt_NRS', type=str, help='path of NRS txt')
parser.add_argument('-sample_name', '--sample_name', type=str, help='sample name')


def main():
    args = parser.parse_args()
    quast_unalign_info = args.quast_unalign_info
    quast_coords_filtered = args.quast_coords_filtered
    out_txt_nrs = args.out_txt_NRS
    sample_name = args.sample_name
    coords_dict = read_coords(quast_coords_filtered)
    info_list, ctg_length_dict = read_unalign(quast_unalign_info)
    write_to_out(coords_dict, info_list, out_txt_nrs, ctg_length_dict, sample_name)


def read_coords(quast_coords_filtered):
    coords_dict = {}
    with open(quast_coords_filtered, 'r') as f:
        for line in f:
            line = line.strip('\n').split(' | ')[:-1]
            line = [i.split(' ') for i in line]
            ctg_name = line[4][1]
            chrom_name = line[4][0]
            ctg_coords1 = int(line[1][0])
            ctg_coords2 = int(line[1][1])
            if ctg_coords1 < ctg_coords2:
                coords_dict[(ctg_name, ctg_coords1 - 1)] = (chrom_name, int(line[0][0]))
                coords_dict[(ctg_name, ctg_coords2 + 1)] = (chrom_name, int(line[0][1]))
            else:
                coords_dict[(ctg_name, ctg_coords1 + 1)] = (chrom_name, int(line[0][0]))
                coords_dict[(ctg_name, ctg_coords2 - 1)] = (chrom_name, int(line[0][1]))
    return coords_dict


def read_unalign(quast_unalign_info):
    ctg_length_dict = {}
    info_list = []
    with open(quast_unalign_info, 'r') as f:
        for line in f:
            if not line.startswith('Contig'):
                line = line.strip('\n').split('\t')
                info = line[4].split(',')
                for per in info:
                    per = per.split('-')
                    length = int(per[1]) - int(per[0]) + 1
                    if length >= 50:  # 调整为200bp，以获得更多的新序列
                        info_list.append((line[0], int(per[0]), int(per[1]), line[3]))
                        # ctg  start  end  type
                        ctg_length_dict[line[0]] = int(line[1])
    return info_list, ctg_length_dict


def write_to_out(coords_dict, info_list, out_txt_nrs, ctg_length_dict, sample_name):
    f_shuchu = open(out_txt_nrs, 'w')
    n = 0
    f_shuchu.write('#NRS_ID\tNRS_TYPE\tNRS_LEN\tCHR1\tCHROM_1\tCHR2\tCHROM_2\tCONTIG\tCONTIG_1\tCONTIG_2\tCONTIG_LEN\n')
    for i in info_list:
        if ctg_length_dict[i[0]] >= 20000:  # 最小的contig，需要大于raw read的长度
            if i[3] == 'full':
                if i[2] - 20000 >= 50:  # 调整为200bp，以获得更多的新序列
                    n += 1
                    nrs_id = '%s_NRS_%s' % (sample_name, n)
                    f_shuchu.write('\t'.join([nrs_id, i[3], str(i[2] - 20000), '.', '.', '.', '.', i[0],
                                              str(i[1] + 10000), str(i[2] - 10000), str(ctg_length_dict[i[0]])]) + '\n')
            else:
                if i[1] >= 10000:
                    if i[2] + 10000 <= ctg_length_dict[i[0]]:
                        n += 1
                        nrs_id = '%s_NRS_%s' % (sample_name, n)
                        f_shuchu.write('\t'.join([nrs_id, i[3], str(i[2] - i[1] + 1), coords_dict[(i[0], i[1])][0],
                                                  str(coords_dict[(i[0], i[1])][1]), coords_dict[(i[0], i[2])][0],
                                                  str(coords_dict[(i[0], i[2])][1]), i[0], str(i[1]), str(i[2]),
                                                  str(ctg_length_dict[i[0]])]) + '\n')
    f_shuchu.close()


if __name__ == '__main__':
    main()

