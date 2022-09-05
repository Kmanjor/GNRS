import argparse
from Bio import SeqIO


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2021.06.18'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='select representative sequence from multiple sequence alignment')
parser.add_argument('-afa', '--afa_kalign', type=str, help='path to fasta format result of kalign')


def main():
    args = parser.parse_args()
    afa_kalign = args.afa_kalign

    total_seq_list_dict = {}
    for i in SeqIO.parse(afa_kalign, 'fasta'):
        total_seq_list_dict[i.id] = list(str(i.seq))

    sort_list = []
    for current_name in total_seq_list_dict:
        score = calculate_score(current_name, total_seq_list_dict)
        sort_list.append((current_name, score))

    sort_list = sorted(sort_list, key=lambda x: x[1], reverse=True)
    print('\t'.join(['%s,%s' % (i[0], i[1]) for i in sort_list]))

    # 优先选择 two end anchor的序列，在得分最高的基础上优化结果


def calculate_score(current_name, total_seq_list_dict):
    score = 0
    current_seq = total_seq_list_dict[current_name]
    for seq_name in total_seq_list_dict:
        if seq_name != current_name:
            seq_seq = total_seq_list_dict[seq_name]
            for idx in range(len(current_seq)):
                current_base = current_seq[idx]
                seq_base = seq_seq[idx]
                # 事实上可以使用 NCBI 的 Multiple Sequence Alignment Viewer 验证结果
                if current_base != '-' and seq_base != '-':
                    if current_base == 'N':
                        score -= 0.5
                    else:
                        if current_base == seq_base:
                            score += 2
                        else:
                            score -= 1
    return score


if __name__ == '__main__':
    main()

