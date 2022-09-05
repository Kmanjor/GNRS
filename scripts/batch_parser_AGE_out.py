import argparse
import os


__author__ = 'Tong Li'
__email__ = 'tli.aioniya@gmail.com'
__data__ = '2021.06.11'


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='batch parse result from AGE')
parser.add_argument('-age_batch_sh', '--AGE_batch_sh', type=str, help='AGE_batch.sh from batch_command_AGE.py')
parser.add_argument('-sample_name', '--sample_name', type=str)


def main():
    args = parser.parse_args()
    age_batch_sh = args.AGE_batch_sh
    sample_name = args.sample_name

    print('#NRS_ID\tNRS_LEN\tCHR1\tCHROM_1\tCHROM_2\tSTRAND\tCONTIG\tCONTIG_1\tCONTIG_2')
    with open(age_batch_sh, 'r') as f:
        for line in f:
            line = line.strip('\n').split('.')[0]
            os.system('python /home/litong/project/NRS_new/src/parser_AGE_out.py -age_result %s/%s.age -nrs_id %s' % (sample_name, line, line))


if __name__ == '__main__':
    main()

