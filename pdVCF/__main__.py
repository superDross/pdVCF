import argparse
import subprocess
import os

# global
here = os.path.realpath(__file__)
test_dir = "/".join(here.split("/")[:-2]) + "/test/"


def cli():
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['test']:
        test = test_dir+"test_filters.py"
        subprocess.call(['python3', test])
        return

    if args['version']:
        print('pdVCF 0.1')
        return

    if args['optimise']:
        activate_line_profiler()
        return


def get_parser():
    parser = argparse.ArgumentParser(description='tool to test pdVCF')
    parser.add_argument('--version', action='store_true', help='display the current version')
    parser.add_argument('--test', action='store_true', help='test installation')
    parser.add_argument('--optimise', action='store_true', help='optimise vcf2dataframe using line_profiler')
    return parser


def activate_line_profiler():
    ''' Execute line_profiler upon vcf2dataframe.py
    '''
    # path to vcf2dataframe.py and command to add to vcf2dataframe for optimisation test
    f = '/'.join(here.split('/')[:-1])+"/vcf2dataframe.py"
    command = "$a" + "vcf2dataframe('" + test_dir + "vcfs/testing4.vcf" + "', UID=True)"

    # add @profile above every function, add command to end of file and execute line_profiler
    subprocess.call(['sed', '-i', '/def /i @profile', f])
    subprocess.call(['sed', '-i', '-e', command, f])
    subprocess.call(['kernprof', '-v', '-l', f])

    # remove @profile and command from file and remove .lprof file
    subprocess.call(['sed', '-i', '/@profile/d', f])
    subprocess.call(['sed', '-i', '$d', f])
    os.remove(f.split("/")[-1] + '.lprof')
    return 


if __name__ == '__main__':
    cli()
