import vcf2dataframe
import argparse
import subprocess
import os


def get_parser():
    parser = argparse.ArgumentParser(description='tool to test pdVCF')
    parser.add_argument('--version', action='store_true', help='display the current version')
    parser.add_argument('--test', action='store_true', help='test installation')
    parser.add_argument('--optimise', action='store_true', help='optimise pdVCF usng line_profiler')
    parser.add_argument('-a')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    here = os.path.realpath(__file__)

    if args['test']:
        test = "/".join(here.split("/")[:-2]) + "/test/test_filters.py"
        subprocess.call(['python3', test])
        return

    if args['version']:
        print(__version__)
        return

    if args['optimise']:
        directory = '/'.join(here.split('/')[:-1])
        f = directory+"/vcf2dataframe.py"
        subprocess.call(['sed', '-i', '/def /i @profile', f])
        subprocess.call(['kernprof', '-v', '-l', f])
        subprocess.call(['sed', '-i', 's/@profile//g', f])
        subprocess.call(['rm', f+'.lprof'])

    if args['a']:
        subprocess.call(['ls', '/home/david/'])

if __name__ == '__main__':
    cli()
