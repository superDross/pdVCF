#!/usr/bin/python3
from pdVCF.pdVCF import Vcf
import pandas as pd
import time


start = time.time()
tsv = pd.read_csv('vcfs/meow.vcf', sep="\t")
end = time.time() - start
print("{} seconds - read_cv()".format(round(end, 2)))


start = time.time()
vcf = Vcf('vcfs/testing3.vcf')
end = time.time() - start
print("{} seconds - vcf2dataframe()".format(round(end, 2)))
