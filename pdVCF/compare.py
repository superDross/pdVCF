''' A collection of functions used to compare VCF objects
'''
from plumbum.cmd import tabix
from pdVCF.vcf2dataframe import vcf2dataframe
import os
import re


def common_variants(vcf1, vcf2):
    ''' Find common variants between two VCF objects
        and return common variants in a list.
    Args:
        vcf1: first VCF object
        vcf2: second VCF object
    '''
    return list(set.intersection(set(vcf1.vcf.index.values), 
                                 set(vcf2.vcf.index.values)))    
 


def multi2bi(df):
    ''' Convert multi-allelic UIDs, deriving from a 
        pdVCF, in a list to bi-allelic UIDs.

    Args:
        variants: a pdVCF dataframe

    Returns:
        list of UIDs from pdVCF with multi-allelic
        ones converted to bi-allelic e.g.

        ['2:1234-G/C,T'] -> ['2:1234-G/C', '2:1234-G/T']
    
    '''
    variants = df.index.tolist()
    result = variants[:]

    for variant in variants:
        if ',' in variant:
            multi = re.split('[,/]', variant)
            bi = ["/".join((multi[0], x)) for x in multi[1:]]
            
            result.pop(result.index(variant)) # Removes multi-allelic variant from list
            result = result + bi

    return result


