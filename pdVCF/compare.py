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
 

def check_vcf_for_variants(pos_snp, vcf):
    ''' Check if a list of variant positions and SNPs
        are present within a given vcf. If so, remove 
        them from the list.

    Args:
        pos_snp: list of UIDs deriving from a VCF object
                    e.g. ['1:123456-G/T', '2:726356-C/A']
        vcf: vcf file to check for pos_snp

    Returns:
        the pos_snp original list with the variants present in the given
        vcf removed

    Notes:
        tabix needs to be installed for this to work
    '''
    pos_snp = [(x.split("-")[0], x.split("-")[1]) for x in pos_snp]
    # make a copy of the orginal list to modify each iteration of the loop
    result = pos_snp[:]

    for pos, snp in pos_snp:
        pos_range = "{}-{}".format(pos, pos.split(":")[1])
        tab = tabix['-h', vcf, pos_range]().rstrip("\n")

        # need to write tabix output to file to allow it to work with pdVCF. I would prefer to just parse it into pdVCF instead.
        if tab:
            with open("temp.vcf", 'w') as f:
                f.write(tab)
            
            # check if the tabix output is empty. Parsing empty vcf into pdVCF results in errors.
            check_empty = [x for x in open("temp.vcf").read().split("\n") if not x.startswith("#")]

            if check_empty:
                # check if variant is present in tabix output
                df = vcf2dataframe('temp.vcf', info_level=False, UID=True, genotype_level=False)
                uid = "{}-{}".format(pos, snp)
                variants = multi2bi(df) # get all variants present in df
                #variants = df.index.tolist()

                if uid in variants:
                    result.pop(result.index((pos, snp)))

            else:
                pass
    # clean up
    if os.path.isfile("temp.vcf"):
        os.remove("temp.vcf")

    return result

  

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


