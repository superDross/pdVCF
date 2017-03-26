''' A collection functions used to compare VCF objects
'''

def common_variants(vcf1, vcf2):
    ''' Find common variants between two VCF objects
        and return common variants in a list.
    '''
    return list(set.intersection(set(vcf1.vcf.index.values), 
                                 set(vcf2.vcf.index.values)))    
        
    
