''' Functions utilised for error handeling 
'''

def check_fields(vcf, conds):
    ''' Ensure the given list of conditions does not apply to both
        the INFO and genotype fields of the VCF object.

    Args:
        vcf: pandas VCF
        conds: list of filtering conditions
    '''
    fields = [x.split(" ")[0] for x in conds]
    info = vcf.INFO.columns.tolist() 
    geno = ['DP', 'GT', 'AD', 'GQ', 'PL', 'AB']
    mandatory = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    
    if len(set(fields).intersection(geno)) > 0 and len(set(fields).intersection(info)) > 0:
        raise ValueError("Elements in {} are present in both the genotype and INFO fields. Filter genotype and INFO fields seperately.".format(fields))
    
    if len(set(fields).intersection(geno)) > 0 and len(set(fields).intersection(mandatory)) > 0:
        raise ValueError("Elements in {} are present in both the genotype and mandatory fields. Filter genotype and mandatory fields seperately.".format(fields))

    if len(set(fields).intersection(info)) > 0 and len(set(fields).intersection(mandatory)) > 0:
        raise ValueError("Elements in {} are present in both the INFO and mandatory fields. Filter mandatory and INFO fields seperately.".format(fields))


