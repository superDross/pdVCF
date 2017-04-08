from pdVCF.vcf2dataframe import vcf2dataframe
import pdVCF.conditions as cond
import pandas as pd
import re
import operator

pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', -1)


class VCF(object):
    ''' A VCF file stored as a Pandas DataFrame
    
    Atrributes:
        vcf: vcf file to be converted to a Pandas DataFrame or a VCF object
        genotype_level: place the genotype information into a second level column index
        info_level: place the info IDs into a second level column index
        UID: rename index to a unique variant identifier
        
    Notes:
        it is not recommended to alter the boolean attributes when initilising
        a VCF object, as it may break method functionality and limit data 
        manipulation of the resulting object.
    '''
    def __init__(self, vcf, genotype_level=True, info_level=True, UID=True):
        
        self.vcf = vcf2dataframe(vcf, genotype_level=genotype_level,
                                 info_level=info_level, UID=UID)


    def get_samples(self):
        ''' Get all sample names within the vcf and return as a list
        '''
        return self.vcf.xs('DP', level=1, axis=1).columns.tolist()


    def get_genotype(self, gen):
        ''' Access specific genotype information across samples
            in the vcf.
        Args:
            gen: genotype attribute of interest in string format e.g 'DP'
        '''
        return self.vcf.xs(gen, level=1, axis=1)

    
    def get_info(self, info):
        ''' Return INFO field of interest e.g. 'AC'
        '''
        return self.vcf['INFO'][info]


    def subset(self, sams, exclude_ref=False, remove_uncalled=True):
        ''' Subset a multisample VCF by a given sample(s).
        Args:
            vcf: Pandas DataFrame VCF
            sams: list of samples to subset the vcf for
            exlude_ref: remove variant if all GT values for subset are 0/0
            remove_uncalled: remove variant if all GT values for subset are ./.

        Returns:
            subsetted Pandas DataFrame VCF
        '''
        # split variant and genotype information 
        sams = sams if isinstance(sams, list) else [sams]
        genotype = self.vcf[sams]
        num_info = self.vcf['INFO'].columns.shape[0]
        variant = self.vcf.ix[:,:8+num_info]

        GT = genotype.xs('GT', level=1, axis=1)
        uncalled= []

        if remove_uncalled:
            uncalled = GT[GT[sams] == './.'].dropna().index.tolist() 

        if exclude_ref:
            uncalled += GT[GT[sams] == '0/0'].dropna().index.tolist() 

        sub = pd.concat([variant, genotype], axis=1)
        self.vcf = sub.drop(uncalled)
        return self.vcf


    def filter_vcf(self, cond_list, op="&", how='any'):
        ''' Filter out variants in the VCF that do not meet
            the given conditions.

        Args:
            cond_list: list of conditions to be met e.g ['DP >= 50', 'AB < 0.3']
            op: opertor to apply to cond_list e.g. "&" means all conditions in cnd_ist must be met
            how: 'all' mean all samples for a given variant must pass the applied conditions to remain in the VCF. 'any' means at least one sample for a given variant must meet the given conditions to remain in the VCF.

        Notes:
            filtering INFO, mandatory and Gentotype fields must be performed seperately.
            does not filter multi-allelic variants.
        '''

        # check for errors
        self.check_fields(cond_list)
        
        # apply all conditions to the vcf
        cond_list = [cond.create_condition(self.vcf, c) for c in cond_list]
        conditions = cond.combine_conditions(cond_list, op)

        if how == 'all':
            filtered = conditions.T.all() 

        elif how == 'any':
            filtered = conditions.T.any()

        self.vcf = self.vcf[filtered]

        return self.vcf


    def check_fields(self, l):
        ''' Ensure the given list of conditions does not apply to both
            the INFO and genotype fields of the VCF object.
        '''
        fields = [x.split(" ")[0] for x in l]
        info = self.vcf.INFO.columns.tolist() 
        geno = ['DP', 'GT', 'AD', 'GQ', 'PL', 'AB']
        mandatory = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
        
        if len(set(fields).intersection(geno) ) > 0 and len(set(fields).intersection(info) ) > 0:
            raise ValueError("Elements in {} are present in both the genotype and INFO fields. Filter genotype and INFO fields seperately.".format(fields))
        
        if len(set(fields).intersection(geno) ) > 0 and len(set(fields).intersection(mandatory) ) > 0:
            raise ValueError("Elements in {} are present in both the genotype and mandatory fields. Filter genotype and mandatory fields seperately.".format(fields))

        if len(set(fields).intersection(info) ) > 0 and len(set(fields).intersection(mandatory) ) > 0:
            raise ValueError("Elements in {} are present in both the INFO and mandatory fields. Filter mandatory and INFO fields seperately.".format(fields))




    def filter_flag(self, flag='PASS', include=True):
        ''' Include or exclude variants with the given flag
            in the VCF FILTER field.
        Args:
            flag: defaults to 'PASS'
            include: include in the vcf if true, otherwise exclude
        '''
        if include:
            self.vcf = self.vcf[self.vcf['FILTER'] == flag]
        else:
            self.vcf = self.vcf[self.vcf['FILTER'] != flag]
    
        return self.vcf
        
    
    def indels(self, include=True):
        ''' Remove or filter for indels from vcf.
        Args:
            include: include in the vcf if true, otherwise exclude
        '''
        alt_mask = (self.vcf.ALT.str.len() == 1) | (self.vcf.ALT.str.contains(','))
        ref_mask = (self.vcf.REF.str.len() == 1) | (self.vcf.REF.str.contains(','))

        if include:
            self.vcf = self.vcf[~alt_mask & ref_mask]
        else:
            self.vcf = self.vcf[alt_mask & ref_mask]

        return self.vcf
    

    def biallelic(self):
        ''' Filter for biallelic variants only.
        '''
        self.vcf = self.vcf[self.vcf.ALT.str.split(',').str.len() == 1]
        return self.vcf
    
    
    def multiallelic(self):
        ''' Filter for multiallelic variants only.
        '''
        self.vcf = self.vcf[self.vcf.ALT.str.split(',').str.len() > 1]
        return self.vcf
    
    
    def positions(self, positions, include=True):
        ''' Include or exclude variants that lie within 
            the given position(s) or position ranges.
        
        Args:
            positions: a position, position ranges or list of the two e.g.
                            ['1:2234385', '1:2235901-2240000']
            include: include in the vcf if true, otherwise exclude
        '''
        positions = positions if isinstance(positions, list) else [positions]
        
        # get the variants named in positions
        selected_variants = []
        
        # get the indexes of the varants within pos from the vcf 
        for pos in positions:
            pos = self.pos2range(pos)
            chrom, start, end = [int(x) for x in re.split(r'[:-]', pos)]
            mask = (self.vcf['CHROM'] == chrom) & (self.vcf['POS'] >= start) & (self.vcf['POS'] <= end)
            variants = self.vcf[mask].index.tolist()
            selected_variants.append(variants)
            
        # flatten
        selected_variants = [y for x in selected_variants for y in x]
        
        if not include:
            selected_variants = list(set(self.vcf.index.tolist()) - set(selected_variants))
            
        self.vcf = self.vcf.loc[self.natural_sort(selected_variants)] 

        return self.vcf
    
    
    def chromosome(self, chrom, include=True):
        ''' Include or exculde variants within a given chromosome.

        Args: 
            chrom: chromosome
            include in the vcf if True, otherwise exclude
        '''
        if include:
            self.vcf = self.vcf[self.vcf['CHROM'] == chrom]
        else:
            self.vcf = self.vcf[self.vcf['CHROM'] != chrom]

        return self.vcf


    @staticmethod
    def natural_sort(l): 
        ''' Sort a list in human natural alphanumerical
            order.
        '''
        convert = lambda text: int(text) if text.isdigit() else text.lower() 
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        return sorted(l, key = alphanum_key)
    
    
    @staticmethod
    def pos2range(pos, num=0): 
        ''' Alter a genomic position to a genomic range
            e.g. 2:1234 becomes 2:1234-1234
        '''
        if "-" not in pos:
            split = pos.split(":")
            pos = "{}:{}-{}".format(split[0], split[1], int(split[1])+num)
            return pos
        else:
            return pos
        
        


# code graveyard below


#    def filter_genotype(self, minDP=None, minGQ=None, minAB=None):
#        ''' Filter for variants in which all the samples in the given vcf 
#            meet the minimum genotype values given.
#        
#        Args:
#            minDP: minimum variant depth
#            minGQ: minimum genotype quality
#            minAB: minimum allele balance
#        
#        Notes:
#            Doesn't handle multiallelic information properly and
#            will filter for this first ALT value e.g. if DP = 12,1,100
#            it will be filtered out even if minDP=30.
#        '''
#        # split variant and genotype information
#        num_info = self.vcf['INFO'].columns.shape[0]
#        variant = self.vcf.ix[:,:8+num_info]
#        genotype = self.vcf.ix[:,9+num_info:]
#
#        # store all variants that don't meet the minimum value given for the args here
#        below_min = []
#        
#        if minDP:
#            DP = genotype.xs('DP', level=1, axis=1).fillna(0)
#            above_min = DP[DP >= minDP] 
#            below_min += DP[above_min.isnull().any(axis=1)].index.tolist()
#            
#        if minGQ:
#            GQ = genotype.xs('GQ', level=1, axis=1).fillna(0)
#            above_min = GQ[GQ >= minGQ]
#            below_min +=  GQ[above_min.isnull().any(axis=1)].index.tolist()
#        
#        if minAB:
#            AB = genotype.xs('AB', level=1, axis=1).fillna(0)
#            above_min = AB[AB >= minAB]
#            below_min +=  AB[above_min.isnull().any(axis=1)].index.tolist()
#        
#        # remove variants that don't meet the requirements from the vcf
#        self.vcf = self.vcf.drop(below_min)
#        return self.vcf
#    
#
#
#
#
#    def filter_info(self, field, value):
#        ''' Filter for variants that are above the given value
#            (if value is number) or are equal to the given value
#            (if value is string).
#            
#        Args:
#            field: INFO field of interest
#            value: string or int value to test the field with
#            
#        Notes:
#            Doesn't handle multiallelic information properly and
#            will filter any thing that has this e.g. if AC = 12,34
#            it will be filtered out even if value=1 as the AC with
#            a comma cant be converted to an int so becomes np.nan
#        '''
#        if isinstance(value, int) or isinstance(value, float):
#            mask = pd.to_numeric(self.vcf['INFO'][field], errors='coerrce') >= value
#            self.vcf = self.vcf[mask]
#            return self.vcf
#        
#        elif isinstance(value, str):
#            mask= self.vcf['INFO'][field] == value
#            self.vcf = self.vcf[mask]
#            return self.vcf
#
#    def passing_variants(self, field, threshold, le=False, info=True): 
#        ''' filter out variants of an info/genotype field in which 
#            any sample in said field is equal to or above the given 
#            threshold.
#        
#        Args:
#            field: genotype of INFO field i.e. 'DP'
#            value: integer
#            le: retrieve variants in which the fields are
#                equal to or less than the given threshold 
#            info: search for the given field within INFO column,
#                  otherwise search within the samples columns
#                  
#        '''
#        if info:
#            fields = pd.to_numeric(self.vcf['INFO'][field].T, errors='coerrce')
#            passed = fields[fields >= threshold]
#            variants = passed.index.tolist()
#        else:
#            fields = self.vcf.xs(field, level=1, axis=1).T
#            if 'INFO' in fields.index.values.tolist():
#                fields = fields.drop('INFO')
#            passed = fields[fields >= threshold].dropna(how='all', axis=1)
#            variants = passed.columns.tolist()
#        
#        if le:
#            variants = list(set(self.vcf.index.tolist()) - set(variants))
#        
#        self.vcf = self.vcf.ix[variants]
#        
#        return self.vcf
#

