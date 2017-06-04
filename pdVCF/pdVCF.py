from .vcf2dataframe import vcf2dataframe
from .Plots import Plot
import pdVCF.conditions as cond
import pdVCF.errors as errors
import pandas as pd
import re

pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', -1)


class Filter(object):
    '''
    '''
    def __init__(self, vcf, genotype_level=True, info_level=True, UID=True):
        
        self.vcf = vcf2dataframe(vcf, genotype_level=genotype_level,
                                 info_level=info_level, UID=UID)

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
        errors.check_fields(self.vcf, cond_list)
        
        # apply all conditions to the vcf
        cond_list = [cond.create_condition(self.vcf, c) for c in cond_list]
        conditions = cond.combine_conditions(cond_list, op)

        if how == 'all':
            filtered = conditions.T.all() 

        elif how == 'any':
            filtered = conditions.T.any()

        self.vcf = self.vcf[filtered]

        return self.vcf


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

    
    def remove_scaffolds(self):
        ''' Remove scaffolding chromosomes from vcf.
        '''
        self.vcf = self.vcf[~self.vcf.CHROM.str.contains('^GL|^KI|^hs')]
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
        


class VCF(Filter):
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
        Filter.__init__(self, vcf, genotype_level, info_level, UID)        

    @property
    def plot(self):
        ''' Gives VCF object access to Plot methods.
        '''
        return Plot(self)

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


 
