import matplotlib.pyplot as plt
import seaborn as sns
import copy

class Plot(object):
    ''' Create plots from a Vcf object.

    Attributes:
        pdvcf: Vcf object

    Note:
        deepcopy is used so filtering doesn't 
        alter the original self.vcf within VCF 
        object
    '''
    def __init__(self, pdvcf):
        self.pdvcf = copy.deepcopy(pdvcf)

    def variants_chrom(self):
        ''' Countplot of number of variants identified
            across all chromosomes.
        '''
        self.pdvcf.remove_scaffolds()
        plt.style.use('seaborn-deep')

        fig, ax = plt.subplots(figsize=(14, 7))
        sns.countplot(data=self.pdvcf.vcf, x='CHROM', palette='GnBu_d')
        
        ax.tick_params(labelsize=15)
        ax.set_ylabel('Variants', fontsize=20)
        ax.set_xlabel('Chromosome', fontsize=20)
        ax.set_title('Variants Identified Across Chromosomes', fontsize=25)
        
        return ax
