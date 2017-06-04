import matplotlib.pyplot as plt
import seaborn as sns

class Plot(object):
    
    def __init__(self, vcf):
        self.vcf = vcf

    def variants_chrom(self):
        self.vcf.remove_scaffolds()
        plt.style.use('seaborn-deep')

        fig, ax = plt.subplots(figsize=(14, 7))
        sns.countplot(data=self.vcf, x='CHROM', palette='GnBu_d')
        
        ax.tick_params(labelsize=15)
        ax.set_ylabel('Variants', fontsize=20)
        ax.set_xlabel('Chromosome', fontsize=20)
        ax.set_title('Variants Identified Across Chromosomes', fontsize=25)
        fig.save('Variants_per_chromosome.png', transparent=False, dpi=80,
                 bbox_inches='tight')
        return ax
