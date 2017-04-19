# pdVCF
Manipulate a VCF file as a Multi-Indexed Pandas DataFrame

## Filtering
Creating a VCF object from a large vcf file (>250Mb) will consume considerable memory resources e.g. a 287Mb vcf file consumes 6Gb of memory. It is therefore advised to perform any major filtering prior to initialisation. A secondary filtering stage can be performed after initialisation if complex filtering is required. 

One of the strengths of filtering via the filtering_vcf() method is its shear flexibility. For example, if one wishes to filter for variants where at least one sample in a given multi-sample vcf has an allele balance between 0.2-0.4, genotype depth above 50, alternative allele depth above 30, genotype quality above 30 and homozygous alternative genotype:

```python3
from pdVCF.pdVCF import VCF
example = VCF('example.vcf')
example.filter_vcf(['AB => 0.2', 'AB <= 0.4', 'DP > 50', 'AD[1] > 30', 'GQ > 30', 'GT = 1/1'], op='&', how='any')
```

## Install
To install and test pdVCF:
```python3
git clone https://github.com/superDross/pdVCF
cd pdVCF/
python3 pdVCF --test
```

### Example Usage
```python3
from pdVCF.pdVCF import VCF

# create a VCF object 
sample = VCF("sample.vcf")

# filter the vcf for variants with a minimum DP of 50 & GQ of 30 across all samples in the VCF object
sample.filter_vcf(['DP >= 50', 'GQ >= 30'], op='&', how='all')

# further filter for variants where any single sample has an allele balance (alt/(alt+ref)) between 0-0.3 
sample.filter_vcf(['AB > 0', 'AB < 0.3'], op="&", how='any')

# further filter for variants where any sample has a GT of 1/1 or 0/1
sample.filter_vcf(['GT == 1/1', 'GT == 0/1'], op="|", how="any")

# only keep variants that have an Allele Frequency above 0.3 for the first alternative allele called in a multi-allelic variant.
sample.filter_vcf(['AF[0] > 0.3'])

# keep A/G variants
sample.filter_vcf(['REF = A', 'ALT = G'])

# remove variants that lie within the below positions
sample.positions(['1:2235000-2235600', '20:3457387'], include=False)

# remove multiallelic variants
sample.biallelic()

# print the filtered vcf
sample.vcf
```

## To Do
- Functionality
  - Implement proper handling of multiallelic variants. Which affects: 
    - Genotype fields: AD, AB
    - Info fields: AF, AC etc.

  - Plotting features 
