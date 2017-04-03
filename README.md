# pdVCF
Manipulate a VCF file as a MultiIndexed Pandas DataFrame

### Example Usage
```python3
from pdVCF.pdVCF import VCF

# create a VCF object 
sample = VCF("sample.vcf")

# filter the vcf for variants with a minimum DP & GQ of 50 & 30 respectively across all samples in the VCF 
sample.filter_vcf(['DP >= 50', 'GQ >= 30'], op='&', how='all')

# further filter for variants where any single sample has an allele balance (alt/(alt+ref)) between 0-0.3 
sample.filter_vcf(['AB > 0', 'AB < 0.3'], op="&", how='any')

# further filter for variants where any sample has a GT of 1/1 or 0/1
sample.filter_vcf(['GT == 1/1', 'GT == 0/1'], op="|", how="any")

# only keep variants that have an Allele Frequency above 0.3
sample.filter_vcf(['AC > 0.3'])

# remove variants that lie within the below positions
sample.positions(['1:2235000-2235600', '20:3457387'], include=False)

# remove multiallelic variants
sample.biallelic()

# print the filtered vcf
sample.vcf
```

## To Do
- Testing
  - Test filtering against vcftools
  - Unit testing

- Functionality
  - Implement proper handling of multiallelic variants. Which affects: 
    - Genotype fields: AD, AB
    - Info fields: AF, AC etc.

  - Add more filtering options
  - Plotting features 
