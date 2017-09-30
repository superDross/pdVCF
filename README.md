# pdVCF
Manipulate a [vcf](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file as a Multi-Indexed Pandas DataFrame

## Install
To install and test pdVCF:
```bash
git clone https://github.com/superDross/pdVCF
export PYTHONPATH=$PYTHONPATH:/path/to/pdVCF/
cd pdVCF/
python pdVCF --test
```

## Filtering
One of the strengths of filtering via the filtering_vcf() method is its shear flexibility. For example, if one wishes to filter for variants where at least one sample in a given multi-sample vcf has an allele balance between 0.3-0.6, alternative allele depth and genotype quality above 30 and homozygous alternative genotype:
```python3
from pdVCF.pdVCF import Vcf
example = Vcf('test/vcfs/testing2.vcf')
example.filter_vcf(['AB >= 0.3', 'AB <= 0.6', 'GT = 1/1', 'GQ > 30', 'AD[1] > 30'], op='&', how='any')
```
![](docs/filtered_vcf2.png?raw=true)

### Example Usage
```python3
from pdVCF.pdVCF import Vcf

# create a Vcf object 
sample = Vcf("sample.vcf")

# filter the vcf for variants with a minimum DP of 50 & GQ of 30 across all samples in the Vcf object
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
## Caveats
Creating a Vcf object from a large vcf file will consume considerable memory resources. It is therefore advised to perform any major filtering prior to initialisation using a tools such as [VcfFilerPy](https://github.com/superDross/VcfFilterPy) or [SnpSift](https://github.com/pcingola/SnpSift). A secondary filtering stage can be performed after initialisation if complex filtering is required. 

## To Do
- Optimisation
  - Use pd.eval() and pd.query() to reduce memory usage

- Functionality
  - Implement proper handling of multiallelic variants. Which affects: 
    - Genotype fields: AD, AB
    - Info fields: AF, AC etc.

  - Plotting features 
