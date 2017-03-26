# pDVCF
Manipulate a VCF as a Pandas DataFrame

# Example
  from pdVCF.pdVCF import VCF, FilterVCF
  
  sample = FilterVCF("sample.vcf")
  sample.filter_genotype(minDP=50, minGQ=30)
  sample.positions('1:2235000-2235600', include=False)
  sample.biallelic

  sample.vcf
