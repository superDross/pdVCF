# pDVCF
Manipulate a VCF as a Pandas DataFrame

# Example Usage
  from pdVCF.pdVCF import VCF
  
  # create a VCF object 
  sample = VCF("sample.vcf")

  # filter the vcf for variants with a minimum DP & GQ of 50 & 30 respectively across all samples in the VCF 
  sample.filter_genotype(minDP=50, minGQ=30)

  # filter out variants that lie within the below positions
  sample.positions(['1:2235000-2235600', '20:3457387'], include=False)

  # filter out multiallelic variants
  sample.biallelic()

  # print the filtered vcf
  sample.vcf
