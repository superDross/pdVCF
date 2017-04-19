from pdVCF.pdVCF import VCF
import unittest
import os

here = os.path.dirname(os.path.realpath(__file__))+"/"
# 'Tested:' within the docstring details commands it was orginally tested against and which produced the same given answer as pdVCF filtering


class TestMandatory(unittest.TestCase):

    def test_pos_ranges(self):
        ''' Filter for variants between range

        Tested:
            SnpSift filter "( CHROM = '1') & ( POS > 2000 ) & ( POS < 2235601 )" testing2.vcf 
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['CHROM = 1', 'POS > 2000', 'POS < 2235601'])
        answer = ['1:2234385-C/T', '1:2235243-C/T', '1:2235501-A/GGT']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_loads(self):
        ''' Test numerous mandatory fields

        Tested:
            SnpSift filter "( FILTER != 'PASS' ) & ( ALT = 'G' ) & ( QUAL > 200 )" testing2.vcf
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['FILTER != PASS', 'ALT = G', 'QUAL > 200'])
        answer = ['1:2235792-A/G']

        self.assertEqual(v.vcf.index.tolist(), answer)


class TestFilterGenotype(unittest.TestCase):
    
    def test_gt_filter(self):
        ''' Simple genotype filter

        Tested:
            SnpSift filter "( GEN[*].GT == '1/1' ) & ( GEN[*].DP > 100 )" testing.vcf
        '''
        v = VCF(here+'vcfs/testing.vcf')
        v.filter_vcf(['GT = 1/1', 'DP > 100'], how='any')
        answer = ['1:2235901-A/G', '1:2235501-A/GGT', '1:2239901-GA/G']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_gt_big_file(self):
        ''' Same test as above but with a larger vcf file

        Tested:
            SnpSift filter "( GEN[*].GT == '1/1' ) & ( GEN[*].DP > 100 )" testing3.vcf | grep -v '^#' | wc -l
        '''
        v = VCF(here+'vcfs/testing3.vcf')
        v.filter_vcf(['GT = 1/1', 'DP > 100'], how='any')

        self.assertEqual(len(v.vcf), 257)


    def test_AB(self):
        ''' Filter AB < 0.09

        Tested:
            I cannot find a tool to filter VCFs by allele balance
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['AB > 0', 'AB < 0.09'])
        answer = ['1:2239999-A/G,T']

        self.assertEqual(answer, v.vcf.index.tolist())



class TestFilterInfo(unittest.TestCase):

    def test_info_depth(self):
        ''' Testing filtering the INFO DP field

        Tested:
            SnpSift filter "( DP > 100 )" testing.vcf
        '''
        v = VCF(here+'vcfs/testing.vcf')
        v.filter_vcf(['DEPTH > 100'])
        answer = ['1:2234385-C/T', '1:2239999-A/G,T']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_info_float_int(self):
        ''' Filter int and float fields

        Tested:
            SnpSift filter "( AF[0] < 0.1 ) & ( AC[0] > 2 )" testing2.vcf 
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['AF[0] < 0.1', 'AC[0] > 2'])
        answer = ['1:2235901-A/G']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_string(self):
        ''' Filter str fields

        Tested:
            SnpSift filter "( MEOW != 'ON' )" testing2.vcf 
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['MEOW != ON'])
        answer = ['1:2234385-C/T', '1:2235901-A/G', '1:2235501-A/GGT']

        self.assertEqual(v.vcf.index.tolist(), answer)


class TestFilterArgs(unittest.TestCase):

    def test_all_genotype(self):
        ''' Test filtering using the all flag in genotype fields

        Tested:
            SnpSift filter "( GEN[?].DP >= 50 ) & ( GEN[?].GQ >= 30 ) " testing.vcf
        '''
        v = VCF(here+'vcfs/testing.vcf')
        v.filter_vcf(['DP >= 50', 'GQ >= 30'], op='&', how='all')
        answer = ['1:2235792-A/G']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_any_filtering(self):
        ''' Test filtering using the any flag

        Tested:
            SnpSift filter "( GEN[*].DP >= 50 ) & ( GEN[*].GQ >= 30 ) " testing.vcf
        '''
        v = VCF(here+'vcfs/testing.vcf')
        v.filter_vcf(['DP >= 50', 'GQ >= 30'], op='&', how='any')
        answer = ['1:2235243-C/T', '1:2235792-A/G', '1:2235901-A/G', 
                  '1:2239999-A/G,T', '1:2235501-A/GGT', '1:2239901-GA/G']

        self.assertEqual(v.vcf.index.tolist(), answer)

    def test_or_operator(self):
        ''' Testing or operator

        Tested:
            SnpSift filter "( MEOW != 'ON' ) | ( LOLZ > 200 )" testing2.vcf
        '''
        v = VCF(here+'vcfs/testing2.vcf')
        v.filter_vcf(['MEOW != ON', 'LOLZ > 200'], op="|")
        answer = ['1:2234385-C/T', '1:2235901-A/G', '1:2239999-A/G,T', '1:2235501-A/GGT']

        self.assertEqual(v.vcf.index.tolist(), answer)


class TestUltimateFilter(unittest.TestCase):
    
    def test_ultimate(self):
        ''' Test limits of filtering

        Tested:
         SnpSift filter "( GEN[*].GT == '0/1' ) & ( GEN[*].DP >= 50 ) & ( GEN[*].GQ >= 30 ) & ( AC[0] > 20 ) & ( CHROM = '1' ) & ( POS > 2235893 ) & ( ID =~ 'rs' )" testing3.vcf
        '''
        m = VCF(here+'vcfs/testing3.vcf')
        m.filter_vcf(['GT = 0/1', 'DP >= 50', 'GQ >= 30'])
        m.filter_vcf(['AC[0] > 20'])
        m.filter_vcf(['CHROM = 1', 'POS > 2235893', 'ID != .'])

        answer = ['1:218519928-A/AAAAC', '1:218578726-ACTCT/A,ACTCTCT,ACT', 
                  '1:218607557-G/T'] 

        self.assertEqual(m.vcf.index.tolist(), answer)

if __name__ == '__main__':
    unittest.main(warnings='ignore')

