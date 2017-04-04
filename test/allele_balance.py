from pdVCF.pdVCF import VCF
import copy

f = '/home/david/projects/scripts_cluster/projects/comparing_fluidigm_variants/vcfs/var.both.taadUkJan2017.filters.vcf'
vcf = VCF(f)


#AB_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1]
AB_list = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

for ab in AB_list:
    v = copy.deepcopy(vcf)

    ab_min = "AB > {}".format( max(0, ab-0.1))
    ab_max = "AB <= {}".format(ab)

    v.filter_vcf([ab_min, ab_max, 'DP >= 50', 'GQ >= 30', 'GT != 0/0', 'GT != ./.'], op='&', how='any')
    not_in_DB = len(v.vcf[v.vcf['ID'] == '.'])

    if v.vcf.shape[0] > 0:
        percent = not_in_DB / v.vcf.shape[0] * 100
        print("{} ({}/{}) of variants are not in dbSNP for AB between {}-{}".format(percent, not_in_DB, v.vcf.shape[0], max(0, round(ab-0.1, 1)), ab))
    else:
        print("0 variants after filtering for AB between {}-{}".format(ab_min, ab_max))


