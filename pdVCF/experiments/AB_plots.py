from pdVCF.pdVCF import VCF
import copy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


f = 'vcfs/vep.dbSNP.var.ug.both.yale.taad.filters.renamed.vcf'
vcf = VCF(f)
vcf_uk = VCF('vcfs/var.both.taadUkJan2017.filters.vcf')


def ab_stuff(vcf):
    ''' Return a nested list containing the number of variants 
        in a given allele balance range from a an AB list.

    '''
    AB_list = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

    data = []

    for ab in AB_list:
        v = copy.deepcopy(vcf)

        ab_min = "AB > {}".format( ab)
        ab_max = "AB <= {}".format(1)

        v.filter_vcf([ab_min, ab_max, 'DP >= 50', 'GQ >= 30', 'GT != 0/0', 'GT != ./.'], op='&', how='any')
        not_in_DB = len(v.vcf[v.vcf['ID'] == '.'])

        if v.vcf.shape[0] > 0:
            percent = not_in_DB / v.vcf.shape[0] * 100
            print("{} ({}/{}) of variants are not in dbSNP for AB between {}-{}".format(percent, not_in_DB, v.vcf.shape[0], 1, ab))
        else:
            print("0 variants after filtering for AB between {}-{}".format(ab_min, ab_max))
            pass

        data.append([ab, v.vcf.shape[0], not_in_DB, percent])




    df = pd.DataFrame(columns=['AB', 'No. Variants', 'Not in dbSNP', 'Percent'],
                      data=data)
    return df



# AB dataframe construction
df_yale = ab_stuff(vcf)
df_uk = ab_stuff(vcf_uk)

# combine the AB dataFrames
var = df_yale['No. Variants'] + df_uk['No. Variants']
no = df_yale['Not in dbSNP'] + df_uk['Not in dbSNP']
percent = no/var * 100
df = pd.concat([df_yale['AB'], var, no, percent.rename('Percent')], axis=1)




# CREATE THE PLOT
sns.set(style="whitegrid")

# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(10, 10))

# all vars
sns.set_color_codes("pastel")
sns.barplot(y="No. Variants", x="AB", data=df,
            label="Total", color="0.9")

# vars NOT in dbSNP
sns.set_color_codes("muted")
sns.barplot(y="Not in dbSNP", x="AB", data=df,
            label="Likely Artifacts", color="0.7")

ax.set(ylabel='Variants', xlabel='Allele Balance')
ax.legend(ncol=2, loc="upper right", frameon=True)


for x in range(len(df)):
    ax.text(x-0.3, ax.patches[x].get_height()+4, "{}%".format(round(df['Percent'].ix[x], 1)))

fig = ax.get_figure()
fig.savefig('AB Series Present in dbSNP.TEST.png')




