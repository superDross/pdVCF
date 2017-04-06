''' Manipulate Multi-Sample VCF files in Pandas (Python3). 
'''
import pandas as pd
import numpy as np
import re
import sys
import os
import collections

# path to this %%file
if sys.platform == "win32":
    file_path = os.path.dirname(os.path.abspath("__file__"))+"\\"
else:
    file_path = os.path.dirname(os.path.abspath("__file__"))+"/"

pd.set_option('display.max_columns', 500)
pd.set_option('display.max_colwidth', -1)


def vcf2dataframe(filename, genotype_level=True, info_level=True, UID=False):
    '''Open a VCF file and returns a MultiIndex pandas.DataFrame.
    Args:
        filename: vcf file to be converted to a dataframe
        genotype_level: place the genotype information into a second level column index
        info_level: place the info IDs into a second level column index
        UID: rename index to a unique variant identifier
    Notes:
        having any of these variables set to True will result
        in the DataFrame being generated very slowly. This is
        especially true for the UID variable.

        INFO DP field renamed to DEPTH
    '''
    if filename.endswith(".gz"):
        raise IOError("pdVCF does not support compressed VCF files.")

    # get INFO fields and Headers as lists
    VCF_HEADER = get_vcf_header(filename)
    INFO_FIELDS = get_info_fields(filename)

    # Count how many comment lines should be skipped.
    comments = count_comments(filename)

    # Return a simple dataframe representative of the VCF data.
    df = pd.read_table(filename, skiprows=comments,
                       names=VCF_HEADER, usecols=range(len(VCF_HEADER)))

    if genotype_level:
        df = get_genotype_data(df)

    if info_level:
        df = get_info_data(df, INFO_FIELDS)

    if UID:
        df = index2UID(df)
    
    # replace empty cells and drop FORMAT field
    df = df.replace("", np.nan)
    
    return df


def get_vcf_header(filename):
    ''' Get all header names from a given VCF file and return as a list.
    '''
    with open(filename) as input_file:
        row = [x for x in input_file if not x.startswith('##')] # skip unwanted headers
        head = next(iter(row))    # generator to deal with the header line only.
        split_head = [re.sub(r'#|\n', '', x) for x in head.split("\t")]
        return split_head


def get_info_fields(filename):
    ''' Get all ID names in the given VCFs INFO field and return as a list.
    '''
    with open(filename) as input_file:
        row = [x for x in input_file if x.startswith('##INFO')]
        info_fields = [x[11:].split(',')[0] for x in row]
        info_fields = [x.replace("DP", "DEPTH") for x in info_fields]
        return info_fields



def count_comments(filename):
    ''' Count all lines in a given VCF file starting with #.
    '''
    comments = 0
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                comments += 1

            else:
                break

    return comments



def replace_series_strings(df, col, dic, substring):
    ''' Replace the the keys with the items of the given
        dictionary for all strings or substrings in a
        given column
    Args:
        col: column name to replace strings
        dic: dictionary where the key is the string to replace with the item
        substrings: search and replace for either substrings (True) or exact strings (False)
    Returns:
        dataframe with the given column having all the
        entries identified as the key in the given dict
        replaced with the item in said dict
    '''
    if not isinstance(substring, bool):
        raise TypeError("substring argument must equal True or False")

    for string, correction in dic.items():
        if substring is True:
            df[col] = df[col].str.replace(string, correction)
        elif substring is False:
            df[col] = df[col].replace(string, correction, regex=True)

    return df



def get_genotype_data(df):
    ''' Give each sample column a second level column for every field
        detailed in the FORMAT column and return as a MultIndex dataframe.
    Args:
        df: DataFrame deriving from a VCF via vcf2dataframe()
    '''
    # contain the variant columns and the sample names in seperate lists
    normal = list(df.iloc[:, :9].columns)
    samples = list(df.iloc[:, 9:].columns)
    form = df['FORMAT'].str.split(":")[0]

    # These columns remain the same
    remain = pd.DataFrame(data=df[normal].values,
                          columns=pd.MultiIndex.from_tuples(
                            [(x, '') for x in normal] ))
    

    # list of dataframes where every sample has sub columns for each genotype info
    sams = [pd.DataFrame(data=[
                                # if all samples gentopes are './.' this fills it with 0
                                [x[0], '0', '0', '0', '0'] if len(x) == 1 else x
                                for x in list(df[col].str.split(':').dropna())
                              ],

                         columns=pd.MultiIndex.from_product([ [col], form ]))

            for col in samples]
    
    # add allele balance to sample genotype information
    sams = [calc_AB(sam) for sam in sams]
    
    # concat all dfs in the list
    df2 = pd.concat([remain] + sams, axis=1)

    return df2



def calc_AB(vcf):
    ''' Calculate allele balance for all samples in a given 
        pdVCF. Also converts DP & GQ to numeric type.
    
    Args:
        vcf: pdVCF with genotype information extracted
        
    Notes:
        ONLY WORKS FOR BIALLELIC VARIANTS
    '''
    sam = vcf.columns.levels[0][0]
    vcf[sam,'DP'] = pd.to_numeric(vcf[sam,'DP'].str.replace('.', '0')) # bcftools places '.' in empty fields
    vcf[sam,'GQ'] = pd.to_numeric(vcf[sam,'GQ'].str.replace('.', '0'))
    AD = vcf.xs('AD', level=1, axis=1).unstack().str.split(",", n=2)
    DP = vcf.xs('DP', level=1, axis=1).unstack()
    AB = round(pd.to_numeric(AD.str[1]) / pd.to_numeric(DP), 2)
    vcf[sam, 'AB'] = AB.tolist()
            
    return vcf



def get_info_data(df, info_fields):
    ''' Transform the INFO IDs into second level column indexes and return
        the df as a MultiIndex dataframe.
    Args:
        df: DataFrame deriving from a VCF via vcf2dataframe()
        info_fields: a list of all the INFO IDs in the given df
    '''
    # Alter Info field for some variables that don't work well
    df['INFO'] = df['INFO'].str.replace(";DB",";DB=1")
    df['INFO'] = df['INFO'].str.replace(";STR",";STR=1")
    df['INFO'] = df['INFO'].str.replace("DP", "DEPTH")

    # identify Info fields not present in each row and fill them with a 0
    for name in info_fields:
        if name == info_fields[0]:
            name = "{}=".format(name)
        else:
            name = ";{}=".format(name)

        not_present = df['INFO'][~df.INFO.str.contains(name)].add("{}0".format(name))
        present = df['INFO'][df.INFO.str.contains(name)]
        df['INFO'] = not_present.append(present).sort_index()

    # reorder INFO fields so they are are all in the same order
    df['INFO'] = df['INFO'].apply(lambda x: ';'.join(elem for elem in sorted(x.split(";"))))

    # remove all info_field names from the info values, starting with the info field with the longest name first
    unwanted = info_fields + ['=']
    unwanted.sort(key=len, reverse=True)
    remove = collections.OrderedDict([(x, '') for x in unwanted])
    df = replace_series_strings(df, col='INFO', dic=remove, substring=True)

    # create a new multi-index df containing only the info fields with the IDs as the second level
    info = pd.DataFrame(data=list(df['INFO'].str.split(';')),
                        columns=pd.MultiIndex.from_product([ ['INFO'], info_fields]))

    if not isinstance(df.columns, pd.MultiIndex):
        # create another multi-index df without the info fields where the second level is nothing
        df = pd.DataFrame(data=df.drop('INFO', axis=1, level=0).values,
                              columns=pd.MultiIndex.from_tuples(
                                [(x, '') for x in list(df.drop('INFO', axis=1, level=0).columns)] ))

    else:
        df = df.drop('INFO', level=0, axis=1)

    variant = df.iloc[:, :8]
    samples = df.iloc[:, 8:]

    # replace the info fields in the original df with the multi-index df created above
    final_df = pd.concat([variant] + [info] + [samples], axis=1)

    # MQ and MQ0 are in the wrong order so name swapping is required
    if 'MQ0' in info_fields and 'MQ' in info_fields:
        final_df = final_df.rename(columns={'MQ0': 'TEMP', 'MQ': 'MQ0'})
        final_df = final_df.rename(columns={'TEMP': 'MQ'})

    return final_df



def index2UID(df):
    ''' Replace the index with a unique variant identifier.
    '''
    if isinstance(df.columns, pd.MultiIndex):
        UID = df.apply(lambda x: "{}:{}-{}/{}".format(x['CHROM'][0], x['POS'][0],
                                                      x['REF'][0], x['ALT'][0]), axis=1)
    else:
        UID = df.apply(lambda x: "{}:{}-{}/{}".format(x['CHROM'], x['POS'],
                                                      x['REF'], x['ALT']), axis=1)

    df['UID'] = UID

    #if df['UID'].value_counts()[0] > 1:
    #    raise ValueError("The UID is not unique.")
        
    # remove UID column (UID now only accessible via index)  
    df = df.drop('UID', axis=1)


    return df.rename(UID)

