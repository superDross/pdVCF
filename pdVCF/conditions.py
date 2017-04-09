''' Collection of functions that create and combine conditons
'''
import pandas as pd
from functools import reduce
import operator
import re

ops = { 
    
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    "&": operator.and_,
    "|": operator.or_ 
}


def create_condition(vcf, string):
    ''' Return a condition from an string representing 
        desired condition.
        
    Args:
        vcf: VCF object
        string: condition in string form
        
    Returns:
        boolean dtype of the given data frame after 
        being tested with the given condition(s)
        
    Example:
        string = "DP >= 40"
        Assesses whether the given vcf has any DP columns 
        with value of or above 40
            
    '''
    field, op_sign, val = string.split(" ")

    # for string indexing fields with ',' in cells
    if '[' in field:
        index = int(field.split("[")[1].replace(']', ''))
        field = field.split("[")[0]

    # get the operator function that matches the op_sign
    op = ops[op_sign]
    # if not in INFO or genotype fields, assume its in first level 
    level = 1 if field in vcf.columns.levels[1] else 0
    # get the cross section of all columns with the given field
    all_fields = vcf.xs(field, axis=1, level=level)
    
    # determine if there is a comma in any of all_fields values
    comma = all_fields.astype(str).apply(lambda x: x.str.contains(',', na=False)).any().any()
    if comma and 'index' in locals():
        all_fields = all_fields.apply(lambda x: x.str.split(",").str[index])

    # convert digits to numeric type. replace needed for identifying floats
    if val.replace("0.", "").isdigit():   
        all_fields = all_fields.apply(pd.to_numeric) 
        all_fields = all_fields.fillna(0)
        val = float(val)

    # apply the operator function to the cross sectioned vcf against the given value
    cond = op(all_fields, val)
    
    return cond
    



def combine_conditions(cond_list, op):
    ''' Combine a list of conditions with the given
        operator.
        
    Args:
        cond_list: list of conditions
        op: operator to combine cond_list with e.g. "&"
    
    '''
    op_func = ops[op]
    combined = reduce(lambda x,y: op_func(x, y), cond_list)
    return combined

