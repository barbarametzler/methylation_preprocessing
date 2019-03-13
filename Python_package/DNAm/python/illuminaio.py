# -*- coding: utf-8 -*-
# For license information, see LICENSE.TXT

#load requirements
import pandas as pd
import sys
import os
import re
import glob

def list_idat(path):
    file_list = [os.path.basename(x) for x in glob.glob(path + "/*.idat")]

    df_ = pd.DataFrame()

    regex = "(([0-9]+)_([RC0-9]+))_(Grn|Red)\\.idat"
    list1 = []

    for file in file_list:
        x = re.split(regex, file)
        list1.append(x)

    df_ = df_.append(pd.DataFrame(list1))
    df_.columns = ['id', 'sample.id', 'chip', 'position', 'channel', 'empty']
    df_ = df_.drop('empty', 1)

    return df_