from samplesheet import singleCellSheet
import pandas as pd
from os import environ

def read_samplesheets():
    basepath = environ.get('GITHUB_WORKSPACE', '')
    samplesheet = pd.read_csv(f'{basepath}/data/singlecell/samplesheet.csv')
    feature_ref = pd.read_csv(f'{basepath}/data/singlecell/feature_ref.csv')
    flex_config = pd.read_csv(f'{basepath}/data/singlecell/flex_config.csv')
    with open('data/singlecell/ctg_samplesheet.csv') as f: 
        ctg_samplesheet = f.read()
    return samplesheet, feature_ref, flex_config, ctg_samplesheet

def test_samplesheet():
    ss, fr, fc, ctg_samplesheet = read_samplesheets()
    samplesheet = singleCellSheet(ss, fc, fr, False) 
    assert samplesheet.dataDf == ctg_samplesheet