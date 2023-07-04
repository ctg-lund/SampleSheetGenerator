from samplesheet import singleCellSheet
import pandas as pd

def read_samplesheets():
    samplesheet = pd.read_csv('data/singlecell/samplesheet.csv')
    feature_ref = pd.read_csv('data/singlecell/feature_ref.csv')
    flex_config = pd.read_csv('data/singlecell/flex_config.csv')
    return samplesheet, feature_ref, flex_config

def test():
    ss, fr, fc = read_samplesheets()
    samplesheet = singleCellSheet(ss, fc, fr, False) 