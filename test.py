from samplesheet import singleCellSheet
import pandas as pd
import os
from os.path import join, normpath, isabs
def abspath(path):
    """Return an absolute path."""
    path = os.fspath(path)
    if not isabs(path):
        if isinstance(path, bytes):
            cwd = os.getcwdb()
        else:
            cwd = os.getcwd()
        path = join(cwd, path)
    return normpath(path)

def read_samplesheets(samplesheet:str, feature_ref:str, flex_config:str, ctg_samplesheet: str) -> tuple:
    samplesheet = pd.read_csv(abspath(samplesheet))
    feature_ref = pd.read_csv(abspath(feature_ref))
    flex_config = pd.read_csv(abspath(flex_config))
    with open(abspath(ctg_samplesheet)) as f:
        ctg_samplesheet = f.read()
    return samplesheet, feature_ref, flex_config, ctg_samplesheet

def test_samplesheet():
    ss, fr, fc, ctg_samplesheet = read_samplesheets('data/singlecell/SampleSheet.csv', 
                                                    'data/singlecell/feature_ref.csv', 
                                                    'data/singlecell/flex_config.csv',
                                                    'data/singlecell/CTG_SampleSheet.csv'
                                                    )
    samplesheet = singleCellSheet(ss, fc, fr, False) 
    assert samplesheet.dataDf == ctg_samplesheet