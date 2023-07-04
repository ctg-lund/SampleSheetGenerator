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

def read_samplesheets():
    samplesheet = pd.read_csv(abspath('data/singlecell/samplesheet.csv'))
    feature_ref = pd.read_csv(abspath('data/singlecell/feature_ref.csv'))
    flex_config = pd.read_csv(abspath('data/singlecell/flex_config.csv'))
    with open(abspath('data/singlecell/ctg_samplesheet.csv')) as f: 
        ctg_samplesheet = f.read()
    return samplesheet, feature_ref, flex_config, ctg_samplesheet

def test_samplesheet():
    ss, fr, fc, ctg_samplesheet = read_samplesheets()
    samplesheet = singleCellSheet(ss, fc, fr, False) 
    assert samplesheet.dataDf == ctg_samplesheet