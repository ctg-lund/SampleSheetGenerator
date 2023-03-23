import pandas as pd
import re

class singleCellSheet():
    def __init__(self, data_csv):
        # attempt to open data with pandas
        self.data = pd.read_csv(data_csv)
        # index kits
        self.index_kits = {
            'NN': pd.DataFrame.read_csv('data/Dual_Index_Kit_NN_Set_A.csv'),
            'NT': pd.DataFrame.read_csv('data/Dual_Index_Kit_NT_Set_A.csv'),
            'TT': pd.DataFrame.read_csv('data/Dual_Index_Kit_TT_Set_A.csv')
            }
    def parse_indeces():
        ...
    def write_data(self):
        ...
    def write_adt(self):
        ...
    def write_10X(self):
        ...
    def write_flex(self):
        ...
    def data_to_string(self):
        ...
class illuminav2():
    def __init__(self, data_csv):
        # attempt to open data with pandas
        self.data = pd.read_csv(data_csv)
        # defaults
        # commas = columns + 1 (should never be less than 5)
        self.commas = self.data.shape[1] + 1
        self.read1cycles = 151
        self.adapters = [
            'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 
            'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
            ]
        self.pipeline = 'SeqOnly'
        # the flowcell is located in RunInfo.xml
        #self.flowcell = 'H7Y2VBCXY'
        self.lab_worker = 'John Doe'
        self.bnf_worker = 'Jane Doe'
        # always validate the input data!
        self.validate()

    def validate(self):
        # we should only have these columns:
        # 'Lane'*, 'Sample_ID', 'Sample_Project', 'index', 'index2'
        # Lane is optional!
        # all other columns need to be removed
        # if these columns do not exist we fail
        valid_columns = ['Lane', 'Sample_ID', 'Sample_Project', 'index', 'index2']
        for col in self.data.columns:
            if col not in valid_columns:
                self.data.drop(col, axis=1, inplace=True)
        # see that the number of columns is now 4
        if self.data.shape[1] < 4:
            raise Exception('No valid columns found in input data!')

    def set_read1cycles(self, read1):
        self.read1cycles = int(read1) + 1 
    
    def set_adapters(self, a1, a2): 
        self.adapters = [a1, a2]
    
    def set_pipeline(self, pipeline):
        self.pipeline = pipeline
    
    def set_lab_worker(self, lab_worker):
        self.lab_worker = lab_worker
    
    def set_bnf_worker(self, bnf_worker):       
        self.bnf_worker = bnf_worker    
    
    def make_full_string(self):
        dstr = ','.join(re.split(r'[ \t]+', (self.data.to_string(index=False))))
        self.string = "[Header]" + "," * self.commas + "\n" +\
        "FileFormatVersion,2" + "," * (self.commas - 1) + "\n" +\
        "," * self.commas + "\n" +\
        "[Reads]" + "," * self.commas + "\n" +\
        "Read1Cycles," + str(self.read1cycles) + "," * (self.commas - 1) + "\n" +\
        "," * self.commas + "\n" +\
        "[Yggdrasil_Settings]" + "," * self.commas + "\n" +\
        "Pipeline" + "," + self.pipeline + "," * (self.commas - 1) + "\n" +\
        "LabWorker" + "," + self.lab_worker + "," * (self.commas - 1) + "\n" +\
        "BNFWorker" + "," + self.bnf_worker + "," * (self.commas - 1) + "\n" +\
        "," * self.commas + "\n" +\
        "[BCLConvert_Settings]" + "," * self.commas + "\n" +\
        "AdapterRead1" + "," + self.adapters[0] + "," * (self.commas - 1) + "\n" +\
        "AdapterRead2" + "," + self.adapters[1] + "," * (self.commas - 1) + "\n" +\
        "," * self.commas + "\n" +\
        "[BCLConvert_Data]" + "," * self.commas + "\n" +\
        f"{dstr.lstrip(',')}" + "\n" +\
        "," * self.commas + "\n"
        # eventually add Yggdrasil_Data
    
    def write_to_file(self, file):
        self.make_full_string() 
        with open(file, 'wb') as f:  
            f.write(self.string.encode('ascii', 'ignore')) 