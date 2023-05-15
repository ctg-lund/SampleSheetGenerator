import pandas as pd
import re

class singleCellSheet():
    def __init__(self, data_csv, flexfile):
        # attempt to open data with pandas
        self.data = pd.read_csv(data_csv)
        self.flexfile = flexfile
        self.write_settings()
        self.write_header()
        self.parse_data()
        # index kits
        self.index_kits = {
            'NN': pd.read_csv('data/Dual_Index_Kit_NN_Set_A.csv'),
            'NT': pd.read_csv('data/Dual_Index_Kit_NT_Set_A.csv'),
            'TT': pd.read_csv('data/Dual_Index_Kit_TT_Set_A.csv'),
            'TS': pd.read_csv('data/Dual_Index_Kit_TS_Set_A.csv'),
            'TotalSeq': pd.read_csv('data/TotalSeq_A_Dual_Index_Primer_Tables.csv'),
            }
        
        self.parse_indeces()

        self.write_flex()

        self.write_data()

        self.write_10X()

        self.join_headers()

        print(self.flex_header)


    def parse_indeces(self):
        for counter, row in enumerate(self.data.itertuples()):
            for index_kit in self.index_kits:
                if row.index in self.index_kits[index_kit]['index_name'].tolist() and index_kit != 'TotalSeq':
                    self.data.loc[counter, 'index'] = self.index_kits[index_kit].loc[self.index_kits[index_kit].index_name == row.index, 'index(i7)'].values[0]
                    self.data.loc[counter, 'index2'] = self.index_kits[index_kit].loc[self.index_kits[index_kit].index_name == row.index, 'index2_workflow_b(i5)'].values[0]
                elif row.index in self.index_kits[index_kit]['index_name'].tolist():
                    print(row)
                    self.data.loc[counter, 'index'] = self.index_kits[index_kit].loc[self.index_kits[index_kit].index_name == row.index, 'index_sequence'].values[0]
                    self.data.loc[counter, 'index2'] = self.index_kits[index_kit].loc[self.index_kits[index_kit].index_name == row.index2, 'index_sequence'].values[0]

    def parse_data(self):
        # Check if all Sample_IDs are unique
        if len(self.data['Sample_ID'].unique()) != len(self.data['Sample_ID']):
            raise Exception('Sample_IDs are not unique!')
        # Check if all the indeces are unique
        indeces = list()
        for row in self.data.itertuples():
            indeces.append((row.index, row.index2))
        if len(indeces) != len(set(indeces)):
            raise Exception('Indeces are not unique!')

    def write_flex(self):
        if self.flexfile is not None:
            flex_columns = ['sample_id','probe_barcode_ids', 'Sample_Project']
            self.flex_header = '[Flex_Config]\n'
            self.flex_header = self.flex_header + self.flexfile[flex_columns].to_csv(index=False)
        else:
            self.flex_header = ''
        
    def write_data(self):
        data_columns = ['Sample_ID', 'index', 'index2','Sample_Project']
        self.data_header = "[BCLConvert_Data]\n"
        self.data_header += self.data[data_columns].to_csv(index=False)

    def write_header(self):
        self.header = '[Header]\n'
        self.header += 'FileFormatVersion,2\n'
    
    def write_settings(self):
        self.settings = '[BCLConvert_Settings]\n'
        self.settings += 'CreateFastqForIndexReads,0\n'
    def write_adt(self):
        self.adt_header = ''

    def write_10X(self):
        tenx_columns = ['Sample_ID','Sample_Project', 'Sample_Species', 
                        'pipeline', 'agg', 'force', 'test', 
                        'hto', 'libtype', 'sample_pair', 'nuclei'
                        ]
        tenx_columns = [x for x in tenx_columns if x in self.data.columns]
        self.tenx_header = '[10X_Data]\n'
        self.tenx_header += self.data[tenx_columns].to_csv(index=False)

    def join_headers(self):
        self.data = self.header + self.settings + self.data_header + self.tenx_header + self.flex_header



class illuminav2():
    def __init__(self, data_csv):
        # attempt to open data with pandas
        self.data = pd.read_csv(data_csv)
        self.indexes = pd.DataFrame()
        self.index_kit = ''
        self.commas = self.data.shape[1] + 1
        self.read1cycles = 151
        self.adapters = [
            'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 
            'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
            ]
        self.pipeline = 'SeqOnly'
        # the flowcell is located in RunInfo.xml
        #self.flowcell = 'H7Y2VBCXY'
        # flag for removal, track in database!
        self.lab_worker = 'John Doe'
        self.bnf_worker = 'Jane Doe'
        self.sequencer = 'NovaSeq'
        # always validate the input data!
        self.validate()

    def validate(self):
        # we should only have these columns:
        # 'Lane', 'Sample_ID', 'Sample_Project'
        # all other columns need to be removed
        # if these columns do not exist we fail
        valid_columns = ['Lane', 'Sample_ID', 'Sample_Project']
        for col in self.data.columns:
            if col not in valid_columns:
                self.data.drop(col, axis=1, inplace=True)
        # see that the number of columns is now 3
        if self.data.shape[1] < 3:
            raise Exception('No valid columns found in input data!')
    
    def get_indexes(self, index_kit, RC=True, reference='data/index_table.csv'):
        """
        Get the indexes for the given index kit.
        The index table is a csv file with the following columns:
        Index_Adapters, Index, Index2_forward_read, Index2_reverse_complement
        See: https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800
        """
        # read in the index table
        index_table = pd.read_csv(reference)
        # get the indexes for the given kit
        select = index_table.loc[index_table['Index_Adapters'] == index_kit]
        # only extract columns Index and Index2_reverse_complement or Index2_forward_read
        # depending on whether i2_rev_trans is True or False
        if RC:
            indexes = select.loc[:,['Index', 'Index2_reverse_complement']]
        else:
            indexes = select.loc[:,['Index', 'Index2_forward_read']]
        self.indexes = indexes
        
        # reshape indexes by removing rows until there are as many rows as in self.data
        self.indexes.drop(self.indexes.index[self.data.shape[0]:], inplace=True)
       
        # reset the index
        self.indexes.reset_index(drop=True, inplace=True)
        # add the indexes to the data
        print(indexes)
        self.data['index'] = self.indexes.loc[:,'Index']
        # variable column name!
        self.data['index2'] = self.indexes.iloc[:,1]
        print(self.data)

    
    def make_full_string(self):
        dstr = ','.join(re.split(r'[ \t]+', (self.data.to_string(index=False))))
        # remove leading commas
        dstr = re.sub(r'\n,','\n',dstr)
        # remove first comma
        dstr = re.sub(r'^,','\n',dstr)
        # trim whitespace
        dstr = dstr.strip()
        self.string = f"""[Header],,
FileFormatVersion,2,
,,
[Reads],,
Read1Cycles,{self.read1cycles},
,,
[Yggdrasil_Settings],,
Pipeline,{self.pipeline},,
LabWorker,{self.lab_worker},,
BNFWorker,{self.bnf_worker},,
,,
[BCLConvert_Settings],,
AdapterRead1,{self.adapters[0]},,
AdapterRead2,{self.adapters[1]},,
,,
[BCLConvert_Data],,
{dstr}
,,
"""

    def write_to_file(self, file):
        self.make_full_string() 
        with open(file, 'wb') as f:  
            f.write(self.string.encode('ascii', 'ignore'))