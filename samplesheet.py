import pandas as pd
import re

class singleCellSheet():
    def __init__(self, data_csv, flexfile, feature_ref):
        # attempt to open data with pandas
        self.data = pd.read_csv(data_csv)
        self.flexfile = flexfile
        self.feature_ref = feature_ref
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

        self.write_feature_ref()

        self.join_headers()


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
            flex_columns = ['sample_id','probe_barcode_ids', 'Sample_Source']
            self.flex_header = '[10X_Flex_Settings]\n'
            self.flex_header = self.flex_header + self.flexfile[flex_columns].to_csv(index=False)
        else:
            self.flex_header = '[10X_Flex_Settings]\n'
        
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

    def write_feature_ref(self):
        if self.feature_ref is not None:
            self.feature_ref_header = '[Feature_Reference]\n'
            self.feature_ref_header += self.feature_ref.to_csv(index=False)
        else:
            self.feature_ref_header = ''

    def write_10X(self):
        tenx_columns = ['Sample_ID','Sample_Project', 'Sample_Species', 
                        'pipeline', 'agg', 'force', 'vdj',
                        'sample_pair', 'nuclei', 'library'
                        ]
        self.chemistries2pipelines()
        tenx_columns = [x for x in tenx_columns if x in self.data.columns]
        self.tenx_header = '[10X_Data]\n'
        self.tenx_header += self.data[tenx_columns].to_csv(index=False)


    def chemistries2pipelines(self):
        chemistries_to_pipelines = {
            '3GEX': 'scrna-10x',
            '3CMO': 'scmulti-10x',
            '3ADT': 'scciteseq-10x',
            '3HTO': 'scciteseq-10x',
            '3CRISPR': 'scciteseq-10x',
            '5FB': 'seqonly',
            '5TCR': 'scmulti-10x',
            '5BCR': 'scmulti-10x',
            'ATAC': 'scatac-10x',
            'FLEX': 'scflex-10x',
            'VISIUM': 'scvisium-10x'
        }
        chemistries_to_libraries = {
            '3GEX': 'gex',
            '3CMO': 'cmo',
            '3ADT': 'adt',
            '3HTO': 'hto',
            '3CRISPR': 'crispr',
            '5FB': 'fb',
            '5TCR': 'tcr',
            '5BCR': 'bcr',
            'ATAC': 'atac',
            'FLEX': 'flex',
            'VISIUM': 'visium'
        }
        # Add pipeline and library columns
        self.data['pipeline'] = self.data['chemistry'].map(chemistries_to_pipelines)
        self.data['library'] = self.data['chemistry'].map(chemistries_to_libraries)
        self.data['vdj'] = ['n' for _ in range(len(self.data))]
        if any(['5BCR'==_ for _ in self.data['chemistry']]) or any(['5TCR'==_ for _ in self.data['chemistry']]):
            vdj_dictionary = {'5BCR': 'bcr', '5TCR': 'tcr', 'n': 'n'}
            self.data['vdj'] = self.data['chemistry'].apply(lambda x: vdj_dictionary[x] if x in vdj_dictionary else 'n')
        else:
            self.data['vdj'] = None
        
        # Checks if the add on libraries are correctly defined
        # Each sample should have a matching sample pair
        # Each sample pair should have matching pipelines
        # Its not beautiful but it's honest work
        for row in self.data.itertuples():
            if row.pipeline in ('scmulti-10x', 'scciteseq-10x'):
                if 'sample_pair' not in self.data.columns:
                    raise Exception('sample_pair column not defined in samplesheet file when using scmulti or scciteseq for sample: ' + row.Sample_ID)
                if row.sample_pair == 'n':
                    raise Exception('Matching sample pair is not defined for sample: ' + row.Sample_ID)
                else:
                    projectDf = self.data[self.data['Sample_Project'] == row.Sample_Project]
                    matching_sample_pairs = projectDf[projectDf['sample_pair'] == row.sample_pair]
                    if len(matching_sample_pairs) <= 1:
                        raise Exception('Matching sample pair is not defined for sample: ' + row.Sample_ID)
                    pipeline = (lambda x: 'scmulti-10x' if 'scmulti-10x' in x else 'scciteseq-10x')(matching_sample_pairs.pipeline.values)
                    self.data.at[row.Index, 'pipeline'] = pipeline
                    matching_sample_pipelines = projectDf[projectDf['sample_pair'] == row.sample_pair].pipeline.values
            elif 'sample_pair' in self.data.columns and row.sample_pair != 'n':
                projectDf = self.data[self.data['Sample_Project'] == row.Sample_Project]
                matching_sample_pipelines = projectDf[projectDf['sample_pair'] == row.sample_pair].pipeline.values
                if 'scmulti-10x' in matching_sample_pipelines:
                    matching_sample_pipeline = 'scmulti-10x'
                elif 'scciteseq-10x' in matching_sample_pipelines:
                    matching_sample_pipeline = 'scciteseq-10x'
                self.data.at[row.Index, 'pipeline'] = matching_sample_pipeline
                
        self.data.fillna('n', inplace=True)

    def join_headers(self):
        self.data = self.header + self.settings + self.data_header + self.tenx_header + self.flex_header + self.feature_ref_header



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