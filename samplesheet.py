import pandas as pd
import re


class singleCellSheet:
    def __init__(self, data_csv, flexfile, feature_ref, singleindex: bool):
        # attempt to open data with pandas
        self.singleindex = singleindex
        self.dataDf = data_csv
        self.flexfile = flexfile
        self.feature_ref = feature_ref
        self.write_settings()
        self.write_header()
        self.parse_data()
        # index kits
        self.dual_index_kits = {
            "NN": pd.read_csv("data/Dual_Index_Kit_NN_Set_A.csv"),
            "NT": pd.read_csv("data/Dual_Index_Kit_NT_Set_A.csv"),
            "TT": pd.read_csv("data/Dual_Index_Kit_TT_Set_A.csv"),
            "TS": pd.read_csv("data/Dual_Index_Kit_TS_Set_A.csv"),
            "TotalSeq": pd.read_csv("data/TotalSeq_A_Dual_Index_Primer_Tables.csv"),
        }
        self.single_index_kits = {"N": pd.read_csv("data/Single_Index_Kit_N_Set_A.csv")}

        self.parse_indeces()

        self.write_flex()

        self.write_data()

        self.write_10X()

        self.write_feature_ref()

        self.parse_feature_ref()

        self.join_headers()

    def parse_indeces(self):
        for counter, row in enumerate(self.dataDf.itertuples()):
            if self.singleindex:
                for index_kit in self.single_index_kits:
                    if (
                        row.index
                        in self.single_index_kits[index_kit]["index_name"].tolist()
                    ):
                        # Get the matching row from the index kit
                        _index = self.single_index_kits[index_kit].loc[
                            self.single_index_kits[index_kit].index_name == row.index
                        ]
                        d = {
                            "Sample_ID": [row.Sample_ID for _ in range(4)],
                            "index": [_index[f"index{i}"].iloc[0] for i in range(1, 5)],
                            "Sample_Project": [row.Sample_Project for _ in range(4)],
                        }
                        newDf = pd.merge(
                            pd.DataFrame(d),
                            self.dataDf[self.dataDf["Sample_ID"] == row.Sample_ID].drop(
                                columns=["index", "Sample_Project"]
                            ),
                            on="Sample_ID",
                        )
                        # Remove the old row by the name in the Sample_ID column
                        self.dataDf = self.dataDf[
                            self.dataDf["Sample_ID"] != row.Sample_ID
                        ]
                        # Add the new rows
                        self.dataDf = pd.concat(
                            [self.dataDf, newDf], ignore_index=True, axis=0
                        )

            else:
                for index_kit in self.dual_index_kits:
                    if (
                        row.index
                        in self.dual_index_kits[index_kit]["index_name"].tolist()
                        and index_kit != "TotalSeq"
                    ):
                        self.dataDf.loc[counter, "index"] = (
                            self.dual_index_kits[index_kit]
                            .loc[
                                self.dual_index_kits[index_kit].index_name == row.index,
                                "index(i7)",
                            ]
                            .values[0]
                        )
                        self.dataDf.loc[counter, "index2"] = (
                            self.dual_index_kits[index_kit]
                            .loc[
                                self.dual_index_kits[index_kit].index_name == row.index,
                                "index2_workflow_b(i5)",
                            ]
                            .values[0]
                        )
                    elif (
                        row.index
                        in self.dual_index_kits[index_kit]["index_name"].tolist()
                    ):
                        self.dataDf.loc[counter, "index"] = (
                            self.dual_index_kits[index_kit]
                            .loc[
                                self.dual_index_kits[index_kit].index_name == row.index,
                                "index_sequence",
                            ]
                            .values[0]
                        )
                        self.dataDf.loc[counter, "index2"] = (
                            self.dual_index_kits[index_kit]
                            .loc[
                                self.dual_index_kits[index_kit].index_name
                                == row.index2,
                                "index_sequence",
                            ]
                            .values[0]
                        )

    def parse_data(self):
        # Check if all Sample_IDs are unique
        if len(self.dataDf["Sample_ID"].unique()) != len(self.dataDf["Sample_ID"]):
            raise Exception("Sample_IDs are not unique!")
        # Check if all the indeces are unique
        indeces = list()
        for row in self.dataDf.itertuples():
            if self.singleindex:
                indeces.append(row.index)
            else:
                indeces.append((row.index, row.index2))
        if len(indeces) != len(set(indeces)):
            raise Exception("Indeces are not unique!")

    def write_flex(self):
        if self.flexfile is not None:
            # Check if all the samples in the Sample_Source column exists in the self.dataDf Sample_ID column
            for row in self.flexfile.itertuples():
                if row.Sample_Source not in self.dataDf["Sample_ID"].tolist():
                    raise Exception(
                        "Sample_Source: "
                        + row.Sample_Source
                        + " does not exist in the Sample_ID column of the samplesheet!"
                    )
            flex_columns = ["sample_id", "probe_barcode_ids", "Sample_Source"]
            self.flex_header = "[FlexConfig_Data]\n"
            self.flex_header = self.flex_header + self.flexfile[flex_columns].to_csv(
                index=False
            )
        else:
            self.flex_header = "[FlexConfig_Data]\n"

    def write_data(self):
        if self.singleindex:
            data_columns = ["Sample_ID", "index", "Sample_Project"]
        else:
            data_columns = ['Sample_ID', 'index', 'index2','Sample_Project']
        self.dataDf_header = "[Data]\n"
        self.dataDf_header += self.dataDf[data_columns].to_csv(index=False)

    def write_header(self):
        self.header = "[Header]\n"
        self.header += "FileFormatVersion,1\n"

    def write_settings(self):
        self.settings = '[Settings]\n'

        if self.singleindex:
            self.settings += "CreateFastqForIndexReads,1\n"
            self.settings += "TrimUMI,0\n"
            self.settings += "OverrideCycles,Y50;I8;U24;Y49\n"
        else:
            self.settings += "CreateFastqForIndexReads,0\n"

    def write_adt(self):
        self.adt_header = ""

    def write_feature_ref(self):
        if self.feature_ref is not None:
            self.feature_ref_header = "[FeatureReference_Data]\n"
            self.feature_ref_header += self.feature_ref.to_csv(index=False)
        else:
            self.feature_ref_header = ""

    def parse_feature_ref(self):
        if self.feature_ref is not None:
            # Check if all the samples in the Sample_ID column exists in the self.dataDf Sample_ID column
            for row in self.feature_ref.itertuples():
                samples = row.Sample_Source.split("|")
                for sample in samples:
                    if sample not in self.dataDf["Sample_ID"].tolist():
                        raise Exception(
                            "Sample_Source: "
                            + sample
                            + " from the feature reference does not exist in the Sample_Source column of the samplesheet!"
                        )
            feature_ref_columns = [
                "id",
                "name",
                "read",
                "pattern",
                "sequence",
                "feature_type",
                "Sample_Source",
            ]
            self.feature_ref_header = "[FeatureReference_Data]\n"
            self.feature_ref_header = self.feature_ref_header + self.feature_ref[
                feature_ref_columns
            ].to_csv(index=False)

    def write_10X(self):
        tenx_columns = [
            "Sample_ID",
            "Sample_Project",
            "Sample_Species",
            "pipeline",
            "agg",
            "force",
            "vdj",
            "sample_pair",
            "nuclei",
            "libtype",
        ]
        self.parse_libraries_pipelines()
        tenx_columns = [x for x in tenx_columns if x in self.dataDf.columns]
        self.tenx_header = "[10X_Data]\n"
        self.tenx_header += (
            self.dataDf[tenx_columns].drop_duplicates().to_csv(index=False)
        )

    def parse_libraries_pipelines(self):
        # Parses the libraries and pipelines columns
        # Checks if the values are allowed
        # Check if the combinations are allowed
        allowed_libraries = [
            "gex",
            "cmo",
            "adt",
            "hto",
            "crispr",
            "tcr",
            "bcr",
            "atac",
            "flex",
            "visium",
            "flex",
        ]
        allowed_pipelines = [
            "scrna-10x",
            "scmulti-10x",
            "scatac-10x",
            "scciteseq-10x",
            "scarc-10x",
            "scvisium-10x",
            "scflex-10x",
        ]
        allowed_combinations = {
            "gex": ["scrna-10x", "scmulti-10x", "scarc-10x", "scciteseq-10x"],
            "cmo": ["scmulti-10x"],
            "adt": ["scmulti-10x", "scciteseq-10x"],
            "hto": ["scmulti-10x", "scciteseq-10x"],
            "crispr": ["scmulti-10x", "scciteseq-10x"],
            "tcr": ["scmulti-10x"],
            "bcr": ["scmulti-10x"],
            "atac": ["scatac-10x", "scarc-10x"],
            "flex": ["scflex-10x"],
            "visium": ["scvisium-10x"],
        }
        # Check if library column and pipeline column consist of allowed values
        if "pipeline" not in self.dataDf.columns:
            raise Exception("pipeline column not found in samplesheet!")
        for row in self.dataDf.itertuples():
            if "libtype" in row._fields and row.libtype not in allowed_libraries:
                raise Exception(
                    "libtype: "
                    + row.libtype
                    + " not allowed! Allowed libtypes are: "
                    + ", ".join(allowed_libraries)
                )
            if row.pipeline not in allowed_pipelines:
                raise Exception(
                    "pipeline: "
                    + row.pipeline
                    + " not allowed! Allowed pipelines are: "
                    + ", ".join(allowed_pipelines)
                )
            if (
                "libtype" in row._fields
                and row.pipeline not in allowed_combinations[row.libtype]
            ):
                raise Exception(
                    "libtype: "
                    + row.libtype
                    + " not allowed for pipeline: "
                    + row.pipeline
                    + "! Allowed libtypes for pipeline: "
                    + row.libtype
                    + " are: "
                    + ", ".join(allowed_combinations[row.libtype])
                )

    def parse_pairs(self):
        ...

    def join_headers(self):
        self.dataDf = (
            self.header
            + self.settings
            + self.dataDf_header
            + self.tenx_header
            + self.flex_header
            + self.feature_ref_header
        )


class illuminav2:
    def __init__(self, data_csv):
        # attempt to open data with pandas
        self.dataDf = pd.read_csv(data_csv)
        self.indexes = pd.DataFrame()
        self.index_kit = ""
        self.commas = self.dataDf.shape[1] + 1
        self.read1cycles = 151
        self.adapters = [
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        ]
        self.pipeline = "SeqOnly"
        # the flowcell is located in RunInfo.xml
        # self.flowcell = 'H7Y2VBCXY'
        # flag for removal, track in database!
        self.lab_worker = "John Doe"
        self.bnf_worker = "Jane Doe"
        self.sequencer = "NovaSeq"
        # always validate the input data!
        self.validate()

    def validate(self):
        # we should only have these columns:
        # 'Lane', 'Sample_ID', 'Sample_Project'
        # all other columns need to be removed
        # if these columns do not exist we fail
        valid_columns = ["Lane", "Sample_ID", "Sample_Project"]
        for col in self.dataDf.columns:
            if col not in valid_columns:
                self.dataDf.drop(col, axis=1, inplace=True)
        # see that the number of columns is now 3
        if self.dataDf.shape[1] < 3:
            raise Exception("No valid columns found in input data!")

    def get_indexes(self, index_kit, RC=True, reference="data/index_table.csv"):
        """
        Get the indexes for the given index kit.
        The index table is a csv file with the following columns:
        Index_Adapters, Index, Index2_forward_read, Index2_reverse_complement
        See: https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800
        """
        # read in the index table
        index_table = pd.read_csv(reference)
        # get the indexes for the given kit
        select = index_table.loc[index_table["Index_Adapters"] == index_kit]
        # only extract columns Index and Index2_reverse_complement or Index2_forward_read
        # depending on whether i2_rev_trans is True or False
        if RC:
            indexes = select.loc[:, ["Index", "Index2_reverse_complement"]]
        else:
            indexes = select.loc[:, ["Index", "Index2_forward_read"]]
        self.indexes = indexes

        # reshape indexes by removing rows until there are as many rows as in self.dataDf
        self.indexes.drop(self.indexes.index[self.dataDf.shape[0] :], inplace=True)

        # reset the index
        self.indexes.reset_index(drop=True, inplace=True)
        # add the indexes to the data
        print(indexes)
        self.dataDf["index"] = self.indexes.loc[:, "Index"]
        # variable column name!
        self.dataDf["index2"] = self.indexes.iloc[:, 1]
        print(self.dataDf)

    def make_full_string(self):
        dstr = ",".join(re.split(r"[ \t]+", (self.dataDf.to_string(index=False))))
        # remove leading commas
        dstr = re.sub(r"\n,", "\n", dstr)
        # remove first comma
        dstr = re.sub(r"^,", "\n", dstr)
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
        with open(file, "wb") as f:
            f.write(self.string.encode("ascii", "ignore"))


class pep2samplesheet:
    """
    The PEP is published here: <github link>
    Right now the PEP is included in this project until we decide on how to host it.
    This class will simply take a PEP and turn it into a samplesheet.
    As decided we will create Illumina v1 samplesheets going forward
    due to the restrictiveness of v2
    """
    def __init__(self, data_csv):
        # init data
        self.data = data_csv
        self.df = pd.read_csv(data_csv)
        # init patterns
        self.project_id_pattern = re.compile(r'^\d{4}_\d{3}$')
        self.sample_name_pattern = re.compile(r'^[0-9A-Za-z_-]+$')
        self.index_pattern = re.compile(r'^[ATCG]+$')
        # map PEP columns to samplesheet columns
        self.column_map : dict = {
            'project_id': 'Sample_Project',
            'sample_name': 'Sample_ID',
        }
        # validate the PEP
        self.validate()

    def validate(self):
        """
        Validate the PEP.
        Eventually this could be done with eidos
        but the regex patterns in the yaml can just
        be copied and validated directly
        without adding another dependency.

        This function will raise an exception if
        the PEP is invalid.

        Returns nothing.
        """
        
        if not self.df['project_id'].str.match(self.project_id_pattern).all():
            raise Exception('Invalid project_id found in PEP!')
        if not self.df['sample_name'].str.match(self.sample_name_pattern).all(): 
            raise Exception('Invalid sample_name found in PEP!')
        if not self.df['index'].str.match(self.index_pattern).all():
            raise Exception('Invalid index found in PEP!')
        
    
    def seq_only(self):
        """
        Assemble dataframe and illumina v1 static string
        return a string that can be written to a file
        """
        ss_1_string = """[Header]
        FileFormatVersion,1,
        [Data]
        """
        # rename columns
        df = self.df.rename(columns=self.column_map)
        # convert to string and append to ss_1_string
        ss_1_string += df.to_string(index=False)
        return ss_1_string
    
    def write_to_file(self, file):
        """
        Write the samplesheet to a file
        returns nothing
        """
        ss_1_string = self.seq_only()
        with open(file, 'wb') as f:
            f.write(ss_1_string.encode('ascii', 'ignore'))
    
