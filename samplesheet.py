import pandas as pd
import re


class singleCellSheet:
    def __init__(
        self,
        data_csv,
        flexfile,
        feature_ref,
        singleindex: bool,
        development_status: bool,
    ):
        # attempt to open data with pandas
        self.singleindex = singleindex
        self.development_status = development_status
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
            "TN": pd.read_csv("data/Dual_Index_Kit_TN_Set_A.csv"),
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
            data_columns = ["Sample_ID", "index", "index2", "Sample_Project"]
        self.dataDf_header = "[Data]\n"
        self.dataDf_header += self.dataDf[data_columns].to_csv(index=False)

    def write_header(self):
        self.header = "[Header]\n"
        self.header += "FileFormatVersion,1\n"
        if self.development_status:
            self.header += "DevelopmentProject,Yes\n"
        else:
            self.header += "DevelopmentProject,No\n"

    def write_settings(self):
        self.settings = "[Settings]\n"

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
            "cytaimage",
            "darkimage",
            "image",
            "slide",
            "slide_area",
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


def replace_missing_num_values_with_default(df, column, default) -> None:
    # check if column exists first
    if column in df.columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
        df[column].fillna(default, inplace=True)


class pep2samplesheet:
    """
    The PEP is published here: <github link>
    Right now the PEP is included in this project until we decide on how to host it.
    This class will simply take a PEP and turn it into a samplesheet.
    As decided we will create Illumina v1 samplesheets going forward
    due to the restrictiveness of v2
    """

    def __init__(self, data_csv, projects_csv):
        # init data
        self.data = data_csv
        self.df = pd.read_csv(data_csv)
        self.projects = pd.read_csv(projects_csv)
        self.lane_divider = False
        self.shared_flowcell = False
        self.single_index = False
        # check if single index
        self.check_single_index()
        # check lanes
        self.check_shared_flowcell()
        # lowercase
        self.lower_case_colnames()
        # check duplicate indexes
        self.check_duplicate_indexes()
        # Params
        self.sequencer = ""
        self.seqonly_project = "No"
        # keep empty for testing
        self.flowcell = ""
        self.dev_project = "No"
        # init patterns
        self.project_id_pattern = re.compile(r"^\d{4}_\d{3}$")
        self.sample_name_pattern = re.compile(r"^[0-9A-Za-z_-]+$")
        self.index_pattern = re.compile(r"^[ATCG]+$")
        # map PEP columns to samplesheet columns
        self.column_map: dict = {
            "project_id": "Sample_Project",
            "sample_name": "Sample_ID",
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
        # check for empty columns
        if self.df.isnull().values.any():
            raise Exception("Empty columns found in samples.csv!")
        # check for invaild project id, sample name and index
        if not self.df["project_id"].str.match(self.project_id_pattern).all():
            raise Exception("Invalid project_id found in samples.csv!")
        if not self.df["sample_name"].str.match(self.sample_name_pattern).all():
            raise Exception("Invalid sample_name found in samples.csv!")
        if not self.df["index"].str.match(self.index_pattern).all():
            raise Exception("Invalid index found in samples.csv!")
        # valid columns for samples.csv are:
        # project_id, sample_name, index, index2, reference, experiment, control, lane
        valid_columns = [
            "project_id",
            "sample_name",
            "index",
            "index2",
            "reference",
            "experiment",
            "control",
            "lane",
            "panel",
            "overridecycles",
            "barcodemismatchesindex1",
            "barcodemismatchesindex2",
        ]
        for col in self.df.columns:
            if col not in valid_columns:
                raise Exception(f"Invalid column found in samples.csv: {col}")
        # if there is a reference column, check that it contains valid
        if "reference" in self.df.columns:
            # values, accepted values are:
            valid_references = [
                "hg19",
                "hg38",
                "Homo sapiens",
                "Human",
                "Mouse",
                "Rat",
                "Mus musculus",
                "Rattus norvegicus",
                "Pig",
            ]
            if not self.df["reference"].isin(valid_references).all():
                raise Exception(
                    f"Invalid reference!\nValid references are:\{valid_references}"
                )
            # check that BarcodeMismatchesIndex1 and BarcodeMismatchesIndex2
            # default to 1 if not specified as a number
            replace_missing_num_values_with_default(
                df=self.df, column="BarcodeMismatchesIndex1", default=1
            )
            replace_missing_num_values_with_default(
                df=self.df, column="BarcodeMismatchesIndex2", default=1
            )

        # only check that columns exist
        # projects df needs at least a project_id and a fastq column
        if not "project_id" in self.projects.columns:
            raise Exception("project_id column not found in projects.csv!")
        if not "fastq" in self.projects.columns:
            raise Exception("fastq column not found in projects.csv!")

    def check_duplicate_indexes(self):
        """
        Illumina machines will start even if a sample has the same
        index pair as another in the same pool. This will cause
        the samples to be mixed up. This function will check for
        duplicate index pairs and raise an exception if any are found.

        If there is a lane divider, the index pairs are allowed to be
        the same across lanes, but not within.
        """
        # check for duplicate index pairs
        if self.lane_divider and not self.single_index:
            # in this case subset by lane
            # then look for duplicates in index + index2
            # report which samples and lanes contain duplicates
            for lane in self.df["lane"].unique():
                subset = self.df[self.df["lane"] == lane]
                duplicated_rows = subset[
                    subset["index"].str.cat(subset["index2"]).duplicated()
                ]
                if duplicated_rows.shape[0] > 0:
                    msg = f"""
                    <h4> Duplicates found in samples.csv!</h4>

                    The following samples have the same index pair in lane {lane}:
                    {duplicated_rows['sample_name'].to_string(index=False)}
                    
                    Table:
                    {duplicated_rows.to_html(index=False)}
                    """
                    raise Exception(msg)

        elif not self.single_index:
            # check for duplicate index pairs (index pair = index + index2)
            # report which samples contain duplicates
            duplicated_rows = self.df[
                self.df["index"].str.cat(self.df["index2"]).duplicated()
            ]
            if duplicated_rows.shape[0] > 0:
                # we can only print in one line
                # so we will print the sample names, project ids and index pairs
                msg = f"""
                <h4>Duplicates found in samples.csv!</h4>

                The following samples have the same index pair:
                {duplicated_rows['sample_name'].to_string(index=False)}
                
                Table:
                {duplicated_rows.to_html(index=False)}    
"""
                raise Exception(msg)

    def check_shared_flowcell(self):
        """
        Check if samples belong to a shared flowcell.
        A shared flowcell has more than one project id in
        the project_id column.
        Also check if there is a lane divider. There is
        a lane divider if there are more than one unique
        value in the lane column.
        """
        # check if there is more than one unique project id
        if len(self.df["project_id"].unique()) > 1:
            self.shared_flowcell = True
        # check if there is more than one unique lane
        # the lane column is optional
        if "lane" in self.df.columns:
            if len(self.df["lane"].unique()) > 1:
                self.lane_divider = True

    def check_single_index(self):
        """
        The input is single index if there is a index column
        and not an index2 column OR
        if there is an index column and an index2 column
        but index2 is empty.
        """
        if not "index2" in self.df.columns:
            self.single_index = True
        elif self.df["index2"].isna().all():
            self.single_index = True

    def lower_case_colnames(self):
        """
        Make all column names lowercase in self.df and self.projects
        """
        self.df.columns = self.df.columns.str.lower()
        self.projects.columns = self.projects.columns.str.lower()

    def make_ss(self):
        """
        Assemble dataframe and illumina v1 static string
        return a string that can be written to a file
        """
        # header
        ss_1_string = "[Header]\nFileFormatVersion,1,\n"
        # sequencer
        ss_1_string += f"Sequencer,{self.sequencer},\n"
        # dev project
        ss_1_string += f"DevelopmentProject,{self.dev_project},\n"
        # flowcell
        ss_1_string += f"Flowcell,{self.flowcell},\n"
        # data
        ss_1_string += "\n[Data]\n"
        # rename columns
        df = self.df.rename(columns=self.column_map)
        # headers
        ss_1_string += ",".join(df.columns) + "\n"
        # Convert the dataframe to a string
        str_df = df.apply(lambda x: ",".join(x.astype(str)), axis=1)
        # Join the rows with a newline
        str_df = "\n".join(str_df)
        # Append the string to the samplesheet
        ss_1_string += str_df
        # rename columns in projects.csv
        self.projects.rename(columns=self.column_map, inplace=True)
        # convert self.projects to a string
        str_projects = self.projects.apply(lambda x: ",".join(x.astype(str)), axis=1)
        # join the rows with a newline
        str_projects = "\n".join(str_projects)
        # append the string to the samplesheet
        ss_1_string += "\n\n[Projects_Data]\n"
        ss_1_string += ",".join(self.projects.columns) + "\n"
        ss_1_string += str_projects
        # return the string
        return ss_1_string

    def reverse_complement(self, nucleotide_string):
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        reverse = nucleotide_string[::-1]
        reverse_complement = "".join(complement.get(base, base) for base in reverse)
        return reverse_complement

    def rc_indexes(self):
        """
        Use reverse_complement() on indexes in self.df['index2]
        """
        self.df["index2"] = self.df["index2"].apply(self.reverse_complement)

    def write_to_file(self, file):
        """
        Write the samplesheet to a file
        returns nothing
        """
        ss_1_string = self.seq_only()
        with open(file, "wb") as f:
            f.write(ss_1_string.encode("ascii", "ignore"))
