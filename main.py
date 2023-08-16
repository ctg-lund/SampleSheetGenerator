from flask import Flask, render_template, request, make_response, send_file
from io import StringIO
from samplesheet import singleCellSheet, pep2samplesheet
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.wrappers import Response
import pandas as pd
import sys
import traceback
from pprint import pprint

app = Flask(__name__)


app.wsgi_app = DispatcherMiddleware(
    Response("Not Found", status=404), {"/samplesheet": app.wsgi_app}
)


@app.errorhandler(Exception)
def handle_error(e):
    etype, value, tb = sys.exc_info()
    print(traceback.print_exception(etype, value, tb))
    return render_template("error.html", e=e), 500


@app.route("/", methods=["GET", "POST"])
def upload():
    if request.method == "POST":
        csvfile = request.files["csvfile"]
        # Do something with the uploaded CSV file...
        csv_data = csvfile.stream.read().decode("utf-8")
        samplesheet = generate_genomics_sheet(csv_data, request.form)
        response = make_response(samplesheet)
        response.headers["Content-Type"] = "text/csv"
        response.headers[
            "Content-Disposition"
        ] = f"attachment; filename=CTG_SampleSheet.csv"
        # under construction
        return response

    else:
        df = pd.read_csv("data/index_table.csv")
        index_kits = df["Index_Adapters"].unique()
        return render_template("forms.html", index_kits=index_kits)


@app.route("/singlecell", methods=["GET", "POST"])
def upload_singlecell():
    if request.method == "POST":
        if request.form.get("singleindex") == "true":
            singleindex = True
            samplesheet_columns = ["Sample_ID", "index", "Sample_Project"]
        else:
            singleindex = False
            samplesheet_columns = ["Sample_ID", "Sample_Project", "index", "index2"]
        # Process the sample info configuration
        samplesheet_info = combine_filestreams(
            request.files.getlist("samplesheets"), samplesheet_columns
        )
        if "pipeline" not in samplesheet_info.columns:
            samplesheet_info["pipeline"] = "seqonly"
        else:
            samplesheet_info["pipeline"] = samplesheet_info["pipeline"].fillna(
                "seqonly"
            )
        samplesheet_info = samplesheet_info.fillna("n")

        # Process the flex configuration
        if "" == request.files["flexfile"].filename:
            flexdata = None
        else:
            flexdata = combine_filestreams(
                request.files.getlist("flexfile"),
                ["sample_id", "probe_barcode_ids", "Sample_Source"],
            )
        # Process the feature reference
        if "" == request.files["feature_ref"].filename:
            feature_ref = None
        else:
            feature_ref = combine_filestreams(
                request.files.getlist("feature_ref"),
                [
                    "id",
                    "name",
                    "read",
                    "pattern",
                    "sequence",
                    "feature_type",
                    "Sample_Source",
                ],
            )

        samplesheet = generate_singlecell_sheet(
            samplesheet_info, flexdata, feature_ref, singleindex
        )
        response = make_response(samplesheet)
        response.headers["Content-Type"] = "text/csv"
        response.headers[
            "Content-Disposition"
        ] = f"attachment; filename=CTG_SampleSheet.csv"
        return response

    else:
        return render_template("singlecell_forms.html")


@app.route("/lab-report", methods=["GET", "POST"])
def upload_lab_report():
    if request.method == "POST":
        # the uploaded file is a single pdf
        lab_report = request.files["lab_report"]
        # Do something with the uploaded PDF file...
        response = "Hippity hoppity, mattis need to finish this property"
        # under construction
        return response

    else:
        return render_template("lab_report.html")


def generate_singlecell_sheet(csv_data, flexfile, feature_ref, singleindex):
    samplesheet = singleCellSheet(csv_data, flexfile, feature_ref, singleindex)
    samplesheet = samplesheet.dataDf
    return samplesheet


def generate_genomics_sheet(csv_data, form):
    samplesheet = pep2samplesheet(StringIO(csv_data))
    # set params
    samplesheet.sequencer = form.get("sequencer")
    # dev project
    if form.get("checkbox_dev"):
        samplesheet.dev_project = "Yes"
    # flowcell serial number
    samplesheet.flowcell = form.get("flowcell")
    # make sure it is filled in
    if not samplesheet.flowcell:
        raise Exception("Flowcell serial number is required!")
    if form.get('checkbox_fastq'):
        samplesheet.fastq = 'Yes'
    if form.get('checkbox_bcl'):
        samplesheet.bcl = 'Yes'
    if form.get('checkbox_bam'):
        samplesheet.bam = 'Yes'
    if form.get('checkbox_vcf'):
        samplesheet.vcf = 'Yes'
        # disable this exception for now
        # we don't know if we will support VCF in the future
        #raise Exception("VCF is not supported yet")

    # rna counts
    if form.get('checkbox_rnacounts'):
        samplesheet.rnacounts = 'Yes'
    if form.get('checkbox_fastqc'):
        samplesheet.fastqc = 'Yes'
    if form.get('checkbox_fastscreen'):
        samplesheet.fastscreen = 'Yes'
    # RC
    if form.get('checkbox_rc'):
        samplesheet.rc_indexes() 
    # generate samplesheet
    ss_string : str = ''
    if form.get("checkbox_seqonly"):
        samplesheet.seqonly_project = 'Yes'
        ss_string = samplesheet.make_ss()
    else:
        ss_string = samplesheet.make_ss()   

    return ss_string



def combine_filestreams(filestreams, allowed_columns):
    file_list = list()
    for file in filestreams:
        file_csv = file.stream.read().decode("utf-8")
        file_csv = pd.read_csv(StringIO(file_csv))
        # Check if the file has the required columns
        if not set(allowed_columns).issubset(file_csv.columns):
            raise Exception(
                f"File {file.filename} does not have all the required columns {allowed_columns}"
            )
        file_list.append(file_csv)
    filestreams = pd.concat(file_list)
    return filestreams


if __name__ == "__main__":
    # localhost:5000
    app.run(debug=True)
