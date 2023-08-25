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
        samples = request.files["samples.csv"]
        projects = request.files["projects.csv"]
        # Do something with the uploaded CSV file...
        samples_data = samples.stream.read().decode("utf-8")
        projects_data = projects.stream.read().decode("utf-8")
        samplesheet, flowcell = generate_genomics_sheet(samples_data, projects_data, request.form)
        response = make_response(samplesheet)
        response.headers["Content-Type"] = "text/csv"
        response.headers[
            "Content-Disposition"
        ] = f"attachment; filename=CTG_SampleSheet_{flowcell}.csv"
        # under construction
        return response




@app.route("/singlecell", methods=["GET", "POST"])
def upload_singlecell():
    if request.method == "POST":
        if request.form.get("singleindex") == "true":
            singleindex = True
            samplesheet_columns = ["Sample_ID", "index", "Sample_Project"]
        else:
            singleindex = False
            samplesheet_columns = ["Sample_ID", "Sample_Project", "index", "index2"]
        if (request.form.get("dev-project") =="true"):
            development_status = True
        else:
            development_status = False
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
            samplesheet_info, flexdata, feature_ref, singleindex, development_status
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


def generate_singlecell_sheet(csv_data, flexfile, feature_ref, singleindex, development_status):
    samplesheet = singleCellSheet(csv_data, flexfile, feature_ref, singleindex, development_status)
    samplesheet = samplesheet.dataDf
    return samplesheet


def generate_genomics_sheet(samples_data, projects_data, form):
    samplesheet = pep2samplesheet(StringIO(samples_data), StringIO(projects_data)
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

    return ss_string, form.get("flowcell")



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
