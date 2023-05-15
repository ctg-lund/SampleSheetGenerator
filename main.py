from flask import Flask, render_template, request, make_response, send_file
import csv
from io import StringIO
from samplesheet import illuminav2, singleCellSheet
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.wrappers import Response
import pandas as pd
from pprint import pprint
app = Flask(__name__)


app.wsgi_app = DispatcherMiddleware(
    Response('Not Found', status=404),
    {'/samplesheet': app.wsgi_app}
)  

@app.errorhandler(Exception)
def handle_error(e):
    return render_template("error.html", e=e), 500

@app.route('/', methods=['GET', 'POST'])
def upload():
    if request.method == 'POST':
        csvfile = request.files['csvfile']
        # Do something with the uploaded CSV file...
        csv_data = csvfile.stream.read().decode('utf-8')
        samplesheet = generate_genomics_sheet(csv_data, request.form)
        response = make_response(samplesheet)
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Disposition'] = f'attachment; filename=CTG_SampleSheet.csv'
        # under construction
        return response
    
    else:
        df = pd.read_csv('data/index_table.csv')
        index_kits = df['Index_Adapters'].unique()
        return render_template('forms.html', index_kits=index_kits)
    

@app.route('/singlecell', methods=['GET', 'POST'])
def upload_singlecell():
    if request.method == 'POST':

        # Process the sample info configuration
        samplesheet_info = combine_filestreams(request.files.getlist('samplesheets'),
                            ['Sample_ID', 'Sample_Project', 'index', 'index2']
                            )
        if 'pipeline' not in samplesheet_info.columns:
            samplesheet_info['pipeline'] = 'seqonly'
        else:
            samplesheet_info['pipeline'] = samplesheet_info['pipeline'].fillna('seqonly')
        samplesheet_info = samplesheet_info.fillna('n')

        # Process the flex configuration
        if '' == request.files['flexfile'].filename:
            flexdata = None
        else:
            flexdata = combine_filestreams(request.files.getlist('flexfile'),
                            ['sample_id','probe_barcode_ids','Sample_Project']
                            )
        # Process the feature reference
        if '' == request.files['feature_ref'].filename:
            feature_ref = None
        else:
            feature_ref = combine_filestreams(request.files.getlist('feature_ref'),
                            ['id','name','read','pattern','sequence','feature_type','Sample_Project']
                            )
            
        samplesheet = generate_singlecell_sheet(samplesheet_info.to_csv(), flexdata, feature_ref)
        response = make_response(samplesheet)
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Disposition'] = f'attachment; filename=CTG_SampleSheet.csv'
        return response
    
    else:
        return render_template('singlecell_forms.html')
    
@app.route('/lab-report', methods=['GET', 'POST'])
def upload_lab_report():
    if request.method == 'POST':
        # the uploaded file is a single pdf
        lab_report = request.files['lab_report']
        # Do something with the uploaded PDF file...

        return response
    
    else:
        return render_template('lab_report.html')
    

def generate_singlecell_sheet(csv_data, flexfile, feature_ref):
    samplesheet = singleCellSheet(StringIO(csv_data), flexfile, feature_ref)
    samplesheet = samplesheet.data
    return samplesheet

def generate_genomics_sheet(csv_data, form):
    samplesheet = illuminav2(StringIO(csv_data))
    samplesheet.sequencer = form['sequencer']
    samplesheet.read1cycles = form['readstructure'].split('-')[0]
    samplesheet.pipeline = form['pipeline']
    samplesheet.lab_worker = form['labworker']
    samplesheet.bnf_worker = form['bnfworker']
    if form['sequencer'] == 'Novaseq':
        RC = True
    else:
        RC = False
    
    
    samplesheet.get_indexes(index_kit=form['indexkit'], RC=RC )
    samplesheet.make_full_string()
    samplesheet = samplesheet.string
    return samplesheet

def combine_filestreams(filestreams, allowed_columns):
    file_list = list()
    for file in filestreams:
        file_csv = file.stream.read().decode('utf-8')
        file_csv = pd.read_csv(StringIO(file_csv))
        # Check if the file has the required columns
        if not set(allowed_columns).issubset(file_csv.columns):
            raise Exception(f'File {file.filename} does not have all the required columns {allowed_columns}')
        file_list.append(file_csv)
    filestreams = pd.concat(file_list)
    return filestreams

if __name__ == '__main__':
    # localhost:5000
    app.run(debug=True)