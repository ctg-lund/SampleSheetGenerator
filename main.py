from flask import Flask, render_template, request, make_response, send_file
import csv
from io import StringIO
from samplesheet import illuminav2, singleCellSheet
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.wrappers import Response
import pandas as pd
from pprint import pprint
app = Flask(__name__)
#Test comment


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
        uploaded_files = request.files.getlist('csvfile')
        samplesheets = list()
        # Do something with the uploaded CSV file...
        for file in uploaded_files:
            csv_data = file.stream.read().decode('utf-8')
            samplesheets.append(pd.read_csv(StringIO(csv_data)))
        csv_data = pd.concat(samplesheets)
        if 'pipeline' not in csv_data.columns:
            csv_data['pipeline'] = 'seqonly'
        else:
            csv_data['pipeline'] = csv_data['pipeline'].fillna('seqonly')
        csv_data = csv_data.fillna('n')
        samplesheet = generate_singlecell_sheet(csv_data.to_csv(), request.form)
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
    

def generate_singlecell_sheet(csv_data, form):
    samplesheet = singleCellSheet(StringIO(csv_data))
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

if __name__ == '__main__':
    # localhost:5000
    app.run(debug=True)