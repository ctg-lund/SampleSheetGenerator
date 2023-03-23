from flask import Flask, render_template, request, make_response, send_file
import csv
from io import StringIO
from samplesheet import illuminav2, singleCellSheet
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.wrappers import Response
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
        return response
    
    else:
        return render_template('forms.html')

@app.route('/singlecell', methods=['GET', 'POST'])
def upload_singlecell():
    if request.method == 'POST':
        csvfile = request.files['csvfile']
        # Do something with the uploaded CSV file...
        csv_data = csvfile.stream.read().decode('utf-8')
        samplesheet = generate_singlecell_sheet(csv_data, request.form)
        response = make_response(samplesheet)
        response.headers['Content-Type'] = 'text/csv'
        response.headers['Content-Disposition'] = f'attachment; filename=CTG_SampleSheet.csv'
        return response
    
    else:
        return render_template('singlecell_forms.html')
    
    
def generate_singlecell_sheet(csv_data, form):
    samplesheet = singleCellSheet(StringIO(csv_data))
    samplesheet = samplesheet.data.astype(str)
    return samplesheet

def generate_genomics_sheet(csv_data, form):
    samplesheet = illuminav2(StringIO(csv_data))
    samplesheet.set_read1cycles(form['readstructure'].split('-')[0])
    samplesheet.set_pipeline(form['pipeline'])
    samplesheet.set_lab_worker(form['labworker'])
    samplesheet.set_bnf_worker(form['bnfworker'])
    samplesheet.make_full_string()
    samplesheet = samplesheet.string
    return samplesheet