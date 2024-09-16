# From https://flask.palletsprojects.com/en/2.3.x/patterns/fileuploads/

import os
from flask import Flask, flash, request, redirect, url_for
from werkzeug.utils import secure_filename

from personal_dna_analyzer.analyzer import main, initialize_all

UPLOAD_FOLDER = '/path/to/the/uploads'
ALLOWED_EXTENSIONS = {'txt', 'csv', 'tsv'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = "personal_dna_analyzer/data/"


initialize_all(False, app.config['UPLOAD_FOLDER'])


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            input_file = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            output_file = os.path.join(app.config['UPLOAD_FOLDER'], "temp.html")
            file.save(input_file)
            main(input_file, output_file, False, app.config['UPLOAD_FOLDER'], False)
            with open(output_file, 'r') as f:
                return f.read()
    return '''
    <!doctype html>
    <title>Upload your DNA file</title>
    <h1>Upload your DNA file to generate the report</h1>
    <p>Note that this process can take several minutes.</p>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''
