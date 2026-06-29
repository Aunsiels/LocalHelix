import os
import time
import json
from flask import Flask, flash, request, redirect, url_for, send_from_directory, abort
from werkzeug.utils import secure_filename

from localhelix.analyzer import analyze_dna_to_json, initialize_all, create_page_from_body, get_html_page, \
    _clean_for_json, generate_llm_prompt

UPLOAD_FOLDER = '/path/to/the/uploads'
ALLOWED_EXTENSIONS = {'txt', 'csv', 'tsv'}

app = Flask(__name__)
# Use an absolute path for the upload folder to avoid ambiguity.
# This constructs a path to the 'data' directory inside the 'localhelix' directory.
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

initialize_all(False, app.config['UPLOAD_FOLDER'])


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/data/<path:filename>')
def serve_data(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

@app.route('/report/<string:json_filename>')
def show_report(json_filename):
    """Serves a specific report page."""
    report_path = os.path.join(app.config['UPLOAD_FOLDER'], json_filename)
    if not os.path.exists(report_path):
        abort(404, "Report not found.")
    
    with open(report_path, 'r', encoding='utf-8') as f:
        all_rss = json.load(f)

    json_data_url = url_for('serve_data', filename=json_filename)
    llm_prompt = generate_llm_prompt(all_rss, n=20)
    return get_html_page(all_rss=all_rss, json_data_url=json_data_url, llm_prompt=llm_prompt)

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
            
            # Generate a unique filename for the report JSON
            base_filename = os.path.splitext(filename)[0]
            timestamp = int(time.time())
            json_filename = f"{base_filename}_{timestamp}.json"

            output_json_path = os.path.join(app.config['UPLOAD_FOLDER'], json_filename)
            file.save(input_file)
            analyze_dna_to_json(input_file, output_json_path, False, app.config['UPLOAD_FOLDER'], False)
            return redirect(url_for('show_report', json_filename=json_filename))

    # List existing reports
    reports = [f for f in os.listdir(app.config['UPLOAD_FOLDER']) if f.endswith('.json')]
    report_links = ''.join([f'<li><a href="{url_for("show_report", json_filename=report)}">{report.replace(".json", "")}</a></li>' for report in sorted(reports, reverse=True)])
    
    report_list_html = f"<h3>Previous Reports:</h3><ul>{report_links}</ul>" if report_links else ""
    upload_form = """
    <center>
    <p>Note that this process can take several minutes.</p>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    </center>
    """
    content = upload_form + report_list_html
    body = """<div class="container-fluid p-5 bg-primary text-white text-center"> <h1>LocalHelix - Your Personal DNA 
        Analyzer</h1> <p>Upload your DNA file to generate the report</p> <form action="https://www.paypal.com/donate" method="post" 
        target="_top"> <input type="hidden" name="hosted_button_id" value="GR3D64Y7S7TU2" /> <input type="image" 
        src="https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif" border="0" name="submit" title="PayPal - 
        The safer, easier way to pay online!" alt="Donate with PayPal button" /> <img alt="" border="0" 
        src="https://www.paypal.com/en_FR/i/scr/pixel.gif" width="1" height="1" /> </form>

        </div>
          
        <div class="container mt-5">
        """ + content + """  
        </div>"""
    return create_page_from_body(body)
