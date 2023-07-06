from flask import Flask, request, jsonify, send_file, render_template
import PeptideProteinParser
import requests
import csv

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route("/parse", methods=['POST'])
def process():
    json_data = request.json
    
    if not json_data:
        return jsonify({'error': 'No parameters specified'})
    if 'input_file' not in request.json:
        return jsonify({'error': "'input_file' field not specified"})

    input_file = json_data.get("input_file")
    species = None
    species = json_data.get("species")

    #determine the name of the output file
    output_filepath = input_file.split(".")[0] + "_PARSED.csv"

    #Parse the file
    PeptideProteinParser.main(["-i", input_file, "-s", species])
    #return send_file(output_filepath, as_attachment=True)

    #Return limited number of lines
    line_count = 5
    rows = []
    with open(output_filepath, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        headers = next(csvreader)
        i = 0
        for row in csvreader:
            if i < line_count:
                rows.append(row)
            i = i+1
            
    # Generate the HTML table
    table_html = '<table>'
    
    # Create the table header
    table_html += '<tr>'
    for header in headers:
        table_html += f'<th>{header}</th>'
    table_html += '</tr>'
    
    # Add the table rows
    for row in rows:
        table_html += '<tr>'
        for value in row:
            table_html += f'<td>{value}</td>'
        table_html += '</tr>'
    
    table_html += '</table>'
    
    # Return the HTML response
    return table_html

@app.route("/upload", methods=['POST'])
def processTwo():
    file_data = request.files
    misc_data = request.form
    
    if not file_data:
        return jsonify({'error': 'File data was not uploaded'})
    elif 'input_file' not in file_data:
        return jsonify({'error': "File data does not include 'input_file' field"})

    file_name = file_data.get("input_file").filename
    species = None
    species = misc_data.get("species")

    #Store the FileStorage
    file_data.get("input_file").save(file_name)

    #Parse the file
    url = 'http://localhost:5000/parse'
    json_inputs = {'input_file': file_name, 'species': species}

    response = requests.post(url, json=json_inputs)
    return response.text

if __name__ == '__main__':
    app.debug = True
    app.run()
