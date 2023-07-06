from flask import Flask, request, jsonify, send_file
import PeptideProteinParser
import requests

app = Flask(__name__)

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
    return send_file(output_filepath, as_attachment=True)

@app.route("/upload", methods=['POST'])
def processTwo():
    file_data = request.files
    misc_data = request.form
    
    if not file_data:
        return jsonify({'error': 'File data was not uploaded'})
    elif not misc_data:
        return jsonify({'error': 'File metadata was not uploaded'})
    elif 'input_file' not in file_data:
        return jsonify({'error': "File data does not include 'input_file' field"})
    elif 'input_file' not in misc_data:
        return jsonify({'error': "File metadata does not include 'input_file' field"})

    file_name = misc_data.get("input_file")
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
