import requests

file = './2/googleDrive/data/test/30833711_2_Su.txt'

#Upload the file
url = 'http://localhost:5000/upload'
file_inputs = {'input_file': open(file, 'rb')}
misc_inputs = {'input_file': '30833711_2_Su.txt', 'species': 'ara'}

response = requests.post(url, files=file_inputs, data=misc_inputs)
print(response.text)
