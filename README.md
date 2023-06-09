# P3DB Peptide-Protein Parser

The python code to parse peptide-protein pairs for P3DB. Development occurs in the beta branch until testing confirms that it works as expected, at which point it is pushed to the main branch.

## Terminal Instructions:
  Call the script using <code>python PeptideProteinParser.py -i <i>input_file</i></code> in the terminal.<br>
  Optionally, use <code>-s <i>species</i></code> to choose which species' fasta file to use while parsing. If no species is selected, the program will only query online databases.<br>
<br>
<b>Input Specifications:</b><br>
  The chosen input file must follow certain formatting restrictions in order to be properly parsed.<br>
  The input file must contain the following input columns: <br>
  &emsp; <b>protein_id</b>, the designated protein ID.<br>
  &emsp; <b>peptide</b>, the peptide to search for.<br>
  Other columns are also supported; columns that appear in the input file will also appear in the output file.<br>
  If there is a column labeled <b>peptide_in_protein</b> in the input file, it will be translated into <b>original_in_protein</b> in the output file. The <b>peptide_in_protein</b> column is used in the output file to identify where in the protein the peptide is found.<br>
<br>
<b>Additional Features:</b><br>
  Use <code>-h</code> to open the help menu.<br>
  Use <code>-t <i>number</i></code> to determine the number of threads (default 10). This can usually be left at its default value.<br>

## API Instructions:
  Place the WebAPI.py script in the same folder as the PeptideProteinParser.py script and run it to activate the server. So long as the server is running, it can be queried at the '/upload' domain. The API follows the same conventions as the terminal usage.<br>
<br>
<b>Input Specifications:</b><br>
  The chosen input file must follow certain formatting restrictions in order to be properly parsed.<br>
  &emsp; The request must include an 'input_file' field uploaded to <b>files</b>, which contains the contents of the file to be parsed.<br>
  &emsp; The request must include an 'input_file' field uploaded to <b>data</b>, which contains the name of the file.<br>
  &emsp; Optionally, the request may include a 'species' field uploaded to <b>data</b>, which describes the species of the file.<br>
  An example query can be found in the WebTest.py file.<br>
