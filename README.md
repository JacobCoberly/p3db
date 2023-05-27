# P3DB Peptide-Protein Parser

The python code to parse peptide-protein pairs for P3DB. Development occurs in the beta branch until testing confirms that it works as expected, at which point it is pushed to the main branch.

<b>Instructions:</b><br>
  Call the script using <code>python PeptideProteinParser.py -i <i>input_file</i></code> in the terminal.<br>
  Optionally, use <code>-s <i>species</i></code> to choose which species' fasta file to use while parsing. If no species is selected, the program will only query online databases.<br>
<br>
<b>Input Specifications:</b><br>
  The chosen input file must follow certain formatting restrictions in order to be properly parsed.<br>
  The program supports the following input columns: <br>
  &emsp; <b>protein_id</b>, the designated protein ID.<br>
  &emsp; <b>peptide</b>, the peptide to search for.<br>
  &emsp; <b>detected_peptide</b>, the phosphorylation sites in the peptide.<br>
  &emsp; <b>peptide_in_protein</b>, the predicted location(s) of the peptide within the protein.<br>
  &emsp; <b>MI</b>.<br>
  &emsp; <b>D</b>.<br>
  &emsp; <b>MZ</b>.<br>
  &emsp; <b>C</b>.<br>
  &emsp; <b>Protein Description</b>.<br>
  If the input file is missing one or more of these columns, the corresponding output file will also lack them.<br>
  The input file must always contain the <b>protein_id</b> and <b>peptide</b> columns.<br>
<br>
<b>Additional Features:</b><br>
  Use <code>-t <i>number</i></code> to determine the number of threads (default 10). This can usually be left at its default value.
