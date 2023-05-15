from Bio import Entrez, SeqIO, SeqRecord, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import os
from ast import literal_eval
from urllib.error import HTTPError
import urllib.request
from math import isnan

#Need to do:
#   Refactor (Rather important)
#   Reformat for terminal


#This function is used to turn the input database fasta files into an actual database.
#In desparate need of refactoring
def createDatabase(file):
    print("Creating database from file", file)
    f = open(file, "r")
    records = []
    indexes = []
    #Add each id to the records frame
    for title, sequence in SimpleFastaParser(f):
        record = []
        title = title.split()
        for i in range(len(title)):
            title[i] = (title[i].split('='))[-1]
            title[i] = title[i].upper()
        if len(title) == 6:
            #MemoryError
            record.append(title[3])
            #Transcript
            record.append(title[2])
            #Glyma
            record.append(title[0])
            #ID
            record.append(title[4])
            #annot-version
            record.append(title[5])
            #PACID
            record.append(title[1])
            indexes.append(title[3])
            record.append(sequence)
            
        elif len(title) == 1:
            record.append(title[0])
            indexes.append(title[0])
            record.append(sequence)
            
        elif len(title) == 4:
            record.append(title[0])
            indexes.append(title[0])
            record.append(sequence)
            
        elif len(title) > 6:
            record.append(title[0])
            record.append(title[0][0:-2])
            indexes.append(title[0][0:-2])
            record.append(sequence)
        records.append(record)
            
    f.close()
    #Turn records into a dataframe
    if len(title) == 6:
        db = pd.DataFrame(records,
            index = indexes,
            columns=['Locus', 'Transcript', 'Glyma', 'ID', 'ANNOT-Version', 'PACID', 'Sequence'])
    elif len(title) == 1 or len(title) == 4:
        db = pd.DataFrame(records,
            index = indexes,
            columns=['Locus', 'Sequence'])
    elif len(title) > 6:
        db = pd.DataFrame(records,
            index = indexes,
            columns=['ID', 'Locus', 'Sequence'])
    return db

#This function searches through the given database and returns a series of records
#A record is returned if it has the proper protein ID
def searchDatabase(proteinID, db, column):
    #Search for the id by Locus
    try:
        options = db.loc[[proteinID]]
    except:
        return None

    #option = options.iloc[0]
    possibleMatches = []
    for i in range(len(options)):
        #Create record
        option = options.iloc[i]
        sequence = option.loc['Sequence']
        #Was originally 'ID' for the GLYMA files
        ident = option.loc['Locus']
        record = SeqRecord.SeqRecord(Seq.Seq(sequence), id=ident)
        possibleMatches.append(record)

    if len(possibleMatches) == 1:
        return record
    else:
        return possibleMatches

def isValid(char):
    return char.isdigit()

#This function takes two strings, string and sub, and returns all of the positions
#in string where sub is located (where the first element is location 1).
def findSubstring(string, sub):
    stringLength = len(string)
    subLength = len(sub)
    maxPos = stringLength-subLength
    locations = []

    for i in range(maxPos+1):
        if string[i : i+subLength] == sub:
            locations.append(i+1)
    return locations

#This function returns the protein ID found in a certain row in a given database
def findProteinID(inFile, row):
    string = str(inFile['protein_id'].iloc[row]).upper().strip()
    #If there are multiple sections, take the second one
    if len(string.split('|')) > 1:
        proteinID = string.split('|')[1]
        #Remove non-numeric characters
        proteinID = list(filter(isValid, proteinID))
        proteinID = ''.join(proteinID)
    #If there are not multiple sections, take the whole thing
    #Includes non-numeric characters
    else:
        proteinID = string
    return proteinID

#This function searches through available databases for matches to the proteinID
#It returns either a record or a list of records that match
def searchDatabases(proteinID):
    Entrez.email = "cx9p9@umsystem.edu"
    #Standard-format ids should be entirely numeric
    if proteinID.isdigit():
        try:
            handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
        except HTTPError:
            print("\tHTTPError: Bad Request. Protein id",
            proteinID,
            "could not be found.")
            return None
    #GLYMA ids begin with the word 'GLYMA'
    elif (proteinID[0:5]).upper() == "GLYMA":
        glymaID = 'GLYMA.'
        proteinList = proteinID.split('.')
        if proteinID[5].isdigit():
            glymaID = glymaID + proteinList[0][5:]
        elif proteinID[5] != '.':
            glymaID = glymaID + proteinList[0][6:]
        else:
            glymaID = glymaID + proteinList[1]
        column = 'Locus'                
        return searchDatabase(glymaID, glymaDatabase, column)
    #Some standard and nonstandard ids are letter-digit combos
    else:
        try:
            #Check custom database first
            column = 'Locus'
            proteinID = proteinID.split('.')[0]
            result = searchDatabase(proteinID, customDatabase, column)
            if result is not None:
                return result
        except:
            pass
        try:
            print("Custom database failed to locate", proteinID, ". Checking Uniprot.")
            handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + proteinID + ".xml")
            record = SeqIO.read(handle, "uniprot-xml")
            handle.close()
            return record
        except:
            #print("Online database failed. Checking standard database.")
            pass
        #Attempt standard database
        try:
            print("Uniprot failed to locate", proteinID, ". Checking Standard.")
            handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
        except HTTPError:
            #print("\tHTTPError: Bad Request. Protein id",
            #proteinID,
            #"could not be found.")
            print("Protein", proteinID, "could not be found.")
            return None
    print("Protein", proteinID, "could not be found.")
    return None

#This function takes an input file path, and input database, and an output database
#Writes everything to appropriate output files
def findLocations(fname, inFile, db):
    f = open(fname.removesuffix(".txt").removesuffix(".xlsx") + "_RECORD" + ".txt", "w")
    unfound = set()
    unparsed = set()
    found = set()
    db = standardizeInput(db, inFile)
    #Iterate through every peptide-protein pair
    for i in range(len(inFile)):
        #Find the protein ID
        proteinID = findProteinID(inFile, i)
        
        #Find and standardize the peptide sequence
        peptideSeq = db['peptide'].iloc[i]
        for character in peptideSeq:
            if not character.isupper():
                peptideSeq = peptideSeq.replace(character, '')
        
        #If it was previously unfound, skip it again
        if proteinID in unfound or tuple([proteinID, peptideSeq]) in found:
            continue
        
        #Find the protein in the database
        record = searchDatabases(proteinID)
        if record is None:
            #print("The protein could not be found in any database.")
            unfound.add(proteinID)
            continue
            #break
        if type(record) != list:
            recordSeq = record.seq
            if searchProtein(f, db, proteinID, record, peptideSeq, i) == False:
                unparsed.add(proteinID)
                print(peptideSeq, "could not be found in", proteinID)
        else:
            #Must check every item in the list until one matches
            flag = False
            for x in record:
                if searchProtein(f, db, proteinID, x, peptideSeq, i) != False:
                    flag = True
            if flag == False:
                unparsed.add(proteinID)
                print(peptideSeq, "could not be found in", proteinID)
            else:
                found.add(tuple([proteinID, peptideSeq]))
           
    f.close()
    
    if len(unfound) != 0:
        f = open(fname.removesuffix(".txt") + "_UNFOUND" + ".txt", "w")
        for x in unfound:
            f.write(x)
            f.write("\n")
        f.close()
    if len(unparsed) != 0:
        f = open(fname.removesuffix(".txt").removesuffix("_UNFOUND") + "_UNPARSED" + ".txt", "w")
        for x in unparsed:
            f.write(x)
            f.write("\n")
        f.close()
    print(len(unfound)+len(unparsed), "proteins failed.")
    return db

#This function updates the peptide locations for a certain protein to the
#appropriate values
def searchProtein(f, db, proteinID, record, peptideSeq, line):
    recordSeq = record.seq
    
    #Check if there is a preexisting number to verify
    if(db['original_in_protein'].iloc[line] != '_'):
        db.loc[line, 'verified'] = verifyLocation(proteinID, peptideSeq, recordSeq, db['original_in_protein'].iloc[line])   

    if len(findSubstring(recordSeq, peptideSeq)) != 0:
        #Matched the peptide
        string = str(findSubstring(recordSeq, peptideSeq))
        string = string.removesuffix(']')
        string = string.removeprefix('[')
        db.loc[line, 'peptide_in_protein'] = string
        db.loc[line, 'verified_id'] = proteinID
        #Add the peptide record to the special file
        f.write(">" + record.id + "|" + str(proteinID) + "\n")
        f.write(str(recordSeq) + "\n")
        return True
    else:
        #The peptide could not be found, make a note
        #print(peptideSeq, "could not be found in", proteinID)
        return False

#This function verifies whether or not a given peptide location is correct
def verifyLocation(proteinID, peptideSeq, recordSeq, location):
    #Cycle through protein versions looking for relevant data
    findSub = findSubstring(recordSeq, peptideSeq)
    string = str(findSub)
    string = string.removesuffix(']')
    string = string.removeprefix('[')
    #Check if the given location matches
    if string == str(location):
        return True
    #Check if the given location is partially right
    elif location in findSub:
        return True
    #The given location is wrong
    return False

#This function runs the program on a specific file in a specific directory
def parseFile(directory, filename):
    inFile = None
    file = createFilePath(directory, filename)
    if(file[-4:] == ".txt"):
        inFile = translateFile(file)
    elif(file[-5:] == ".xlsx"):
        inFile = translateTable(file)
        inFile = fixColumns(inFile)
    
    db = pd.DataFrame(inFile, columns=['protein_id', 'verified_id', 'version', 'peptide','detected_peptide', 'original_in_protein', 'peptide_in_protein', 'verified', 'MI', 'D', 'MZ', 'C', 'Protein Description'])
    db = findLocations(file, inFile, db)
    print("\tFinished parsing " + file + "\n")
    #Remove needless information
    remaining_lines = clean(db)
    #Send to file
    file = file.removesuffix(".txt") + ".xlsx"
    db.to_excel(file, index=False)
    return remaining_lines


def main():
    if(filename != '*'):
        parseFile(directory, filename)
    else:
        #Parse the entire directory
        files = []
        for r, d, f in os.walk(os.path.join('./data/', directory)):
            files = f
        #Files identified
        count = len(files)
        y = 0
        #Remove bad files
        for x in range(count):
            if files[y][-11:] == "UNFOUND.txt" or \
               files[y][-12:] == "UNPARSED.txt" or \
               files[y][-10:] == "RECORD.txt" or \
               files[y][-5:] == ".xlsx" or \
               files[y] == ".DS_Store":
                files.pop(y)
                y = y-1
            y = y+1
        #Perform the calculations
        results = open("table_sizes.txt", 'w')
        for f in files:
            remaining_lines = parseFile(directory, f)
            results.write(f.split('_')[0])
            results.write(" ")
            results.write(str(remaining_lines))
            results.write("\n")
        results.close()
            
    print("All files parsed.")

#This function creates a file path given the directory and filename
def createFilePath(directory, filename):
    return "./data/"+directory+"/"+filename

#Takes a file path and turns the file at that location into a dataframe
def translateFile(file):
    try:
        print("\tNow parsing " + file)
        return pd.read_table(file)
    except FileNotFoundError:
        print("The file could not be found.")
        return None
    except ValueError:
        print("The file could not be parsed.")
        return None

#Takes a file path and turns the file at that location into a dataframe
def translateTable(file):
    try:
        print("\tNow parsing " + file)
        return pd.read_excel(file)
    except FileNotFoundError:
        print("The file could not be found.")
        return None
    except ValueError:
        print("The file could not be parsed.")
        return None

#Takes a database and removes empty columns and incomplete rows
def clean(db):
    flag = False
    #Remove excess columns
    for column in db:
        if(db[column].isnull().values.all() == True or \
           db[column].isna().values.all() == True or \
           (column == 'verified' and db[column].values.all() == False) or \
           ("None" in set(db[column]))):
            db.drop(column, axis=1, inplace=True)
            if column == 'verified_id':
                flag = True
    if flag:
        #There were no verified ids
        return 0;
    #Remove excess rows
    depth = 0 #actual position in the table
    index = 0 #index name in the table
    pairs = set()
    while depth < len(db):
        if pd.isna(db.iloc[depth].loc['verified_id']) or \
           tuple([db.iloc[depth].loc['protein_id'], db.iloc[depth].loc['peptide']]) in pairs:
            db.drop(index, axis=0, inplace=True)
            depth=depth-1
        else:
            pairs.add(tuple([db.iloc[depth].loc['protein_id'], db.iloc[depth].loc['peptide']]))
        depth=depth+1
        index=index+1
    return depth

#This function takes a dataframe, read in from a file, that has nonstandard
#column layout. It "fixed" the data to a more standardized format.
def fixColumns(inFile):
    db = pd.DataFrame(columns=['protein_id', \
                               'verified_id', \
                               'version', \
                               'peptide', \
                               'detected_peptide', \
                               'original_in_protein', \
                               'peptide_in_protein', \
                               'verified', \
                               'MI', \
                               'D', \
                               'MZ', \
                               'C', \
                               'Protein Description'])
    for i in range(len(inFile)):
        proteins = inFile.iloc[i].loc['Proteins'].split(';')
        #originalLocs = str(inFile.iloc[i].loc['Positions within proteins']).split(';')
        #Standardize peptide
        peptideSeq = inFile['Modified sequence'].iloc[i]
        ph = peptideSeq.replace('(ph)', '*')
        for character in peptideSeq:
            if not character.isupper():
                peptideSeq = peptideSeq.replace(character, '')
        for character in ph:
            if not character.isupper() and not character == '*':
                ph = ph.replace(character, '')

        for j in range(len(proteins)):
            proteins[j] = proteins[j].split('|')[-1]
            #originalLocs[j] = int(float(originalLocs[j]))

            row = {'protein_id':proteins[j], \
                   'peptide':peptideSeq, \
                   'detected_peptide':ph, \
                   'original_in_protein':'_', \
                   'C':inFile.iloc[i].loc['Charge']}
            db.loc[len(db.index)] = row
    
    return db

#Simple standardization function that takes the fields in a database and
#strips them
def standardizeInput(db, inFile):
    #Standardize the fields
    for i in range(len(inFile)):
        #Strip the peptide
        try:
            db.loc[i, 'peptide'] = str(inFile.loc[i, 'peptide']).strip()
        except:
            db.loc[i, 'peptide'] = '_'
        #Strip the peptide in protein
        try:
            db.loc[i, 'original_in_protein'] = str(inFile.loc[i, 'peptide_in_protein']).strip()
        except:
            print("None 3")
            db.loc[i, 'original_in_protein'] = 'None'
        #Strip the protein id
        proteinID = findProteinID(inFile, i)
        if(inFile.loc[i, 'protein_id'] != proteinID):
            db.loc[i, 'protein_id'] = proteinID
        
    return db

glymaDatabase = createDatabase("Gmax_508_Wm82.a4.v1.protein.fa")
customDatabase = createDatabase("Athaliana_447_Araport11.protein.fa")

directory = "Jatropha curcas"
filename = '30626061.xlsx'

main()
