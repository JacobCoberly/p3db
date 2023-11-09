from Bio import Entrez, SeqIO, SeqRecord, Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import os
from ast import literal_eval
from urllib.error import HTTPError
import urllib.request
from math import isnan
import sys
import getopt
import concurrent.futures
import threading

fasta_directory = "./2/googleDrive/"

#This function is used to turn the input database fasta files into an actual database.
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
        if title[0].split('.')[0] != 'GLYMA':
            record.append(title[0].split('.')[0])
            #indexes.append(title[0].split('.')[0])
        else:
            record.append('GLYMA.' + title[0].split('.')[1])
            #indexes.append('GLYMA.' + title[0].split('.')[1])
        record.append(sequence)
        record.append(file.split('/')[-1])
            
        records.append(record)
            
    f.close()
    
    #Turn records into a dataframe
    db = pd.DataFrame(records,
                      #index = indexes,
                      columns=['Locus', 'Sequence', 'Origin'])

    return db

def combineFiles(path):
    db = pd.DataFrame([], columns=['Locus', 'Sequence', 'Origin'])
    for fname in os.listdir(path):
        if(fname[-4:] == ".txt" or fname[-3:] == ".fa" or fname[-6:] == ".fasta"):
            db = pd.concat([db, createDatabase(path+fname)], ignore_index=True)
    return db
        

#This function searches through the given database and returns a series of records
#A record is returned if it has the proper protein ID
def searchDatabase(proteinID, db):
    #Search for the id by Locus
    #try:
        #options = db.loc[[proteinID]]
        #print(proteinID, len(options), options, "\n")
    #except Exception as e:
        #print(e)
        #return []

    options = db[db['Locus'] == proteinID]

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
        possibleMatches.append(option.loc['Origin'])

    return possibleMatches

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

#Helper function for findProteinID
def isValid(char):
    return char.isdigit()

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
#It returns a list of records that match (could be empty)
def searchDatabases(proteinID, customDatabase):
    Entrez.email = "cx9p9@umsystem.edu"
    result = []
    
    #Ensure GLYMA ids are properly formatted
    #GLYMA ids aren't stored in either online database
    if (proteinID[0:5]).upper() == "GLYMA":
        if type(customDatabase) == None:
            return result
        
        glymaID = 'GLYMA.'
        proteinList = proteinID.split('.')
        if proteinID[5].isdigit():
            glymaID = glymaID + proteinList[0][5:]
        elif proteinID[5] != '.':
            glymaID = glymaID + proteinList[0][6:]
        else:
            glymaID = glymaID + proteinList[1]                
        proteinID = glymaID
        result = searchDatabase(glymaID, customDatabase)
        if len(result) > 0:
            return result
        return []
    
    #Check databases
    if type(customDatabase) != None:
        try:
            #Check custom database first
            proteinID = proteinID.split('.')[0]
            result = searchDatabase(proteinID, customDatabase)
            if len(result) > 0:
                return result
        except:
            pass
    
    try:
        #Check Uniprot
        handle = urllib.request.urlopen("http://www.uniprot.org/uniprot/" + proteinID + ".xml")
        record = SeqIO.read(handle, "uniprot-xml")
        handle.close()
        result.append(record)
        result.append("Uniprot")
        return result
    except:
        pass
    
    try:
        #Check Entrez
        handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        result.append(record)
        result.append("Entrez")
        return result
    except HTTPError:
        pass
    print("Protein", proteinID, "could not be found.")
    return result

#This function takes an input file path, and input database, and an output database
#Writes everything to appropriate output files
def findLocations(fname, inFile, db, customDatabase, thread_count):
    f = open(fname.split(".")[0] + "_RECORD" + ".txt", "w")
    unfound = set()
    unparsed = set()
    found = set()
    db = standardizeInput(db, inFile)
    #Iterate through every peptide-protein pair
    lines = ((inFile, db, i, found, unfound, unparsed, f, customDatabase) for i in range(len(inFile)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=thread_count) as executor:
        executor.map(parseLine, lines)
           
    f.close()

    #Write the file for unfound proteins
    if len(unfound) != 0:
        f = open(fname.split(".")[0] + "_UNFOUND" + ".txt", "w")
        for x in unfound:
            f.write(x)
            f.write("\n")
        f.close()
    #Write the file for unparsed proteins
    if len(unparsed) != 0:
        f = open(fname.split(".")[0].removesuffix("_UNFOUND") + "_UNPARSED" + ".txt", "w")
        for x in unparsed:
            f.write(x)
            f.write("\n")
        f.close()
    print(len(unfound)+len(unparsed), "proteins failed.")
    return db

lock = threading.Lock()

#Parses a single line in a database.
#Used for threading in findLocations (above).
def parseLine(args):
    inFile = args[0]
    db = args[1]
    i = args[2]
    found = args[3]
    unfound = args[4]
    unparsed = args[5]
    f = args[6]
    customDatabase = args[7]
    #print("Beginning to parse line")

    #Find the protein ID
    proteinID = findProteinID(inFile, i)
    #print(proteinID)
        
    #Find and standardize the peptide sequence
    peptideSeq = db['peptide'].iloc[i]
    if peptideSeq == None or peptideSeq == "" or peptideSeq == "nan":
        return
    #peptideSeq = peptideSeq.upper()
    for character in peptideSeq:
        if not character.isupper():
            peptideSeq = peptideSeq.replace(character, '')
    #print(proteinID, peptideSeq)

    #If it was previously unfound, skip it again
    if proteinID in unfound or tuple([proteinID, peptideSeq]) in found:
        #print("Duplicate")
        return

    #Find the protein in the database
    record = searchDatabases(proteinID, customDatabase)

    lock.acquire()
    if len(record) == 0:
        #The protein could not be found in any database
        unfound.add(proteinID)
        lock.release()
        return
    else:
        #The protein was found at least once
        #Must check every item in the list until one matches
        flag = False
        for x in range(len(record)):
            if x % 2 == 0 and \
               searchProtein(f, db, proteinID, record[x], record[x+1], peptideSeq, i) == True:
                flag = True
        #There was no match
        if flag == False:
            unparsed.add(proteinID)
            print(peptideSeq, "could not be found in", proteinID)
        #There was a match
        else:
            found.add(tuple([proteinID, peptideSeq]))
    lock.release()

#This function updates the peptide locations for a certain protein to the
#appropriate values
def searchProtein(f, db, proteinID, record, source, peptideSeq, line):
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
        db.loc[line, 'reference_full_sequence'] = recordSeq
        db.loc[line, 'reference_source'] = source
        #Add the peptide record to the output fasta file
        f.write(">" + record.id + "|" + str(proteinID) + "\n")
        f.write(str(recordSeq) + "\n")
        return True
    else:
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

#This function runs the program on a specific file
def parseFile(file, customDatabase, thread_count):
    inFile = None
    #Attempt to read the given input file
    try:
        if(file[-4:] == ".txt"):
            inFile = pd.read_table(file)
        elif(file[-5:] == ".xlsx" or file[-4:] == ".xls"):
            inFile = pd.read_excel(file)
        elif(file[-4:] == ".csv"):
            inFile = pd.read_csv(file)
        else:
            print("Unrecognized file format.")
            return 0
    except FileNotFoundError:
        print(file, "could not be found.")
        return 0
    except ValueError:
        print(file, "could not be parsed.");
        return 0

    #Parse the input file
    db = inFile
    if 'protein_id' not in db:
        print("'protein_id' column could not be found.")
        return 0
    elif  'peptide' not in db:
        print("'peptide' column could not be found.")
        return 0
    
    db = findLocations(file, inFile, db, customDatabase, thread_count)
    print("\tFinished parsing " + file + ".")
    #Remove needless information
    remaining_lines = clean(db)
    #Send to file
    file = file.split(".")[0] + "_PARSED.csv"
    #Move the output location to the location of the script
    file = file.split("/")[-1]
    db.to_csv(file, index=False)
    #Return successful lines
    return remaining_lines

#Takes a database and removes empty columns and incomplete rows
def clean(db):
    flag = False
    #Remove excess columns
    for column in db:
        if((db[column].isnull().values.all() == True or \
           db[column].isna().values.all() == True or \
           #(column == 'verified' and db[column].values.all() == False) or \
           ("None" in set(db[column]))) and \
           column != 'verified'):
            db.drop(column, axis=1, inplace=True)
    if 'verified_id' not in db:
        #There were no verified ids
        print("This table has no valid pairs.")
        #Clear the file
        db.drop(db.index, inplace=True)
        return 0;
    if 'verified' not in db:
        db['verified'] = None

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
    return 0

#Simple standardization function that takes the fields in a database and
#strips/formats them
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
            db.loc[i, 'original_in_protein'] = 'None'
        #Strip the protein id
        proteinID = findProteinID(inFile, i)
        if(inFile.loc[i, 'protein_id'] != proteinID):
            db.loc[i, 'protein_id'] = proteinID
        
    return db

def helpMenu():
    print("Instructions:\n")
    print("  --Input (abbreviated -i):")
    print("    Required parameter that specifies the path to the file to be parsed.\n")
    print("  --Species (abbreviated -s):")
    print("    Optional parameter that specifies the species file to compare against.")
    print("    Recognized species:")
    print("      Arabidopsis (ara)")
    print("      Gossypium (gos)")
    print("      Glycine (gmax)")
    print("      Brachypodium (bra)")
    print("      Citrus (cit)")
    print("      Hordeum (hor)")
    print("      Jatropha (jat)")
    print("      Cicer (cic)\n")

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:s:t:h", ["Input=", "Species=", "Threads=", "Help"])
    except getopt.GetoptError:
        print("Unrecognized parameter.")
        return
    
    input_file = None
    customDatabase = None
    fasta = None
    path = None
    thread_count = 10
    
    #Apply options
    for opt, arg in opts:
        if opt in ("-h", "--Help"):
            helpMenu()
            return
        elif opt in ("-i", "--Input"):
            input_file = arg
        elif opt in ("-s", "--Species"):
            fasta = arg.lower()
            #Turn the species into an actual fasta file
            if fasta == 'arabidopsis' or fasta == 'ara':
                #fasta = fasta_directory + "Athaliana_447_Araport11.protein.fa"
                fasta = fasta_directory + "Arabidopsis.fa"
                path = "./2/googleDrive/arabidopsis/"
            elif fasta == 'gossypium' or fasta == 'gos':
                fasta = fasta_directory + "Gossypium.fasta"
                path = "./2/googleDrive/gossypium/"
            elif fasta == 'glycine' or fasta == 'gmax':
                fasta = fasta_directory + "Gmax_508_Wm82.a4.v1.protein.fa"
                path = "./2/googleDrive/gmax/"
            #Some species only use the online databases
            elif fasta == 'brachypodium' or fasta == 'bra' or \
                 fasta == 'citrus' or fasta == 'cit' or \
                 fasta == 'hordeum' or fasta == 'hor' or \
                 fasta == 'jatropha' or fasta == 'jat' or \
                 fasta == 'cicer' or fasta == 'cic':
                fasta = None
            else:
                print("Unrecognized species.")
                return
        elif opt in ("-t", "--Threads"):
            try:
                thread_count = int(arg)
            except ValueError:
                print("Thread count must take an int.")
                return

    #Ensure an input file is specified
    if input_file == None:
        print("Input file must be specified. Use '-i' to choose an input file.")
        return
    if thread_count < 1:
        print("Must have at least 1 thread.")
        return

    #Create the specified database
    if(fasta != None):
        #customDatabase = createDatabase(fasta)
        customDatabase = combineFiles(path)
    
    #Parse the input file
    print("Parsing", input_file)
    remaining_lines = parseFile(input_file, customDatabase, thread_count)
    print("Successfully parsed", remaining_lines, "lines.")

if __name__ == "__main__":
   main(sys.argv[1:])
