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

#Need to do:
#   Refactor
#   Make fasta selection work better


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
            indexes.append(title[0].split('.')[0])
            record.append(sequence)
        else:
            record.append('GLYMA.' + title[0].split('.')[1])
            indexes.append('GLYMA.' + title[0].split('.')[1])
            record.append(sequence)
            
        records.append(record)
            
    f.close()
    
    #Turn records into a dataframe
    db = pd.DataFrame(records,
                      index = indexes,
                      columns=['Locus', 'Sequence'])

    return db

#This function searches through the given database and returns a series of records
#A record is returned if it has the proper protein ID
def searchDatabase(proteinID, db):
    #Search for the id by Locus
    try:
        options = db.loc[[proteinID]]
    except:
        return []

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
        if customDatabase == None:
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
        return result
    except:
        pass
    
    try:
        #Check Entrez
        handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        result.append(record)
        return result
    except HTTPError:
        pass
    print("Protein", proteinID, "could not be found.")
    return result

#This function takes an input file path, and input database, and an output database
#Writes everything to appropriate output files
def findLocations(fname, inFile, db, customDatabase, thread_count):
    f = open(fname.removesuffix(".txt").removesuffix(".xlsx") + "_RECORD" + ".txt", "w")
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
        f = open(fname.removesuffix(".txt") + "_UNFOUND" + ".txt", "w")
        for x in unfound:
            f.write(x)
            f.write("\n")
        f.close()
    #Write the file for unparsed proteins
    if len(unparsed) != 0:
        f = open(fname.removesuffix(".txt").removesuffix("_UNFOUND") + "_UNPARSED" + ".txt", "w")
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

    #Find the protein ID
    proteinID = findProteinID(inFile, i)
        
    #Find and standardize the peptide sequence
    peptideSeq = db['peptide'].iloc[i]
    for character in peptideSeq:
        if not character.isupper():
            peptideSeq = peptideSeq.replace(character, '')

    #If it was previously unfound, skip it again
    if proteinID in unfound or tuple([proteinID, peptideSeq]) in found:
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
        for x in record:
            if searchProtein(f, db, proteinID, x, peptideSeq, i) == True:
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
        else:
            print("Unrecognized file format. Try a '.txt' or '.xlsx' file.")
            return 0;
    except FileNotFoundError:
        print(file, "could not be found.")
        return 0;
    except ValueError:
        print(file, "could not be parsed.");
        return 0;
    #Ensure the input file is properly formatted
    #inFile = fixColumns(inFile)

    #Parse the input file
    db = pd.DataFrame(inFile, columns=['protein_id',
                                       'verified_id',
                                       'version',
                                       'peptide',
                                       'detected_peptide',
                                       'original_in_protein',
                                       'peptide_in_protein',
                                       'verified',
                                       'MI',
                                       'D',
                                       'MZ',
                                       'C',
                                       'Protein Description'])
    db = findLocations(file, inFile, db, customDatabase, thread_count)
    print("\tFinished parsing " + file + ".")
    #Remove needless information
    remaining_lines = clean(db)
    #Send to file
    file = file.removesuffix(".txt").removesuffix(".xlsx") + "_PARSED" + ".xlsx"
    db.to_excel(file, index=False)
    #Return successful lines
    return remaining_lines

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
        print("This table has no valid pairs.")
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
#This function only exists because of Jatropha Curcas.
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

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "i:f:t:", ["input=", "fasta=", "threads="])
    except getopt.GetoptError:
        print("Unrecognized parameter.")
        return
    
    input_file = None
    customDatabase = None
    fasta = None
    thread_count = 10
    
    #Apply options
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
        elif opt in ("-t", "--threads"):
            try:
                thread_count = int(arg)
            except ValueError:
                print("Thread count must take an int. Using default value.")

    #Ensure an input file is specified
    if input_file == None:
        print("Input file must be specified. Use '-i' to choose an input file.")
        return
    if thread_count < 1:
        print("Must have at least 1 thread.")
        return

    #Create the specified database
    if(fasta != None):
        customDatabase = createDatabase(fasta)
    
    #Parse the input file
    print("Parsing", input_file)
    remaining_lines = parseFile(input_file, customDatabase, thread_count)
    print("Successfully parsed", remaining_lines, "lines.")

if __name__ == "__main__":
   main(sys.argv[1:])
