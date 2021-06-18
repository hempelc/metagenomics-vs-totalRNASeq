#!/usr/bin/python3

# Script from David Ryder (CEFAS, Weymouth, England)

# Needs ete3 being installed and folder .etetoolkit being present in the home directory
# If folder .etetoolkit is not present in home directory, you need to run ONLY
# the ete3 command once first, delete the output, and then set -e to the
# .etetoolit/taxa.sqlite folder that has now been generated in your folder.

import getopt,sys,sqlite3,os

# Retrieve file names from the command line
def retrieveArguments():
  blastFileName=''
  lineageFileName=''
  outputFileName=''
  etetoolkit=''
  taxonColumn=13

  myopts, args = getopt.getopt(sys.argv[1:],"b:l:o:t:e:")

  helpMessage = "-b blast alignments -l taxonomic lineages -o output -t 13 -e etetoolkit"

  if len(myopts) != 5:
    print("Usage: {} {}".format(sys.argv[0],helpMessage))
    exit();

  for o, a in myopts:
    if o == '-l':
      lineageFileName = a
    elif o == '-b':
      blastFileName = a
    elif o == '-o':
      outputFileName = a
    elif o == '-t':
      taxonColumn = int(a)
    elif o == '-e':
      etetoolkit = a
    else:
      print("Usage: {} {}".format(sys.argv[0],helpMessage))
      exit();

  print("Lineage file: {}\nBlast alignment file: {}\nOutput file: {}\nColumn used to store taxonomic id: {}\netetoolkit/taxa.sqlite location: {}".format(lineageFileName,blastFileName,outputFileName,taxonColumn,etetoolkit))

  return(lineageFileName,blastFileName,outputFileName,taxonColumn,etetoolkit)

lineageFileName,blastFileName,outputFileName,taxonColumn,etetoolkit = retrieveArguments()

# Read file containing results
lineageFile = open(lineageFileName,'r')

# Open file for writting out results
outputFile = open(outputFileName,'w')

# Create a connection to the database
conn = sqlite3.connect(etetoolkit.format(os.getlogin()))

c = conn.cursor()
c.execute("PRAGMA foreign_keys=ON")

taxonomicLineage = dict()
taxon2SpeciesName = dict()

for line in lineageFile:
  line = line.strip('\n')
  splitLine = line.split('\t');
  taxonID = splitLine[0];
  SpeciesName = splitLine[1];
  LowestLevel = splitLine[2];
  lineage = splitLine[3];
  lineageSplit = lineage.split(',');
  taxonLineage = splitLine[4];
  taxonLineageSplit = taxonLineage.split(',');

  taxonLineageRank = []

  if taxonID not in taxon2SpeciesName:
    taxon2SpeciesName[taxonID] = [LowestLevel,SpeciesName]

  #print(line)
  for taxon in taxonLineageSplit:
    c.execute('SELECT rank FROM species WHERE taxid = ?',[taxon])
    #print("\t{}".format(taxon))
    result = c.fetchone()
    if result:
      taxonRank = result[0].strip('\t')
    else:
      taxonRank = 'No Result'
    #print("TaxonID:{}\nTaxon Rank:{}\n".format(taxon,taxonRank))

    taxonLineageRank.append(taxonRank)

  for i in range(1,len(taxonLineageSplit)):
    #outputFile.write("{}\t{}\t{}\t{}\n".format(taxonID,LowestLevel,lineageSplit[i],taxonLineageRank[i]))

    if taxonID in taxonomicLineage:
      taxonomicLineage[taxonID][taxonLineageRank[i]] = lineageSplit[i]
    else:
      taxonomicLineage[taxonID] = {taxonLineageRank[i]: lineageSplit[i]}

blastFile = open(blastFileName,'r')
taxonomicLevels = ['superkingdom','kingdom','phylum','subphylum','class','subclass','order','suborder','infraorder','family','genus']

for line in blastFile:
  line = line.strip('\n')
  splitLine = line.split('\t');
  taxonID = splitLine[taxonColumn-1].split(';')[0]

  outputFile.write(line)

  if (taxonID in taxonomicLineage):
    outputFile.write("\t{}\t{}".format(taxon2SpeciesName[taxonID][0],taxon2SpeciesName[taxonID][1]))

    for level in taxonomicLevels:
      if (level in taxonomicLineage[taxonID]):
        outputFile.write("\t{}".format(taxonomicLineage[taxonID][level]))
      else:
        outputFile.write("\tNA")
  else:
    for i in range(0,len(taxonomicLevels) + 2):
      outputFile.write("\tUnknown")
  outputFile.write("\n")
