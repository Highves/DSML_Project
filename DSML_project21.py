import sys
import subprocess

# Import modules needed to download and install Biopython using script

subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'Biopython'])

# Installs Biopython when you run the script. You can also install Biopython as a package in the environment. In Anaconda Navigator, go to "Environments" and select the appropriate environment (base or your own) and click "not installed". Scroll down to Biopython click the box and then install.

from Bio import SeqIO
from Bio import Entrez

# Imports the modules SeqIO and Entrez in order to read and download GenBank files, respectively. SeqIO can be used to input and output assorted sequence file formats.

#your_genome = Entrez.efetch(db="nucleotide", id="AM747721.1", idtype="acc", rettype="gb", retmode="text")
#genome = SeqIO.parse(your_genome, "gb")
#SeqIO.write(genome, "your_genome.gbk", "gb")
genome_gbk = "your_genome.gbk"

# Here you download the desired GenBank file, in this case an annotated genome, using a given GenBank accession number and save it in the current working directory as "your_genome.gbk".

faa_filename = "your_genome_in_ff.faa"
input_handle  = open(genome_gbk, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    for seq_feature in seq_record.features :
        try:
            if seq_feature.type=="CDS" :
                output_handle.write(">%s|%s\n" % (
                    seq_feature.qualifiers['gene'][0],
                    seq_feature.location,
                    ))
                pass
        except:
                continue
output_handle.close()
input_handle.close()

# The first loop accesses the attributes within the SeqRecord object (such as other objects).

# The SeqIO.parse function turns the "your_genome.gbk" file into a SeqRecord object so we can access the information in this record. SeqRecord object holds a sequence (as a Seq object) and identifiers (ID and name), description and optionally annotation and sub-features. 

# The second loop accesses the attributes in the "features" object.

# The "features" object is iterable and contains data for each annotation, such as gene name or location. However, if an annotation does not contain such information the iteration stops. The try, except and continue statement ensures iteration continues even if exceptions to the loop conditions are encountered.

# Running this code extracts the names of annotated genes and their locations within the genome sequence. These details are saved in the current working directory under a file named "your_genome_in_ff.faa" in FASTA format.

records = SeqIO.parse("your_genome_in_ff.faa", "fasta")
descriptions = list([r.description for r in records])

# The "your_genome_in_ff.faa" is returned as a SeqRecord object assigned to "records". FASTA files contain several key attributes. The description attribute contains the gene names and their locations. The description attribute is extracted and all the data stored in a list.

gene_of_interest = ['esmR']
match_in_list = [s for s in descriptions if any(xs in s for xs in gene_of_interest)]

# Here you enter your gene of interest which is matched in the list. The match is extracted and stored as an object and assigned to the variable match_in_list.

edit_list = ([s.replace('(+)', '') for s in match_in_list])
final_list = ([s.replace('(-)', '') for s in edit_list])
split_string = ([s.split("|", 1) for s in final_list])  
location = [i[1] for i in split_string] 
genes = [i[0] for i in split_string]

# Characters are removed and the gene and its location split into two seperate strings.

location_list = ([s.replace(']', '') for s in location])
f_list = ([s.replace('[', '') for s in location_list])

# Characters are further removed from the string containing the locations.

Gene = f_list[0]
x = Gene.split(":", 1)
location_start = int(x[0])
location_end = int(x[1])

# The start and stop location of the gene is split into two strings and converted to integers.

record = SeqIO.read("your_genome.gbk", "genbank")
sub_record = record[location_start:location_end]
SeqIO.write(sub_record, "your_gene.faa", "fasta")

# The SeqIO.read function returns "your_genome.gbk" file as a SeqRecord object, but uses a different set of rules than the SeqIO.parse function when converting data. The start and stop location of your gene of interest can now be used to extract the gene sequence, its identifiers (ID and name), description and sub-features. The gene sequence and its associated data are written to "your_gene.faa" which is stored in the current working directory.

from Bio.Seq import Seq
from pathlib import Path
from Bio.SeqRecord import SeqRecord

# The Path module allows Python to read data inside text files as a string. The module Seq is needed to convert a string to a Seq object. The SeqRecord module allows Seq object to be stored to a SeqRecord object.

with open("MEP3Bassembly.fasta") as infile, open("output.txt", 'w') as outfile:
    for line in infile:
        if line.startswith(">"): continue
        outfile.write(line.strip())

string = Path('output.txt').read_text()

# Opens your sequencing file of choice and concatenates every single contig and stores it as a text file. The text file is read and the sequence stored as a string.

record = SeqRecord(
    Seq(string),
)

SeqIO.write(record, "search_genome.faa", "fasta")

# The string containing your sequence is converted to a Seq object and stored in the SeqRecord object. The SeqRecord object is written to the file "search_genome.faa" and stored in the current working directory.

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

# The module NcbiblastnCommandline is needed to access the BLAST server using Python and run a BLAST analysis. The module NCBIXML is needed to handle files in the XML data format.

output = NcbiblastnCommandline(query="your_gene.faa", subject="search_genome.faa", outfmt=5, out="blast_result.xml")()[0]

#Your gene of interest is alligned with your concatenated sequencing data and the results stored in the current working directory as "blast_result.xml". 

result = open("blast_result.xml")
blast_result_record = NCBIXML.read(result)

for alignment in blast_result_record.alignments:
    for hsp in alignment.hsps:
        print ('                    BLAST Alignment \n',)
        print ("Query sequence details", "\n", "\n", sub_record.features[1])
        print ("Alignment results \n")
        print ("Query length:", hsp.align_length)
        print ('e value:', hsp.expect)
        print('identity', (hsp.identities/ hsp.align_length)*100, '%')
        print('identities', hsp.identities, "/", hsp.align_length)
        print ("Query sequence:\n", hsp.query)
        print (hsp.match)
        print ("Sequence data:\n", hsp.sbjct)

# Annotation information from your gene of interest, contained within the "your_gene.faa" file, is printed to the terminal window including whether the sequence is the reverse complement (-) or complementary gene sequence (+). Some of the data contained within the XML file is also printed to the terminal window. More information on how to access all the different data classes in the XML file can be found here: http://biopython.org/DIST/docs/tutorial/Tutorial.html#fig%3Ablastrecord

print("Done, time to celebrate!")

# If you're code has run succesfully this message should display in the terminal window.