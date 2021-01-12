# DSML_Project
SWBio Data Science and Machine Learning Project - 2021

Extract and align a gene from a GenBank file to sequence data using web-based BLAST

Scientific papers conventionally include a gene name and an accession number of the gene/genome that is being investigated. The following script was written in Python and extracts a gene sequence from an annotated genome in the GenBank database using only the name and accession number. Furthermore, the script utilises the web-based version of BLAST to align this gene sequence with your sequence data. If you want to determine the presence of only one gene, the script allows you to screen your sequence data without having to run it through an annotation program. To demonstrate the functionality of the script the esmR (epidemic strain marker regulator) gene was chosen for alignment. The esmR resides on a pathogenicity island (cci) which encodes genes linked to enhanced virulence and transmissibility in Burkholderia cenocepacia strains.

The script can be adapted to extract and align any annotated gene from a GenBank record and requires the following input:

    -	A GenBank accession number
    -	Name of an annotated gene
    -	A FASTA file containing sequencing data

Additionally, the gene sequence extracted from the GenBank record, and associated data such as name, description, protein sequence etc., are stored in a FASTA file in the working directory and can be used for further analysis.

To run:

    1	Download the script to your working directory.

For sample data download “MEP3Bassembly.fasta” to your working directory.
                
    2	Choose any GenBank file and use its accession number and the Entrez module to download the file from GenBank.

The cci was first characterized in B. cenocepacia J2315 and its accession number (AM747721) has been used for the example.
                
    3	Choose a gene from your organism and insert the name into the script.

For the example, the gene esmR is used which is assigned to the variable “gene_of_interest”. 

    4	Choose your sequencing file and copy its file path.

For the example, the sample data should be stored in the working directory so copying the file path will not be necessary.
                
    5	Run the script

Information from the gene extracted from the GenBank record and the BLAST result will be displayed in the terminal window. Additionally, results from the BLAST alignment will be written to “blast_result.xml” and stored in the current working directory.
