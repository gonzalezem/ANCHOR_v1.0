<h1>ANCHOR: A 16S rRNA gene amplicon pipeline for microbial analysis of multiple environmental samples</h1>

<h2>What is needed to run it?</h2>

<b><i>Environment</i></b>: ANCHOR currently runs on linux-like machines

<b><i>Dependencies:</i></b>
 - Mothur (used in assembling contigs. See: https://www.mothur.org/wiki/Installation)
 - BLAST (see: https://www.ncbi.nlm.nih.gov/books/NBK279671)
 - usearch9 (used for chimera detection. See: https://drive5.com/usearch/download.html)
 - python 2.7 should be already installed on the machine (https://docs.python-guide.org/starting/install/linux)


<b><i>Python libraries required:</i></b>
- numpy
- pandas
- matplotlib
- seaborn
- openpyxl
- biopython
 


<h2>All seems ok, what should I do now?</h2>


<b><i>Download:</i></b>
1. Download (or clone) ANCHOR_v1.0 in github (https://github.com/gonzalezem/ANCHOR/tree/master/ANCHOR_v1.0)
2. If not already within your system, create a link (or copy) of mothur main file into ANCHOR_v1.0/pipelineScripts/mothur. ANCHOR will look for a file called simply mothur within ANCHOR_v1.0/pipelineScripts/mothur/.
3. Create a link (or copy) of usearch9 main file (usearch9) into ANCHOR_v1.0/pipelineScripts/usearch9. ANCHOR will look for a file called simply usearch9 within ANCHOR_v1.0/pipelineScripts/usearch9/
4. Build (or link) database(s) BLAST index into ANCHOR_v1.0/db folder. Note that NCBI 16S microbial database index is included in ANCHOR download. The name should be: databasename_index (ex: 16SMicrobial_index, nt_index, rdp_index, silva_index) (see how to build an index: https://www.ncbi.nlm.nih.gov/books/NBK279688)
 

<b><i>Experiment files:</i></b>
ANCHORS needs a few files and folders:
 -  A folder containing Illumina reads (ex. PEread1_R1.fastq.gz, PEread1_R2.fastq.gz, etc.)
 -  A design file containing at least 2 columns: Samples and any_Condition_Name. Example:

| Samples  | myCondition |
| ------------- | ------------- |
| PEread1  |   Condition1  | 
| PEread2  |   Condition1  | 
| PEread3  |   Condition1  | 
| PEread4  |   Condition2  | 
| PEread5  |   Condition2  | 
| PEread6  |   Condition2  | 


<h2>Almost there:</h2>

Before running ANCHOR, prepare some room for it. The script <i>preparation_script.sh</i> from within the ANCHOR folder will do this. This script will check for dependencies and required files. It needs 3 arguments to be able to run:
 -  argument 1: raw read location (full path)
 -  argument 2: folder from where ANCHOR will be run (full path)
 -  argument 3: design file (full path)

Example:
```
cd mycomputer/myfolder/ANCHOR_v1.0
bash preparation_script.sh myIlluminaFiles/my_raw_reads mycomputer/myExperiment mycomputer/myfolder/myconditions.txt
```

The script can be run multiple times until there is no more error message.


<h2>Running ANCHOR:</h2>
If running the script preparation_script.sh didn't retrun an error, you're good to go. The last output lines from preparation_script.sh run will tell you what to do (basically customizing ANCHOR to your needs and running the main script).

```
bash bashMe.sh
```

<h2>Output:</h2>
<p>
When anchor is done a folder <b>Results_a_b_c_d</b> will be created (a-d values depend on user's input from metadata/pipe.ini).

A few folders are produced:
 -  Summary (some summary files from ANCHOR run)
 -  STAMP (inut for STAMP software)
 -  Phyloseq (input for Phyloseq)
 -  MicrobiomeAnalyst (input for microbiomeanalyst.ca)
 -  metagenomeSeq (input for metagenomeSeq)
 -  Excel (OTU table in excel format)
and files:
 -  OTU and anchor sequences (fasta files)
 -  OTU and anchor tables (txt files)
</p>

<h2>Test run</h2>
 Go inside <b>test_run</b> folder and run the following command (it takes around 2 minutes to run):
 
 ```
 bash run_test.sh
 ```
