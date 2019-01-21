Usage : 

perl get_fasta_stats.pl -g -T Trinity.fasta > Assembly_Stats_.txt

-----------------------------------
Longest Contig from output file :

awk 'b<$2{b=$2} END{print b}' Assembly_Stats.txt


--------------------------

get_fasta_stats - Get statistics for dna sequences in Fasta format.
The contig names optionally may be shortened by removing everything
before the word "Contig" (-s).  Summary statistics (totals) may be
displayed (-t or -T).  A fasta quality file can also be read to give
error and quality statistics (-q).  A minimum Phred/Phrap score
(-p phred_score) can be specified if quality scores are read.  Base
pairs with a Phred/Phrap quality value of 98 or 99 are counted as
edited bases, but are no used to compute other quality statistics.

If quality scores are not used, only the contig name, length, GC%, and
the number of bases that are not A, C, G, or T are listed for each
contig.  The GC% uses only A, C, G, and T bases in the calculation.
If totals are requested (-t), the number of contigs, total length of
all contigs, lengths of the shortest and longest contigs, average
contig length, average GC%, and total number of non-ACGT bases are
also output.  Column headings are displayed for the contigs, if totals
are requested.

If quality scores are used (-q), the Consed Errors/10Kb, the number
of edited bases (those bases with a quality score of 98 or 99), the
number of Phred/Phrap(#) bases, and the minimum and maximum
Phred/Phrap scores are listed for each contig along with the contig
name and length.  Phred/Phrap(#) bases is the count of the number of
bases with a Phred/Phrap base quality of (#) or greater.  This
Phred/Phrap score (#) defaults to 20, but a different value can be
specified as '-p new#'.

If both quality scores are used (-q) and totals are requested (-t),
the Consed Errors/10Kb, the number of Phred/Phrap(#) bases, and the
minimum and maximum Phred/Phrap scores are output for the entire file
along with the total contig count and size information listed above.

If extended statistics are requested (-d or -x), then mono-, di-, and
tri-nucleotide statistics are computed for each contig, and for the
entire input file if (-t) is also specified.

If a sequence length histogram is requested (-g), then it will follow
the rest of the output.


USAGE: get_fasta_stats [-a] [-d] [-e] [-g] [-G num_bars] [-n]
          [-q [-p phred_score]] [-r] [-s] [-t] [-T] [-x] [-4]
	  [fasta_input_file]
               or
       get_fasta_stats -h           <== What you are reading

  where 'num_bars' is the maximum number of frequency bars in the
            sequence length histogram.  Fewer bars may be printed
	    in order to produce bin boundaries with integer values.

	'phred_score' is a threshhold phred/Phrap score used to
	    indicate a minimum "good" quality score.

        'fasta_input_file' is the name of the input sequence file in
            Fasta format, which contains the contigs to be processed.
            Get_fasta_stats will also read a Fasta quality file named:
            "'fasta_input_file'.qual".  If 'fasta_input_file' is
	    omitted, standard input will be used, but a quality file
	    may not be read.


OPTIONS:

  -a  Produce contig assembly stats. -a may not be used with -g, -G,
      or -r.  -a computes similar stats to -g and -G, but contig
      assembly stats are printed instead of a bar chart.

  -d  Compute extended dna statistics on both strands.  -d may not be
      used with -q or -4.  If both -d and -x are specified, then -d is
      used.

  -e  Produce minimal stats on an empty input file instead of an error.

  -g  Produce a histogram graph with up to 50 bars of contig lengths
      frequencies for the entire file.  A histogram will not be
      produced if there are fewer then 3 contigs or if all contigs are
      of the same length.

  -G num_bars  Produce a histogram graph with up to 'num_bars' bars of
      contig length frequencies for the entire file.  A histogram will
      not be produced if there are fewer then 3 contigs or if all
      contigs are of the same length.

  -n  Include non-ACGT bases in extended statistics if -d or -x is
      specified;  otherwise this option is ignored.

  -p  Specify a threshhold phred/Phrap score.  For example, the
      default value is 20.  Get_fasta_stats will then display a
      Phred/Phrap20 score, which is a count of the number of bases
      with a Phred/Phrap score of 20 or better.  '-p 30' would specify
      that a Phred/Phrap30 score should be computed instead, as a
      count of the number of bases with a base quality score of 30 or
      better.  -p is ignored if a Fasta quality file is not read.

  -q  Compute and output quality statistics, using a Fasta quality
      file named: "'fasta_input_file'.qual".  If that file does not
      exist, and 'fasta_input_file' ends in ".fa", ".fna", or
      ".fasta", then a second try is made by replacing the final
      ".fa", ".fna", or ".fasta" with ".qual".  If the Fasta quality
      file cannot be opened or the 'fasta_input_file' is read from
      standard input, then -q is ignored.  Neither -d nor -x may be
      used with -q.

  -r  Output messages are to refer to input sequences as "reads",
      instead of "contigs".  May not be used with -a.
      
  -s  The contig names will be shortened by removing any prefix before
      the word "Contig", i.e., "gono.fasta.screen.Contig26" becomes
      "Contig26".

  -t  Output total fasta file statistics, as well as individual contig
      statistics.  Column headings for the contigs are also displayed.

  -T  Output total fasta file statistics, but not individual contig
      statistics.

  -x  Compute extended dna statistics on the forward strands.  -x may
      not be used with -q or -4.  If both -d and -x are specified,
      then -d is used.

  -4  Print number of 454 Newbler assembled reads for each contig.  If
      the contigs do not have the 454 Newbler " numReads=" comment,
      the number printed will be zero.  -4 may not be used with -d or
      -x.


EXAMPLES:

$ get_fasta_stats b121i21.fasta.screen.contigs

will read the fasta sequence file 'b121i21.fasta.screen.contigs' and
display only the full contig names, lengths, and GC percentages, and
number of bases that are not A, C, G, or T.

   b121i21.fasta.screen.Contig1	68	44.1%	0
   b121i21.fasta.screen.Contig2	3217	53.5%	3
   b121i21.fasta.screen.Contig3	12452	46.2%	0
   b121i21.fasta.screen.Contig4	29839	46.5%	0
   b121i21.fasta.screen.Contig5	65793	46.8%	0

$ get_fasta_stats -q -s 454AllContigs.fna

will read the fasta sequence file '454AllContigs.fna' and the fasta
quality file '454AllContigs.fna.qual'.  If the file
'454AllContigs.fna.qual' does not exist, then the program will try
'454AllContigs.qual' instead.  If either of the two quality file names
can be read, then the program computes contig quality statistics.

   contig00000	1461	4.27	0	1458	5	97	35.9%	0	1439
   contig00001	148	51.98	0	144	4	97	25.0%	0	115
   contig00002	285	6.96	0	283	8	97	34.0%	0	278
   ...
   contig05261	642	76.88	0	617	4	97	58.4%	0	204
   contig05262	146	133.81	0	138	4	97	45.9%	0	118
   contig05263	123	193.92	0	114	3	97	67.5%	0	46


$ get_fasta_stats -q -t -s mtgsp_001c04.fasta.screen.contigs

will read the fasta sequence file 'mtgsp_001c04.contigs' and the
fasta quality file 'mtgsp_001c04.contigs.qual', then compute
contig quality statistics, and summary (total) length and quality
statistics.  Displayed contig names will be shortened by removing any
prefix before the word "Contig".


   get_fasta_stats - Last Modified: October 20, 2010

   Contig statistics for fasta file:'mtgsp_001c04.fasta.screen.contigs'
         and for fasta quality file:'mtgsp_001c04.fasta.screen.contigs.qual'

   Contig   	Contig	Consed	Edited	Phred20	Minimum	Maximum		NonACGT	Longest
   Name    	Length	Err/10K	Bases	Bases	Quality	Quality	GC%	bases  	Phred20
   Contig1 	71	0.16	5	66	33	78	46.5%	0	71
   Contig2 	784	3277	0	488	0	78	34.6%	1	230
   Contig3 	436	5711	0	91	0	49	35.6%	0	7
   Contig4 	1527	1196	0	1204	0	90	46.4%	0	687
   Contig5 	48867	126.15	10	48159	0	90	36.1%	0	8411
   Contig6 	63359	31.79	0	62712	0	90	34.0%	0	12483

   Number of Contigs=6, Total bp=115044, Shortest=71, Longest=63359,
   Average length=19174, Average GC%=35.1%, Non-ACGT bases=1,
   Errors/10K=130.95, Edited bp=15, Total Phred/Phrap20 bp=112720,
   Minimum Contig Quality=0, Maximum Contig Quality=90,
   Longest Run of Phred/Phrap20 bp=12483

$ get_fasta_stats -q -t -s -p 30 mtgsp_001c04.fasta.screen.contigs

is the same as the previous example, except a Phred/Phrap30 score is
computed instead.

   get_fasta_stats - Last Modified: October 20, 2010

   Contig statistics for fasta file:'mtgsp_001c04.fasta.screen.contigs'
         and for fasta quality file:'mtgsp_001c04.fasta.screen.contigs.qual'

   Contig   	Contig	Consed	Edited	Phred30	Minimum	Maximum		NonACGT	Longest
   Name    	Length	Err/10K	Bases	Bases	Quality	Quality	GC%	bases	Phred30
   Contig1 	71	0.16	5	66	33	78	46.5%	0	71
   Contig2 	784	3277	0	436	0	78	34.6%	1	218
   Contig3 	436	5711	0	41	0	49	35.6%	0	5
   Contig4 	1527	1196	0	1103	0	90	46.4%	0	686
   Contig5 	48867	126.15	10	47879	0	90	36.1%	0	5842
   Contig6 	63359	31.79	0	62143	0	90	34.0%	0	12395

   Number of Contigs=6, Total bp=115044, Shortest=71, Longest=63359,
   Average length=19174, Average GC%=35.1%, Non-ACGT bases=1,
   Errors/10K=130.95, Edited bp=15, Total Phred/Phrap30 bp=111668,
   Minimum Contig Quality=0, Maximum Contig Quality=90,
   Longest Run of Phred/Phrap30 bp=12395

$ get_fasta_stats -4 -t 454AllContigs.fna

displays newbler assembled 454 contigs without qualities on the file
'454AllContigs.fna'.

   get_fasta_stats - Last Modified: October 20, 2010

   Contig statistics for fasta file:'454AllContigs.fna'

   Contig          Contig          NonACGT Number
   Name            Length  GC%     bases   Reads
   contig00001     1791    40.4%   0       157
   contig00002     224     44.6%   0       8
   contig00003     297     36.0%   0       5
   contig00004     126     43.7%   0       3
   contig00005     848     39.5%   0       62
   contig00006     664     50.2%   0       15
   contig00007     110     40.0%   0       2
   contig00008     114     56.1%   0       2
   contig00009     126     42.9%   0       2
   contig00010     218     45.9%   0       2
   contig00011     208     42.8%   0       3
   contig00012     186     49.5%   0       3
   contig00013     440     40.9%   0       6
   contig00014     250     51.6%   0       3
   contig00015     159     42.8%   0       2
   contig00016     272     41.2%   0       7
   contig00017     268     41.4%   0       6
   contig00018     205     47.8%   0       2
   contig00019     117     44.4%   0       2
   contig00020     226     37.6%   0       36
   contig00021     3294    40.1%   0       350
   contig00022     991     37.4%   0       107
   contig00023     263     51.7%   0       5
   contig00024     221     51.1%   0       4
   contig00025     629     49.6%   0       13
   contig00026     191     50.3%   0       5
   contig00027     972     42.7%   0       40
   contig00028     1057    41.4%   0       65
   contig00029     1224    39.7%   0       210

   Number of Contigs=29, Total bp=15691, Shortest=110, Longest=3294,
   Average length=541.1, Average GC%=42.2%, Non-ACGT bases=0,
   Total Newbler Assembled 454 Reads=1127

$ get_fasta_stats -x -t Contigs3

computes extended DNA mono-, di-, and tri-nucleotide statistics on the
sequence file 'Contigs3'.

   get_fasta_stats - Last Modified: October 20, 2010

   Contig statistics for fasta file:'Contigs3'

   contig=Contig1

   Non-ACGT bases = 0 = 0.00% of total

   Seq         Count Percent     Seq         Count Percent            Sum    Sum%

   a             771  30.49%     t             653  25.82%           1424  56.31%
   c             554  21.91%     g             551  21.79%           1105  43.69%

   total1       2529

   aa            304  12.03%     tt            227   8.98%            531  21.00%
   ac            153   6.05%     gt            132   5.22%            285  11.27%
   ag            102   4.03%     ct             83   3.28%            185   7.32%
   at            211   8.35%
   ca            177   7.00%     tg            157   6.21%            334  13.21%
   cc            120   4.75%     gg            117   4.63%            237   9.38%
   cg            174   6.88%
   ga            144   5.70%     tc            123   4.87%            267  10.56%
   gc            158   6.25%
   ta            146   5.78%

   total2       2528

   aaa           107   4.23%     ttt            63   2.49%            170   6.73%
   aac            71   2.81%     gtt            46   1.82%            117   4.63%
   ...
   taa            50   1.98%     tta            71   2.81%            121   4.79%
   tca            44   1.74%     tga            53   2.10%             97   3.84%

   total3       2527


   contig=Contig2

   Non-ACGT bases = 0 = 0.00% of total

   ...

$ get_fasta_stats -g -T spiro.contigs

get_fasta_stats - Last Modified: October 20, 2010

Contig statistics for fasta file:'spiro.contigs'

Number of Contigs=192, Total bp=1693521, Shortest=108, Longest=166615,
Average length=8820.4, Average GC%=26.2%, Non-ACGT bases=8421
Median=1111, Mode=647 occurs 3 times

Histogram of sequence lengths

Maximum          (Each '*' represents 2 occurrences,
K Bases Count     '+' represents 1/2 '*')
0       67      |*********************************+
9       99      |*************************************************+
18      1       |+
27      5       |**+
36      4       |**
45      3       |*+
54      4       |**
63      4       |**
81      1       |+
99      1       |+
108     0       |
117     1       |+
126     0       |
135     0       |
144     0       |
153     0       |
162     1       |+
171     1       |+



DATE LAST MODIFIED: October 20, 2010

