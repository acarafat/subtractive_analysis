# Subtractive analysis
This tool use a robust subtractive approach to find novel therapeutic targets in bacteria. Starting with whole proteome of a pathogenic bacteria, it can find orthologs and paralogs of human and gut-microflora, filter those match out and reports protein targets that are unique to the pathogenic bacteria only. Those targets can be validated using in the down-stream wet-lab analysis. 

## Publication (Please Cite)
[Chakrabarty RP, Alam AR, Shill DK, Rahman A. Identification and qualitative characterization of new therapeutic targets in Stenotrophomonas maltophilia through in silico proteome exploration. Microbial Pathogenesis. 2020 Jun 10:104293.](https://www.sciencedirect.com/science/article/abs/pii/S0882401020306598)


## Requirements
Python libraries that must be installed includes: Biopython 1.68 <br />


## Pipeline and how-to
1. Download or prepare proteome sequence of target bacterium in fasta format.
2. Filter out sequence of less than 100 amino acid length. sa_utilis.py contains function filter_by_length for that. It will create CDHit_input.fasta file. Use following command in terminal to run this program:<br />
`python sa_utilis.py fasta_file sequence_length`
3. Submit the filtered proteome file CDHit_input.fasta in CDHit [web-server](http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/index.cgi?cmd=cd-hit) 
4. Download and extract output files. 
5. CDHit analysis has several outputs. This program work on files with *.fas.1.clstr.sorted extension. This output from CDHit report a cluster of fasta ids from fasta file submited to CDHit. This code takes one sequence id from each paralog cluster, and then extract the sequence of the selected id from fasta file submitted to CDHit. The following code mines on CDHit_input.fasta using the sequence ids from input.fas.1.clstr.sorted file and generate 1.no_paralogs.fasta file as output. <br />
`python paralog_filtering.py input.fas.1.clstr.sorted CDHit_input.fasta` 
6. Then run subtractive_analysis.py on 1_no_paralog.fasta It will generate file 2_non_gut_orthologs.fasta and 3_non_human_orthologs.fasta. Running following code on this script will run query of a set of amino acid sequences on NCBI database to find non-gut microbial and non-human orthologous sequences. The default parameters are env_nr database in blastp program with 0.0001 e-value. This is the most time consuming step which depends on how busy is server. It also creates blast_tmp_record.xml file as a temporary file. <br />
`python subtractive_analysis.py input.fasta`



## Discussion and bug reports
Bug Issues: https://github.com/acarafat/subtractive_analysis/issues <br />


## License announcement

>Copyright 2019 Arafat Rahman
>
>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
>
>The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
>
>THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A >PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN >ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
