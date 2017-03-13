## Subtractive analysis
This project has been created for drug target finding project collaborating with Rajan, JNU. 

One key paper we are follwoing is "Finding Potential Therapeutic Targets against 
Shigella flexneri through Proteome Exploration" by Mohammad Uzzal Hossain et al. from MBSTU.

### Project structure
- All scripts are kept in `/Analysis/Scripts` folder.


### Direction for doing subtractive analysis
1. Create a file containing proteome of an organism.
2. Filter out sequence of less than 100 amino acid length.
`sa_utilis.py` contains function `filter_by_length` for that.
It will create `CDHit_input.fasta` file.
3. Submit proteome file `CDHit_input.fasta` to CD hit. 
4. On CD-hit output, run `paralog_filtering.py`
It will generate file `1_no_paralog.fasta`
5. Then run `subtractive_analysis.py` on `1_no_paralog.fasta`
It will generate file `2_non_gut_orthologs.fasta` and `3_non_human_orthologs.fasta`
6. `3_non_human_orthologs.fasta` is ready for further analysis. 


