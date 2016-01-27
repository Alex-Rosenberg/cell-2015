###Code From the *Learning the Sequence Determinants of Alternative Splicing from Millions of Random Sequences*

In an effort to make my results reproducible and aid others in their own data analysis, I've included all of the code used in the paper. I hope this will be useful and feel free to conact me with suggestions/comments/questions.

-Alex<br><br>
If you want to work with the raw FASTQ files, you can download them from SRA. The links are in the [GEO page](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74070). Put them into a folder called fastq and everything should work.<br>
If you just want to work with the processed isoform counts, unzip all of the files in the data_gz/ folder and change the name of the folder to data/.<br>
Some of the intermediate results are very large, so I have not included them. This means you will have to go through each notebook sequentially and make the results yourself. If you want to go through a specific notebook, without going through all of the previous ones, please contact me at abros "AT" uw "DOT" edu and I can get you the required intermediate files.<br>

The plasmid files from the paper are available in in the plasmids/ directory.<br>
#### ipython notebooks:
[Notebook 0: Download SRA Files and Convert to Fastq](ipython.notebooks/Cell2015_N0_Download_Fastq_Files.ipynb)<br>
[Notebook 0A: Fastq to Isoform Counts (Alt. 3SS)](ipython.notebooks/Cell2015_N0A_A3SS_Fastq_to_Spliced_Reads.ipynb)<br>
[Notebook 0B: Fastq to Isoform Counts (Alt. 5SS)](ipython.notebooks/Cell2015_N0B_A5SS_Fastq_to_Spliced_Reads.ipynb)<br>
[Notebook 1: Library Statistics](ipython.notebooks/Cell2015_N1_Library_Statistics.ipynb)<br>
[Notebook 2: Splice Site Analysis](ipython.notebooks/Cell2015_N2_Splice_Site_Analysis.ipynb)<br>
[Notebook 3: Splicing Frame and Nonsense Mediated Decay Analysis](ipython.notebooks/Cell2015_N3_A5SS_Splicing_Frame_Analysis.ipynb)<br>
[Notebook 4: Estimating Motif Effects](ipython.notebooks/Cell2015_N4_Motif_Effect_Sizes.ipynb)<br>
[Notebook 5: Combinatorial Motif Effects](ipython.notebooks/Cell2015_N5_Combinatorial_Motif_Effects.ipynb)<br>
[Notebook 6: Learning Curves for Models Predicting Alt. 5SS Usage](ipython.notebooks/Cell2015_N6_A5SS_Model_Learning_Curves.ipynb)<br>
[Notebook 7: A Model of Alternative 5\prime Splicing](ipython.notebooks/Cell2015_N7_A5SS_Model.ipynb)<br>
[Notebook 8: Predictions on Alt. 5SS Events in the Genome](ipython.notebooks/Cell2015_N8_HAL_Genome_Predictions.ipynb)<br>
[Notebook 9: Training a Joint Model with Both the Alt. 5SS and Alt. 3SS Libraries](ipython.notebooks/Cell2015_N9_Training_Joint_A5SS_A3SS_Model.ipynb)<br>
[Notebook 10: Predicting the Effects of SNPs on Alternative 5\prime Splicing](ipython.notebooks/Cell2015_N10_A5SS_SNP_Prediction.ipynb)<br>
[Notebook 11: Predicting the Effects of SNPs on Exon Skipping](ipython.notebooks/Cell2015_N11_Predicting_Cassette_Exon_SNP_Effects.ipynb)<br>







