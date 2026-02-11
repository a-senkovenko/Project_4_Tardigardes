# Tardigrades in space 

Tardigrades are microscopic extremofile ecdysozoans capable to survive after desiccation, freezing,
and severe osmotic stress by entering a ametabolic state called cryptobiosis. 
Exact pathways providing tardigrade extremotolerance need to be found.

Gene prediction is one of the principal approaches to assess the properties
of undescribed genes. Generally deployed techniques include ab initio prediction and similarity-based search. 

Therefore, in this project analysis of [*Ramazzottius varieornatus*](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=947166) genome was performed with following search of
proteins associated with radiation damage negation.

Before the analyzis download [genome](http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz), [annotation](https://drive.google.com/file/d/1hCEywBlqNzTrIpQsZTVuZk1S9qKzqQAq/view?usp=sharing) and prebuilt protein [fasta](https://drive.google.com/file/d/12ShwrgLkvJIYQV2p1UlXklmxSOOxyxj4/view?usp=sharing).


## Step 1. Structural annotation

If we have built genome annotation in gff format use it to restore protein fasta file via numerous methods.
The script [getAnnoFast.pl](http://augustus.gobics.de/binaries/scripts/getAnnoFasta.pl) can be used for this task:

```bash
chmod +x getAnnoFast.pl
perl getAnnoFast.pl augustus.whole.gff
```

After the protein fasta is acquired count proteins in it:

```bash
grep -c ">" augustus.whole.aa
```
*Output: 16435*

How can we narrow it down to a list of potential candidates that we can then verify by experimentation?

The answer is to try to find DNA-associated proteins that involved in DNA-reparation. 
In this case it is recommended to check these proteins by their physical localization and collect only nuclear.
Then try to find a functional annotation for the remaining sequences.

## Step 2. Physical localization

The list of DNA-associated peptides obtained via tandem mass spectrometry.
Therefore, need to find relevant proteins from protein fasta.

Make new environment and add [blast](https://anaconda.org/channels/bioconda/packages/blast/overview) and [diamond](https://github.com/bbuchfink/diamond) there.

```bash
conda install blast diamond
```
These tools need to prepare database before the first run and then perform the search.

```bash
# blast+
makeblastdb -in augustus.whole.aa -dbtype prot  -out ./blast/blastdb
blastp -db blastdb -query peptides.fa -outfmt 6 -out ./blast/blastout.txt 

# diamond
mkdir ./diamond
diamond makedb --in augustus.whole.aa  --db ./diamond/diamonddb.dmnd
diamond blastp -d diamonddb.dmnd -q peptides.fa  -f 6 -o ./diamond/diamondout.txt --very-sensitive
```
Blast is not just protein-specific tool, therefore blastp use is needed. Blast database consists of 7 files.
Diamond is more compact.

Output by blast contains 118 entries and the needed ones should be extracted.

```bash
awk '{print $2}'./blast/blastout.txt | sort -u > ids.txt
samtools faidx augustus.whole.aa $(cat ids.txt) > hits.fa
```
There are 34 sequences from original protein list that contain these peptides.

Output of **diamond** is more precise:

```bash
8	g4106.t1	100	18	0	0	1	18	222	239	5.55e-07	40.8
20	g12510.t1	100	18	0	0	1	18	425	442	1.43e-06	39.7
21	g12510.t1	100	22	0	0	1	22	443	464	2.46e-09	47.8
29	g4106.t1	100	18	0	0	1	18	222	239	5.55e-07	40.8
```
Only two entries from original protein fasta are present in diamond output. However, for the next steps better use the list with 34 candidate sequences.

## Step 3. Localization prediction

To evaluate the obtained data it is recommended to predict the localization of these proteins within the cells.
Focus on the candidates present within nucleus.
Use file hits.fa from previous step for each search

### [WoLF PSORT](https://wolfpsort.hgc.jp/)

This web-available tool predicts localization of proteins and possible organisms they belong.

To begin search select an organism type (Animal) and input method (From file), press submit.

<details>
    <summary>Output of Wolf PSORT</summary>

```
#Warning sequence g16318.t1 does not start with M 
#Warning sequence g16368.t1 does not start with M g10513.t1 details nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1
g10514.t1 details nucl: 19, cyto_nucl: 15, cyto: 9, extr: 3, mito: 1
g11320.t1 details plas: 24.5, extr_plas: 16, extr: 6.5, lyso: 1
g11513.t1 details cyto: 17, cyto_nucl: 12.8333, cyto_mito: 9.83333, nucl: 7.5, E.R.: 3, mito: 1.5, plas: 1, pero: 1, golg: 1
g11806.t1 details nucl: 18, cyto_nucl: 11.8333, mito: 5, extr: 4, cyto: 3.5, cyto_pero: 2.66667, cysk_plas: 1
g11960.t1 details nucl: 32
g12388.t1 details extr: 
25, plas: 4, mito: 2, lyso: 1
g12510.t1 details plas: 29, cyto: 3
g12562.t1 details extr: 30, lyso: 2
g1285.t1 details extr: 25, plas: 5, mito: 1, lyso: 1
g13530.t1 details extr: 13, nucl: 6.5, lyso: 5, cyto_nucl: 4.5, plas: 3, E.R.: 3, cyto: 1.5
g14472.t1 details nucl: 28, plas: 2, cyto: 1, cysk: 1
g15153.t1 details extr: 32
g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
g16318.t1 details nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1
g16368.t1 details nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1
g2203.t1 details plas: 29, nucl: 2, golg: 1
g3428.t1 details mito: 18, cyto: 11, extr: 2, nucl: 1
g3679.t1 details extr: 26, mito: 2, lyso: 2, plas: 1, E.R.: 1
g4106.t1 details E.R.: 14.5, E.R._golg: 9.5, extr: 7, golg: 3.5, lyso: 3, pero: 2, plas: 1, mito: 1
g4970.t1 details plas: 32
g5237.t1 details plas: 24, mito: 8
g5443.t1 details extr: 28, nucl: 3, cyto: 1
g5467.t1 details extr: 27, plas: 4, mito: 1
g5502.t1 details extr: 31, lyso: 1
g5503.t1 details extr: 29, plas: 1, mito: 1, lyso: 1
g5510.t1 details plas: 23, mito: 7, E.R.: 1, golg: 1
g5616.t1 details extr: 31, mito: 1
g5641.t1 details extr: 31, lyso: 1
g5927.t1 details nucl: 30.5, cyto_nucl: 16.5, cyto: 1.5
g702.t1 details extr: 29, plas: 2, lyso: 1
g7861.t1 details nucl: 16, cyto_nucl: 14, cyto: 8, plas: 5, pero: 1, cysk: 1, golg: 1
g8100.t1 details nucl: 16.5, cyto_nucl: 12.5, cyto: 7.5, plas: 5, extr: 2, E.R.: 1
g8312.t1 details nucl: 15.5, cyto_nucl: 15.5, cyto: 12.5, mito: 2, plas: 1, golg: 1
```
</details>

```bash
grep -c "nucl" ./search_res/wolfpsort_out.txt
#17
```
Therefore, we have 17 candidates by this round of search.

### [TargetP](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)

TargetP tool predicts subcellular localization of eukaryotic proteins by their N-terminal presequences.

Use the following settings:
Organism group: Non-plant
Output format: Short output (no figures)

<details>
    <summary>TargetP output</summary>

```
# TargetP-2.0	Organism: Non-Plant	Timestamp: 20260207172436
# ID	Prediction	OTHER	SP	mTP	CS Position
g10513.t1	OTHER	0.999999	0.000001	0.000000	
g10514.t1	OTHER	0.999543	0.000349	0.000107	
g11320.t1	SP	0.000184	0.999816	0.000000	CS pos: 20-21. AYS-AG. Pr: 0.7236
g11513.t1	OTHER	0.999434	0.000401	0.000164	
g11806.t1	OTHER	0.998977	0.000887	0.000136	
g11960.t1	OTHER	0.999996	0.000002	0.000002	
g12388.t1	SP	0.000490	0.999481	0.000029	CS pos: 16-17. ASA-SS. Pr: 0.6485
g12510.t1	OTHER	0.999738	0.000099	0.000163	
g12562.t1	SP	0.000076	0.999923	0.000001	CS pos: 16-17. SYA-AN. Pr: 0.7910
g1285.t1	SP	0.003029	0.996798	0.000173	CS pos: 16-17. ASA-TS. Pr: 0.7127
g13530.t1	SP	0.116007	0.883840	0.000153	CS pos: 19-20. TIP-FT. Pr: 0.3552
g14472.t1	OTHER	0.999999	0.000001	0.000000	
g15153.t1	SP	0.000014	0.999986	0.000000	CS pos: 16-17. AYA-AN. Pr: 0.8378
g15484.t1	OTHER	0.999980	0.000010	0.000010	
g16318.t1	OTHER	0.997047	0.002953	0.000000	
g16368.t1	OTHER	0.996693	0.003307	0.000000	
g2203.t1	OTHER	0.999869	0.000031	0.000100	
g3428.t1	OTHER	0.999903	0.000033	0.000064	
g3679.t1	SP	0.001755	0.998023	0.000222	CS pos: 18-19. TFA-AR. Pr: 0.5523
g4106.t1	OTHER	0.729658	0.266917	0.003425	
g4970.t1	OTHER	0.999996	0.000003	0.000001	
g5237.t1	OTHER	0.999545	0.000345	0.000111	
g5443.t1	OTHER	0.952853	0.043784	0.003363	
g5467.t1	SP	0.000096	0.999845	0.000059	CS pos: 16-17. ASA-GS. Pr: 0.6543
g5502.t1	SP	0.001134	0.998823	0.000043	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5503.t1	SP	0.001222	0.998720	0.000058	CS pos: 16-17. ASA-GS. Pr: 0.6833
g5510.t1	OTHER	0.999108	0.000016	0.000876	
g5616.t1	SP	0.000067	0.999933	0.000000	CS pos: 16-17. ACA-AN. Pr: 0.5270
g5641.t1	SP	0.000130	0.999869	0.000001	CS pos: 16-17. ACA-AS. Pr: 0.4873
g5927.t1	OTHER	0.999995	0.000001	0.000004	
g702.t1	SP	0.000347	0.999652	0.000001	CS pos: 16-17. ALA-AN. Pr: 0.8153
g7861.t1	OTHER	0.999975	0.000004	0.000022	
g8100.t1	OTHER	0.999955	0.000024	0.000021	
g8312.t1	OTHER	0.999930	0.000065	0.000004	
```
</details>

SP means Signal peptide and OTHER refers to uncategorized ones.

```bash
grep -c "OTHER" ./search_res/output_protein_type.txt
#22
```

Therefore, before the next round non-nuclear proteins from the first search and signal from the second should be excluded.

<details>
    <summary>The list of remaining IDs is:</summary>
```
>g10513.t1
>g10514.t1
>g11513.t1
>g11806.t1
>g11960.t1
>g12510.t1
>g14472.t1
>g15484.t1
>g16318.t1
>g16368.t1
>g5443.t1
>g5927.t1
>g7861.t1
>g8100.t1
>g8312.t1
```
</details>

There are 15 of them in hits_nuclear.fa prepared for the next steps.

## Step 4. Function prediction

It is possible to predict protein functions even if no orthologous sequences from blast were found.

Perform the search with [HMMER](https://www.ebi.ac.uk/Tools/hmmer/) web-tool on the filtered hits.
Pfam 37.2 was used within hmmscan window.

Results are available [here](https://www.ebi.ac.uk/Tools/hmmer/results/899ffef8-e533-4c36-ad59-72125da6ef66/score)

The remaining 9 proteins without described function are:

```
>g10513.t1
>g10514.t1
>g11806.t1
>g12510.t1
>g14472.t1
>g16318.t1
>g16368.t1
>g5443.t1
>g5927.t1
```

## Step 5. BLAST search

Perform BLAST search on the remaining sequences.

![BLAST search](images/blast_search.png)

Search in UniProtKB/Swiss-Prot database via blastp algorithm.
Use name for species *Ramazzottius varieornatus* (optionally).

Only one result found:

| Protein ID  | Accession number | E-value | % Ident | % Query coverage | Score | Annotation |
| :-------: | :---: | :----: | :----: |:----: | :----: | :----: |
| g14472.t1 | P0DOW4.1 | 0 | 100.000 | 100.000 | 814 | Damage suppressor protein |

## Integrated results

Overall results on 34 candidate sequences are placed within summary.csv table

## Discussion

    In our search for the nuclear proteins taking part in DNA protection, we have analyzed 
34 protein candidates. After the first round of search by physical localisation we have filtered 17 candidates present in nucleus.
The next round of prediction by function helped to cut off another 2 candidates. Finally in BLAST we have identified the only protein g14472.t1 aligned with tardigrade protein (100% Ident, e-value = 0) and this protein is actually a damage suppressor.
This protein was reported earlier as [Dsup1](https://www.nature.com/articles/ncomms12808). 

The discussion about tardigrade survivability in extreme conditions has revealed many protective mechanisms. In one of the
most recent studies specific amplifications of several genes, including MRE11 and XPC, and numerous missense variants exclusive of R.varieornatus in CHEK1, POLK, UNG and TERT were described ([D. Carrero et al.](https://www.nature.com/articles/s41598-019-51471-8)). These genes combined influence DNA reparation and increase overall radiotolerance. How tardigrades gained these protection mechanisms? Evidence of horizontal gene transfer (HGT) in a tardigrade genome (17.5% of genes have foreign origin) was found in a freshwater species Hypsibius dujardini [T. Boothby et al.](https://doi.org/10.1073/pnas.1510461112). However, discussion on this topic has a strong counterargument about cotaminating sequences in the samples. The size difference between estimated genome size (around 100 Mbp) and assembly (212,3 Mbp) was highlighted, and the size of foreign genes is believed to be around 1-2% [G. Koutsovoulos et al.](https://doi.org/10.1073/pnas.1600338113). The discussion still continues, however, evidence of HGT were found recently, DODA1 gene is involved in radiotolerancy of waterbears and this gene is present in genomes of fungi and bacteria. In the same article, ardigrade-specific radiation-induced disordered protein (TRID1) and two mitochondrial pathways were described as part of radioprotection system of H. dujardini [Lei Li et al.](https://doi.org/10.1126/science.adl0799). 






