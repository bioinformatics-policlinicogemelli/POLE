# PathoPoleAnalyzer 

PathoPoleAnalyzer is a tool for calculating a pathogenicity score for Pole mutations, inspired by the findings presented in the Castillo's article (Castillo et. all, 2020)  ["Interpretation of somatic POLE mutations in endometrial carcinoma"] (https://pubmed.ncbi.nlm.nih.gov/31829442/). 

The score is an integer ranging from 0 to 6, indicating whether the POLE mutation is pathogenic (Score >=4), of uncertain significance (Score = 3), or non-pathogenic (Score<3),based on the presence of specific genomic alterations.

To calculate this score, PathoPoleAnalyzer evaluates several parameters:
- Tumor Mutational Burden (TMB)
- Presence of recurrent Pole mutations
- Frequency of C>A substitutions
- Frequency of C>G substitutions
- Frequency of T>G substitutions
- Frequency of Indels (insertions and deletions)

## Requires
- python == 3.9
- PyVCF == 0.6.8
- numpy == 1.23.5
- pandas == 1.5.2

## Input 


### Filter vcf
To reduce the time required for the analysis, it is suggested to first filter the vcf file.
This can easily be done with the [vcf_filter.py]([https://gitlab.com/gstep-bioinformatics-core-facility-research/varan-2.0/-/blob/main/vcf_filter.py](https://github.com/bioinformatics-policlinicogemelli/POLE/blob/main/filter_VCF.py)) script.
Use the following command:
>$ python filterVCF.py -i vcfFilePath -o destinationFilteredVcf

### Input Data
PathoPoleAnalyzer requires a VCF file and a TMB value as input.



|COMMAND OPTION                |DESCRIPTION                          |TYPE                         |REQUIRED                         |EXAMPLE
|----------------|-------------------------------|-----------------------------|-----------------------------|-----------------------------|
|-f <br>--output_folder| <p align="justify">Add this option to insert the path where to save the output folder| String | Yes | /Path/to/file
|-t <br>--TMB| <p align="justify">Add this option to insert the path to the input tsv file where the vcf files are listed| Integer | Yes | 100
|-c <br>--Category| <p align="justify">Category of mutations| String | No | "hotspot", "wt", "exo" or "vus"


Use the following command to receive the full output of the script:
>$ python pathopoleanalyzer.py -f "YOUR_filtered_VCF_FILE.vcf" --TMB (value of the TMB)
or
>$ python pathopoleanalyzer.py -f "YOUR_filtered_VCF_FILE.vcf" -t (value of the TMB)

You can save the result of your analysis in a text file by typing the following command in the Command Prompt:
>$ pathopoleanalyzer.py -f "YOUR_filtered_VCF_FILE.vcf" > VCF_Results.txt

If you prefer to have an output featuring additional details regarding mutation frequency and the conditions which raised the score, try to run the 'POLE_SCORE(additional_output_details).py' with the same process.


## What raises the Score?

The conditions for which the Score is raised are based on the data available on ["Interpretation of somatic POLE mutations in endometrial carcinoma"](https://pubmed.ncbi.nlm.nih.gov/31829442/)(Castillo et. all, 2020).
Takin into consideration the different size of the exome of the FPG500 input data featured in our study, the threshold the allowed a raise of the score were modified with a conversione via Median of Medians technique-.
Given a filtered VCF file, the requirement to find a POLE mutation are:

	C>A substitutions > 6.6%, 
	T>G substitutions > 4%, 
	C>G substitutions < 4.6%, 
	Indels < 4%, 
	TMB > 100 mut/Mb.
	Presence of Recurrent Mutations


<div class="center">

  

```mermaid

graph TD
A[if C>A over 6.6%] --> B(Score +1)
F[if T>G over 4%] --> B(Score +1)
G[if Indels below 4%] --> B(Score +1)
H[if C>G below 4.6%] --> B(Score +1)
I[if TMB over 100 mut/Mb] --> B(Score +1)
J[if Recurrent Mutations] --> B(Score +1)
K[if C>T over 30%] --> B(Score +1)
B --> C(If Score >=4:) --> L(Pathogenic POLE mutation)
B --> D(If Score =3:) --> M(Variant of Unknown Significance)
B --> E(If Score <3:) --> N(Non-pathogenic POLE mutation)
```





## What is defined as a 'Recurrent Mutation'?

The [over-mentioned Paper](https://pubmed.ncbi.nlm.nih.gov/31829442/)(Table 1 and Table 3 of the Paper) takes in consideration whether POLE mutations were recurrent in ECs within the COSMIC or TCGA databases, as recurrent mutations are more likely to be pathogenic.

As mentioned before, if at least one of these next mutations are featured in the analized VCF, the score will rise by +1.

Below the mutations that are considered recurrent in POLE, along wIth chromosomal position, nucloetide substitution and protein change:

	CHROM: POS			NUCL. SUB.		PROT. CHANGE
	chr12:133253184		c.857C>G		P286R
	chr12:133250289 	c.1231G>T		V411L
	chr12:133253151		c.890C>T  		S297F
	chr12:133249847		c.1376C>T		S459F
	chr12:133249857		c.1366G>C		A456P
	chr12:133252327		c.1100T>C		F367S
	chr12:133250250		c.1270C>A		L424I
	chr12:133253157		c.884T>G		M295R
	chr12:133250213		c.1307C>G		P436R
	chr12:133250189		c.1331T>A		M444K
	chr12:133252325		c.1102G>T		D368Y
 	chr12:133250250		c.1270C>G		L424V
	chr12:133253208		c.833C>T		T278M
 	chr12:133249829		c.1394C>T		A465V
	chr12:133256623      	c.340C>T		R114*
	chr12:133252023		c.1187A>G		E396G
	chr12:133257828		c.100C>T		R34C


</div>
