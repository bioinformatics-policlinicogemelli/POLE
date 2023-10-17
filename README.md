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


Use the following command to launch PathoPoleAnalyzer:
>$ python pathopoleanalyzer.py -folder vcfFilePath --TMB TMBvalue
or
>$ python pathopoleanalyzer.py -f vcfFilePath -t TMBvalue

You can save the result of your analysis in a text file, using the option ```-o/--output``` 
>$ pathopoleanalyzer.py -f vcfFilePath -t TMBvalue -o outputFile

## What raises the Score?

The thresholds for the different parameters are readjusted compared to Castillo's article and are as follows:

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

According to [over-mentioned Paper](https://pubmed.ncbi.nlm.nih.gov/31829442/)(Table 1 and Table 3) :

	Recurrent mutations were defined as those present in two or more cancer samples in the COSMIC and TCGA databases 
 	combined (cases present in both databases were counted only once). 
  	A mutation was considered non-recurrent if it was found only one.


Below the mutations that are considered recurrent in POLE, along with chromosomal position, nucleotide substitution and protein change:

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
