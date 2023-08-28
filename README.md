# POLE
Repository for pathogenic estimation of POLE mutations via development of a score


# POLE functions

## The Script

The script overall generates a Dictionary that is composed of the counts of types of mutations found in a filtered VCF file given as input.

After this it generates a score: an integer number whose value depends on the data of the dictionary and returns a specific output depending on the score. This putput is mostly a sentence which suggests if the vcf presents a POLE mutation which is pathogenic (Score >= 4), non-pathogenic(Score < 3) or  a variant of unknown isgnificance (score = 3)
To generate this score are taken into account not only the data of the Dictionary, but also the list Indels (which comprehend the amoun to observes events of insertion and mutation) and the total number of different mutations observed in the VCF file.

The output files of this script are based on the data available on the following article  ["Interpretation of somatic POLE mutations in endometrial carcinoma"](https://pubmed.ncbi.nlm.nih.gov/31829442/). 
(Castillo et. all, 2020)



## Input Files and Requirements
All the functions of the code take a file / folder as an input
The required input files are .vcf
It is suggested to first filter the vcf that needs to be analyzed.
This can easily be done with the [vcf_filter.py](https://gitlab.com/gstep-bioinformatics-core-facility-research/varan-2.0/-/blob/main/vcf_filter.py) script.
Use the following command, which requires a file to input and the preferred path for the output file:
>$ python filter_VCF.py -i "path the the VCF file" -o "name of the output VCF file"






## Functions of the Script

In total the Script has 6 functions, shown in the table below:


|FUNZIONE                |INPUT                          |OUTPUT                         |	EXAMPLE OUTPUT
|----------------|-------------------------------|-----------------------------|-----------------------------|
|'dizionario'						|VCF			|A Dictionary with the number of observed events of substitutions. Keys are the types of mutatioons; Values are the number of said mutations observed.			|{'A>C': 0, 'A>G': 13, 'A>T': 2, 'C>A': 1, 'C>G': 1, 'C>T': 8, 'G>A': 14, 'G>C': 0, 'G>T': 2, 'T>A': 1, 'T>C': 8, 'T>G': 0}
|'listaindels'        			|VCF           |List including both insertion and deletion events ("Indels").            |	['[CT]', '[AT]', '[AT]', '[GC]', 'TC', 'CG']
|'recurrentmutations'        			|VCF           |Dictionary  with the "Recurren Mutations" found in the VCF; Dictionary includes chromosome positions as the keys and mutations as the values.            |	 Numero Recurrent Mutations trovate:  2, Recurrent Mutations =  {12345678: 'C>G', 11223344: 'C>A'}
|'eventimutazionetotali'          |VCF			|Total number of mutations, including substitutions, deletions and insertions.| 12345
|'valorepercentuale'          	|VCF				|Percentage frequency of events of mutation compared to the total VCF data for every event of substitution.| A>G mutations = 26.0%; C>T mutations = 3.0 %, ....
||||
|'calcoloscore'          			|Dictionary, TMB, MSI |Analyzing the given data, the function will increase the score (which starts at "0") of +1 each time one of the given conditions are satisfied(ex. percentage of frequency of events 'C>T' greater than 30%). The output is an integer which represents the Score.	Depending on the value of the Score obtained with the function above, a different comment will pop up as an estimation of the type of observed mutation.	| 3, 4 or 5...... example of the message printed: "Pathogenic POLE mutation"
||||
|'dictoplot'          	|VCF				|Percentage frequency of events of mutation in the form of a Bar Plot. Additional Bar Plot analysis relegated to comparison of the mutation betweeen two data (one 'OLD' and one 'NEW'; the 'NEW one given by degault')| Bar Plot including mutations events, both substitutions and indels.


## Diagram of the Workflow

UML flow chart available made with [Mermaid](https://mermaidjs.github.io/). 
<div class="center">

```mermaid
graph TD
A[VCF INPUT FILE] --> B((Dizionario))
A --> D((Tot. Mutazioni))
A --> C((ListaIndels))
A --> R((RecurrentMutations))
B --> F(ValorePercentuale)
B -- Dictionary--> E{CalcoloScore}
C -- List --> E
R -- Dictionary --> E
D -- Integer --> E --> Score+Comment
```
</div>




## Call the Script in the Command Prompt

An ArgParse has been implemented which allows to call the function in the Terminal.
Make sure to enter the folder with both the Python script and your filtered VCF File.
The required input data include the folder where the filtered VCF is positioned and the TMB (Tumor Mutational Burden)

Use the following command to receive the full output of the script:
>$ python script_dictionary_e_score.py -f "YOUR_filtered_VCF_FILE.vcf" --TMB (value of the TMB)
or
>$ python script_dictionary_e_score.py -f "YOUR_filtered_VCF_FILE.vcf" -t (value of the TMB)

You can save the result of your analysis in a text file by typing the following command in the Command Prompt:
>$ python script_dictionary_e_score.py -f "YOUR_filtered_VCF_FILE.vcf" > VCF_Results.txt



## What raises the Score?

The conditions for which the Score is raised are based on the data available on ["Interpretation of somatic POLE mutations in endometrial carcinoma"](https://pubmed.ncbi.nlm.nih.gov/31829442/)(Castillo et. all, 2020).
Given a filtered VCF file, the requirement to find a POLE mutation are:

	C>A substitutions > 20%, 
	T>G substitutions > 4%, 
	C>G substitutions < 0.6%, 
	Indels < 5%, 
	TMB > 100 mut/Mb.
	Presence of Recurrent Mutations

Analysis of the C>T mutations (>= 30%) has been added in reference to the following paper:
["The mutational signatures of formalin fixation on the human genome | Nature Communications"](https://www.nature.com/articles/s41467-022-32041-5).
(Graham et. al, 2022).

<div class="center">

  

```mermaid

graph TD
A[if C>A over 20%] --> B(Score +1)
F[if T>G over 20%] --> B(Score +1)
G[if Indels below 20%] --> B(Score +1)
H[if C>G below 0.6%] --> B(Score +1)
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
	chr12:132676598		c.857C>G		P286R
	chr12:132673703 	c.1231G>T		V411L
	chr12:132676565		c.890C>T  		S297F
	chr12:132673261		c.1376C>T		S459F
	chr12:132673271		c.1366G>C		A456P
	chr12:132675741		c.1100T>C		F367S
	chr12:132673664		c.1270C>A		L424I
	chr12:132676571		c.884T>G		M295R
	chr12:132673627		c.1307C>G		P436R
	chr12:132673603		c.1331T>A		M444K
	chr12:132675739		c.1102G>T		D368Y



</div>
	
# Filter ClinVar

1. Si sono filtrate prima tutte le benigne, in questo caso rimangono solo le VUS che verranno caricate sul cbioportal, il data_mutation_extended.txt file è nella rispettiva cartella della neoplasia denominata NoBenign

2. Il secondo filtro rimuove tutte le VUS a partire dai maf presenti nella cartella NoBenign da cui sono state eliminate le benigne. In questo caso sono tenute le occorrenze **pathogenic, likely_pathogenic, risk_factor, drug_response** Il data_mutation_extended.txt file è nella rispettiva cartella della neoplasia denominata NoVus

## workflow dei comandi

**Filtro su benigne**

```
python3 filter_clinvar.py \
-f /data/hpc-share/cbioportal/ovary/maf/ \
-o /data/hpc-share/cbioportal/ovary/
```

**Creato file da caricare su cBioportal**

```
python3 concatenate.py \
-f /data/hpc-share/cbioportal/ovary/NoBenign/ \
-e maf \
-o /data/hpc-share/cbioportal/ovary/NoBenign/data_mutations_extended.txt
```

**Filtrato da NoBenign le VUS**
```
python3 filter_clinvar.py \
-f /data/hpc-share/cbioportal/ovary/NoBenign/ \
-o /data/hpc-share/cbioportal/ovary/ \
-v
```

**Creato file da caricare su cBioportal**

```
python3 concatenate.py \
-f /data/hpc-share/cbioportal/ovary/NoVus/ \
-e maf \
-o /data/hpc-share/cbioportal/ovary/NoVus/data_mutations_extended.txt
```

**Creato le tabelle patient e sample clinical data**

```
python3 walk.py -i SAMPLE_TOTAL_OVARY.tsv 	
-t snv 	
-m 	
-o /data/hpc-share/cbioportal/ovary/
```

**Creato tabella fusioni**

```
python3 walk.py \
-i SAMPLE_TOTAL_OVARY.tsv \
-t snv \
-u \
-o /data/hpc-share/cbioportal/ovary/
```


## Datasets 


- [ ] breast
- [ ] colangio
- [x] colon
- [x] endometrium
- [ ] gist
- [x] lung
- [ ] melanoma
- [x] ovary
- [ ] pancreas
- [ ] prostata
- [ ] tiroide







