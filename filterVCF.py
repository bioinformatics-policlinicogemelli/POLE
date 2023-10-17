#funzione per filtrare le colonne ALT e FILTER di un file vcf e creare un nuovo
#vcf filtrato

import argparse
import pandas as pd
import vcf

def parsing_vcf(file_input,file_output):
    
    #lista contenente i dati da #CHROM in poi (no meta-information)
    with open (file_input) as f:
        correct_data=[line.strip().split("\t") for line in f if not line.startswith("##")]

    #creazione dataframe contenente i dati
    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])
    
    #filtro per dati in cui ALT!="."
    data=data[data["ALT"]!="."]

    #filtro per dati in cui FILTER == "PASS"
    data=data[data["FILTER"]=="PASS"]
    # print(type(data))
    #per salvare il file così filtrato
    data.to_csv(file_output,sep="\t", mode="a", index=False)


def write_vcf_meta(file_path, metadata):
    # input file path output - output write meta data file
    with open(file_path, 'w') as meta_file:
            for key, value in metadata.items():
                if type(value) is str:
                    meta_file.write('##'+ key + ': ' + value + '\n')
                elif type(value) is list:
                    meta_file.write('##'+ key + ': '.join(value) + '\n')


def read_vcf_meta(file_path):
    # input vcf file path - return vcf metadata
    with open(file_path, 'r') as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)
        metadata = vcf_reader.metadata
    return metadata
    
    

def main(INPUT, OUTPUT):
    metadata = read_vcf_meta(INPUT)

    write_vcf_meta(OUTPUT, metadata)
    parsing_vcf(INPUT,OUTPUT)

if __name__ == '__main__':


# parse arguments
    parser = argparse.ArgumentParser(description="Filter for a vcf",
                                            epilog="Version: 1.0\n\
                                            Author: Bioinformatics Facility GSTeP'\n\
                                            email: luciano.giaco@policlinicogemelli.it")

    # arguments
    parser.add_argument('-i', '--input', help="<input.vcf>\
                                            VCF file to filter",
                                            required=True)
    parser.add_argument('-o', '--output', help="<output-file.tab>\
                                            file path of the Table output",
                                            required=True)

    args = parser.parse_args()
    INPUT = args.input
    OUTPUT = args.output

    main(INPUT, OUTPUT)




''' possibili ulteriori funzioni
#ecco una possibile nuova funzione che prende come input una delle chiavi del dizionario (quindi una mutazione) e ti ritorna da solo il valore associato (quindi il numero di mutazioni osservate)
def select(vcf):
    vcf=dizionario(vcf)
    print ("Scegli la mutazione da vedere:")
    richiesta= input()
    for key, value in vcf.items():
        if key == str(richiesta):
            print('The amount of', str(richiesta), 'is: ', value)#per i valori presenti nel dizionario, ti ritorna il valore corrispendete alla chiava data come input
            
#in caso parlane con Luciano potrebbe essere figa da mettere

#------------------------------


#possibile nuova funzione per comparare 2 dizionari magari presi da 2 vcf
# basato sulla richiesta del dottor. Angelo Minucci che ha chiesto se è possibile trovarne le differenze tra 2 file di input

def compare(file1, file2):
    file1=dizionario(file1)
    file2 = dizionario(file2)
    valori_uguali = {k: file1[k] for k in file1 if k in file2 and file1[k] == file2[k]}
    print(len(valori_uguali))
#finora questo da i valori uguali presenti nei 2 vcf, sia per le key che per i values. a noi frega poco se sono uguali, dato che è poco probabile che accade una cosa simile.
    valori_differenti = {k: file1[k] for k in file1 if k in file2 and file1[k] != file2[k]}
    print(len(valori_differenti))



#TMB(Tumor mutational burden)= The total number of mutations (changes) found in the DNA of cancer cells
#per trovarlo devi fare le mutazioni sul totale
#TMB > 100 mut/Mb.
#nell'articolo "we used 38 Mb as the estimate of the exome size"
#For comparison, 321 microsatellite-stable (MSS)

all 321 MSS, POLE wild-type ECs scored ≤3 points and all 127 MSI-H POLE wild-type ECs scored ≤2.
COSMIC signature 10 was identified in 11 (8.7%) MSI-POLEwt and 96 (29.9%) MSS-POLEwt ECs (mean signature 10 contribution 0.002, range 0.000–0.048, and mean 0.017, range 0.000–0.218, respectively).


aldilà dei microsatelliti, l'articolo parla di "recurrent variant in EC" quindi chiede se ci sono immagino delle mutazioni più frequenti rispetto ad altre
quindi, MAGARI, lo score deve aumentare di +1 se per ipotesi una dei tipi di mutazioni super un certo numero nello specifico
si ma quale. controlla sull'articolo

def findsatellites(file):
    lista_ALT=[] #i microsatelliti sono visibili nel ALT come quei valori che sono sequenze tipo AATTCCGG e simili
    import vcf
    vcf_reader = vcf.Reader(open(file, 'r')) #leggi il vcf di input
    for record in vcf_reader:
            lista_ALT.append(str(record.ALT)) #mettimi nella lista tutti i valori ALT
            
    lunghezza_minima_microsatelliti = 4
    #lunghezza minima presente nei dati di ALT per trovare i valori dei microsatelliti
    #qui sotto la funzione inserisce in una nuova lista tutte le ALT che non hanno la lunghezza data
    lista_microsatelliti=[]
    lista_microsatelliti=([x for x in lista_REF if len(x) >= lunghezza_minima_microsatelliti])
    

'''
