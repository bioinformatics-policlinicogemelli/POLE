#Script che genera un dizionario che conta i tipi di mutazioni presenti in un VCF filtrato
#Dopo questo genera uno score in base ai risultati presenti nel dictionary e ti riporta varie voci di output in base al risultato dello score.
#le richieste di output presenti in questo script sono basate sull'articolo "Interpretation of somatic POLE mutations in endometrial carcinoma" (Castillo et. all, 2020)

import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 

def printheader(file):
    vcf_reader = vcf.Reader(open(file, 'r'))
    for line in vcf_reader:
        if line.startswith('##'):# This is a metadata line
            header = line.strip()
        #elif line.startswith('#'):# This is the header line
        #    header = line.strip()
        return header


#il risultato della funzione def dizionario(file) riporta un dizionario che in cui le mutazioni sono le keys e i valori sono la frequenza in cui queste mutazioni appaiono nel VCF
def dizionario(file):
   
    #formazione del dizionario vuoto
    mutations= ('A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G')
    dictionary=dict.fromkeys(mutations)
    number=0
    dictionary=dict.fromkeys(mutations,number)
    # il dizionario creato qui sarà {'A>C': 0, 'A>G': 0, 'A>T': 0, 'C>A': 0, 'C>G': 0, 'C>T': 0, 'G>A': 0, 'G>C': 0, 'G>T': 0, 'T>A': 0, 'T>C': 0, 'T>G': 0}
    vcf_reader = vcf.Reader(open(file, 'r'))
    #per ogni possibile mutazione aggiungiamo i dati: se "valore in REF" = 'A' e "valore in INT" = 'C', aggiungi un 1 alla key specifica; stessa condizione ripetuta per ogni possibile evento di sostituzione osservabile
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
           
        # Mutazioni ora NON sono considerate doppie se entrambi gli alleli sono mutati, aggiungiamo +1 a prescindere
            if record.REF == 'A' and record.ALT[alt]=='C':
                dictionary['A>C']+=1
            if record.REF == 'A' and record.ALT[alt]=='G':
                dictionary['A>G']+=1
            if record.REF == 'A' and record.ALT[alt]=='T':
                dictionary['A>T']+=1    
            if record.REF == 'C' and record.ALT[alt]=='A':
                dictionary['C>A']+=1
            if record.REF == 'C' and record.ALT[alt]=='G':
                dictionary['C>G']+=1
            if record.REF == 'C' and record.ALT[alt]=='T':
                dictionary['C>T']+=1
            if record.REF == 'G' and record.ALT[alt]==['A']:
                dictionary['G>A']+=1
            if record.REF == 'G' and record.ALT[alt]=='C':
                dictionary['G>C']+=1
            if record.REF == 'G' and record.ALT[alt]=='T':
                dictionary['G>T']+=1      
            if record.REF == 'T' and record.ALT[alt]=='A':
                dictionary['T>A']+=1
            if record.REF == 'T' and record.ALT[alt]=='C':
                dictionary['T>C']+=1
            if record.REF == 'T' and record.ALT[alt]=='G':
                dictionary['T>G']+=1
    return dictionary



#funzione che crea un dizionario con posizioni relevant come keys e il tipo di mutazione in quella posizione, dati basati sulla tabella del analisi di POLE nel articolo: "Interpretation of somatic POLE mutations in endometrial carcinoma"
def recurrentmutations(file):
    #positions=[]
    #mutations=[]
    vcf_reader = vcf.Reader(open(file, 'r'))
    totalmuts={} #dictionary che includerà POS come keys e REF>ALT come values
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if record.CHROM == 'chr12':
                x=record.POS
                y= str(record.REF)+">"+str(record.ALT[alt])
                totalmuts[x] = y
                #positions.append(x) #aggiungi ogni POS a una lista
                #mutations.append(y) #aggiungi ogni REF+ALT a una lista   
    #print("Dictionary Mutazioni = ", totalmuts, "\n")
    #print("Lunghezza dict Mutazioni = ", len(totalmuts), "\n")

    #print("Posizioni delle mutazioni = ", positions, "\n", "Numero Posizioni = ", len(positions), "\n")
    #print("Mutazioni trovate = ", mutations, "\n", "Numero Mutazioni = ", len(mutations), "\n")
    
    #questa è la lista di mutazioni che secondo l'articolo sono definibili come "Recurrent". Le keys sono le posizioni trovate tramite sito Varsome per le rispettive mutazioni, e i values sono gli eventi di mutazione chiave
    #se tali valori coincidono con quelli visti nel VCF filtrato, allora ci sono "current mutations in EC", condizione che aumenta di +1 lo score generale

    ##########ATTENZIONE
    # dictionary errato dei recurrent, dato che usa come riferimento genoma hg38:
    #listrecurrent= {132676598: 'G>C', 132673703: 'C>A', 132673622: 'C>A', 132676565: 'G>A', 132676143: 'G>A', 132673261: 'G>A', 132673180: 'G>A', 132673271: 'C>G', 132675741: 'A>G', 132673664: 'G>T', 132673583: 'G>T', 132676571: 'A>C', 132676149: 'A>C', 132673627: 'G>C', 132673249: 'G>C', 132673603: 'A>T', 132673225: 'A>T', 132675739: 'C>A', 132675441: 'C>A', 132668416: 'G>A'}
    #nuova list recurrent per hg19
    
    
    
 #(OLD)listrecurrent = {133253184: 'G>C', 133250289: 'C>A', 133250208: 'C>A', 133253151: 'G>A', 133249349: 'G>A', 133252729: 'G>A', 133249847: 'G>A', 133249766: 'G>A', 133249857: 'C>G', 133252327: 'A>G', 133250250: 'G>T', 133250169: 'G>T', 133253157: 'A>C', 133249355: 'A>C', 133225944: 'A>C', 133252735: 'A>C', 133250213: 'G>C', 133249835: 'G>C', 133250189: 'A>T', 133249811: 'A>T', 133252325: 'C>A', 133248833: 'C>A', 133252027: 'C>A'}
    listrecurrent = {133253184: 'G>C', 133250289: 'C>A', 133253151: 'G>A', 133249847: 'G>A', 133249857: 'C>G', 133252327: 'A>G', 133250250: 'G>T', 133253157: 'A>C', 133250213: 'G>C', 133250189: 'A>T', 133252325: 'C>A', 133250250: 'G>C', 133250238: 'C>T', 133250250: 'G>C', 133253208: 'G>A', 133249829: 'G>A', 133256623: 'G>A', 133252023: 'T>G', 133257828: 'G>A'}
    
    #valori alternativi di V411L 132673703: 'C>G', 132673622: 'C>G' da aggiungere una volta risolto il problema del dictionary: tali chiavi annullerebbero i valori precendentemente assegnati: inserirli direttamente non va bene perchè non ci possono essere duplicati nelle chiavi di un dizionario
    # il dictionary qui include tutti valori presi da Varsome, sia quelli NM_006231.4 che quelli ENST00000535270
    #print("Mut ricorrenti = ", listrecurrent)
    shared_items = {k: listrecurrent[k] for k in listrecurrent if k in totalmuts and listrecurrent[k] == totalmuts[k]}
    if len(shared_items) == 0:
     #   print("     No Recurrent Mutations found", "\n") 
        return shared_items   
    else:
      #  print("     Numero Recurrent Mutations trovate: ", len(shared_items)) 
        return shared_items



#funzione che ritorna la lista di eventi di inserzione e delezione
def listaindels(file):
    lista_indels=[]
    lista_inserzioni=[]
    lista_delezioni=[]
    #lista_doppietti=[]
    vcf_reader = vcf.Reader(open(file, 'r')) #leggi il vcf di input
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            #if len(record.REF)!=1:                         # Debug: Capire come mai c'erano degli eventi in più. Risposta: Dinucleotidi in ALT e REf mal classificati nel 1°metodo
             #   if len(record.REF)==len(record.ALT[alt]):
              #      lista_doppietti.append(record.REF)
            if len(record.ALT[alt])!=len(record.REF):
                lista_indels.append(record.REF)
                if len(record.ALT[alt])>len(record.REF):
                    lista_inserzioni.append(record.ALT[alt])
                if len(record.ALT[alt])<len(record.REF):
                    lista_delezioni.append(record.REF)
    
    return lista_indels



#numero totale di elementi nel vcf filtrato, siano essi indels o non indels.
def eventimutazionetotali(file):
    totalemutazioni=0
    #totalemutazioni è il numero totale di elementi nel vcf filtrato, siano essi indels o non indels.    
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
            totalemutazioni+=1

    return totalemutazioni



######calcolo frequenza BED file
'''       
def analyze_fasta(file_path):
    # Initialize variables
    seq_count = 0
    seq_lengths = []
    a_count = 0
    c_count = 0
    g_count = 0
    t_count = 0
    # Open the file and read its contents
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # This is a header line
                seq_count += 1
            else:
                # This is a sequence line
                seq_lengths.append(len(line.strip()))
                total_length += len(line.strip())
                a_count += line.count('A')
                a_count += line.count('a')
                c_count += line.count('C')
                c_count += line.count('c')
                g_count += line.count('G')
                g_count += line.count('g')
                t_count += line.count('T')
                t_count += line.count('t')
    # Calculate statistics
    min_length = min(seq_lengths)
    max_length = max(seq_lengths)
    avg_length = total_length / seq_count
    a_content = a_count / total_length * 100
    c_content = c_count / total_length * 100
    g_content = g_count / total_length * 100
    t_content = t_count / total_length * 100
    # Print the results
    print(f"Number of sequences: {seq_count}")
    print(f"Minimum length: {min_length}")
    print(f"Maximum length: {max_length}")
    print(f"Average length: {avg_length:.2f}")
    print(f"Total length: {total_length}")
    print(f"A content: {a_content:.2f}%")
    print(f"C content: {c_content:.2f}%")
    print(f"G content: {g_content:.2f}%")
    print(f"T content: {t_content:.2f}%")
'''

#ADESSO CALCOLA LO SCORE, cioè la frequenza delle specifiche mutazioni
def valorepercentuale(file):
    dict = dizionario(file)
    totalsos = sum(dict.values())
    #totalsos è la conta totale degli eventi di sostituzione osservati (eventi come gli indels sono scartati da questa conta)
    listpercent = []
   # print('PERCENTUALE DELLE MUTAZIONI OSSERVATE RISPETTO AL TOTALE: ')
    for key, value in dict.items():
        percent = round((int(value) / totalsos),2)*100 #frequenza percentuale "classico" 
        #percent2 = percent*1.98/ 38 # frequenza percentuale modificata in relazione alle proporzioni usate da Castillo; 38 : frequenza_castillo = 1,98 : frequenza_TSO500
    #########    if user_input.lower() == 'exo':
       # print("Frequenza in Percentuale di mutazione:", key, "=", percent, "%")
        listpercent.append(percent)
    #########    elif user_input.lower() == 'panel':
    #########        print("Frequenza in Percentuale di mutazione:", key, "=", percent2, "%")
    #########        listpercent.append(percent2)
    #return listpercent
    return ' '



#CALCOLO SCORE IN BASE AI NOSTRI PARAMETRI PRESI DA DICTIONARY (LISTA_INDELS E NUMERO DI EVENTI TOTALI)
#def calcoloscore(file, TMB):
def calcoloscore_CASTILLO(file, TMB):
  
    score=0
    
    # print ("CALCOLO DELLO SCORE", "\n")
    # print("A seguire sono mostrate le condizioni rispettate, ognuna di esse aumenta lo Score di +1.")
    # print("(La mancanza di Recurrent Mutations non porta a un aumento dello score)", "\n")
    dictionary = dizionario(file)
    lista_indels= listaindels(file)
    totalemutazioni = eventimutazionetotali(file)
    totalsos = sum(dictionary.values())
    #totalsos anche qui è la conta totale degli eventi di sostituzione osservati (eventi come gli indels sono scartati da questa conta)
    

    
    #if int(TMB)>= 100: #TMB over 100mut/Mb
    #    print ("+    Tumour Mutational Burden maggiore di 100 mut/Mb. TMB =", TMB, "\n")
    #    score +=1
    #if int(MSI) >= ((int(totalemutazioni))*30)%100: # MSI-H risultano da valori di MSI >=30% di loci MSI instabili (>2 or more of the 5 loci)
       # print ("+    DNA Microsatellite Instability (MSI) Detected.   MSI site = ", MSI, "\n")
       # score +=1
####user_input = input('Seleziona modalità di analisi (exo / panel): ')    
####if user_input.lower() == 'panel':
        ###############################
        ##modificate le percentuali
        ###############################

    #print('Percentuali utilizzate per il calcolo dello score modificate rispetto al paper di Castillo et al.', '\n', 'Modifica in base alla proporzione: [38 : frequency_Castillo = 1,98 : frequency_TS0500]', '\n')
    if float(TMB)>= 100: #TMB over 100mut/Mb
      #  print ("+    Tumour Mutational Burden maggiore di 100 mut/Mb. TMB =", TMB, "\n")
        score +=1
    for key, value in dictionary.items():

        if key=="C>A": #C>A over 20%

            if round((int(value) / totalsos),4)*100 >= 20: 
              #  print('+    Mutazioni C>A uguali o superiori al valore richiesto (maggiore del 20% in Castillo)')
                score +=1
                
            # print("                    Quantità di C>A: ", dictionary["C>A"])
            # print("                    Frequenza C>A: ", round((int(value) / totalsos),4)*100, " %")
            # print("                    Frequenza C>A rispetto alle mutazioni totali: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            # print("\n")   
            ca_perc=round((int(value) / totalsos),4)*100 
        if key=="T>G": #T>G over 4%

            if round((int(value) / totalsos),4)*100 >= 4:
           #     print('+    Mutazioni T>G uguali o superiori al valore richiesto (maggiore del 4% in Castillo)')
                score +=1
                
            # print("                    Quantità di T>G: ", dictionary["T>G"])
            # print("                    Frequenza T>G: ", round((int(value) / totalsos),3)*100, " %")
            # print("                    Frequenza T>G rispetto alle mutazioni totali: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            # print("\n")  
            tg_perc=round((int(value) / totalsos),4)*100 
        if key=="C>G": #C>G below 0.6%
            if round((int(value) / totalsos),4)*100 <= 0.6: 
           #     print('+    Mutazioni C>G uguali o inferiori al valore richiesto (minore del 0.6% in Castillo)')
                score +=1
                
            # print("                    Quantità di C>G: ", dictionary["C>G"])
            # print("                    Frequenza C>G: ", round((int(value) / totalsos),4)*100, " %")
            # print("                    Frequenza C>G rispetto alle mutazioni totali: ", round((int(value) / totalemutazioni), 3)*100, " %", "\n")
            # print("\n") 
            cg_perc=round((int(value) / totalsos),4)*100 
            
    ratio_indels=round((len(lista_indels)/ totalemutazioni),4) * 100
    #Qui abbiamo la frequenza delle indels rispetto al totale
    
    

    if ratio_indels <= 5:  # ##(0.26 = 5*1.98/ 38)  Indels below 5%
       # print('+    Quantità di indels uguale o inferiore al valore richiesto (minore del 5% in Castillo)','\n')
        score+=1
        
    # print ("        ###################################################")
    # print ("                    Numero indels =", len(lista_indels))    
    # print ("                    Ratio indels =", ratio_indels, "%")
    # print ("        ###################################################")
    # print ("\n")
    
    #elif user_input.lower() == 'exo':
            ###############################
            ##percentuali originali
            ###############################
     #   print('Percentuali utilizzate per il calcolo dello score basate sul paper di Castillo et al.', '\n')
      #  if int(TMB)>= 100: #TMB over 100mut/Mb
        #    print ("+    Tumour Mutational Burden maggiore di 100 mut/Mb. TMB =", TMB, "\n")
       #     score +=1        
        #for key, value in dictionary.items():

         #   if key=="C>A": #C>A over 20%
          ##      if round((int(value) / totalsos),2)*100 >= 20:
           #         print('+    Mutazioni C>A uguali o superiori al 20%','\n')
             #       score +=1
            #if key=="T>G": #T>G over 4%
              #  if round((int(value) / totalsos),2)*100 >= 4:
               #     print('+    Mutazioni T>G uguali o superiori al 4%','\n')
                #    score +=1
            #if key=="C>G": #C>G below 0.6%
             #   if round((int(value) / totalsos),2)*100 <= 0.6:
              #      print('+    Mutazioni C>G uguali o inferiori al 0.6%','\n')
               #     score +=1

       # ratio_indels=round((len(lista_indels)/ totalemutazioni),2) * 100
        #Qui abbiamo la frequenza delle indels rispetto al totale
       # if ratio_indels <= 5: #Indels below 5%
        #    print('+    Quantità di indels uguale o inferiore al 5%','\n')
         #   score+=1

    
    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
       # print(response, "\n")
        score+=1
    #se il dictionary recurrent non è vuoto, sono state identificate Recurrent mutations, perciò lo Score aumenta di +1

    
    # print('--------------------------', "\n")
    
    # print('Current score (CASTILLO) =', score,'\n')
    
    #Commento basato sul valore finale dello score
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return score,round(ca_perc,2), round(tg_perc,2),round(cg_perc,2), round(ratio_indels,2)


#CALCOLO SCORE IN BASE AI NOSTRI PARAMETRI PRESI DA DICTIONARY (LISTA_INDELS E NUMERO DI EVENTI TOTALI)
#def calcoloscore(file, TMB):
def calcoloscore_TS0500(file, TMB):
    score=0
    dictionary = dizionario(file)
    lista_indels= listaindels(file)
    totalemutazioni = eventimutazionetotali(file)
    totalsos = sum(dictionary.values())
    #totalsos anche qui è la conta totale degli eventi di sostituzione osservati (eventi come gli indels sono scartati da questa conta)

    if float(TMB)>= 5.21: #(5.21 = 100*1.98/ 38) #TMB over 100mut/Mb
        score +=1
    for key, value in dictionary.items():

        if key=="C>A": #C>A over 20%
            if round((int(value) / totalsos),4)*100 >= 1.0421: ## (1.0421 = 20*1.98/ 38)
                score +=1
               
        if key=="T>G": #T>G over 4%
            if round((int(value) / totalsos),4)*100 >= 0.208: ## (0.208 = 4*1.98/ 38)
                score +=1
            
        if key=="C>G": #C>G below 0.6%
            if round((int(value) / totalsos),4)*100 <= 0.0312: ## (0.0312 = 0.6*1.98/ 38)
                score +=1

    ratio_indels=round((len(lista_indels)/ totalemutazioni),4) * 100
    if ratio_indels <= 0.26:  # ##(0.26 = 5*1.98/ 38)  Indels below 5%
        score+=1

    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
        score+=1
    
    #print('Current score (TS0500) =', score,'\n')
    
    #Commento basato sul valore finale dello score
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return score


'''
###### DICTOPLOT per aggiungere bar plots e hopefully fisher
#funzione che mostra tramite bar plot la frequenza di mutazione 
def dictoplot(file):
    #formazione del dizionario vuoto
    dictionary = dizionario(file)
    #indels da aggiungere al plot
    lista_indels=[]
    lista_inserzioni=[]
    lista_delezioni=[]
    vcf_reader = vcf.Reader(open(file, 'r')) #leggi il vcf di input
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if len(record.ALT[alt])!=len(record.REF):
                lista_indels.append(record.REF)
                if len(record.ALT[alt])>len(record.REF):
                    lista_inserzioni.append(record.ALT[alt])
                if len(record.ALT[alt])<len(record.REF):
                    lista_delezioni.append(record.REF)
    dictio = {"Ins": len(lista_inserzioni), "Dels": len(lista_delezioni)}
    dictionary.update(dictio)
    # creazione del dataset
    courses = list(dictionary.keys())
    values = list(dictionary.values()) 
    fig = plt.figure(figsize = (10, 5))
    # creazione del bar plot
    plt.bar(courses, values, color ='green', width = 0.4)
    plt.xlabel("Mutation events")
    plt.ylabel("Mutation Frequency")
    plt.title(str(file))
    plt.savefig('my_.png')
    plt.show()
 
    return "Bar Plot of the mutations shown."


'''

import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                    prog = 'Parser e score per dati di VCF', 
                    description= 'Codice che genera un dizionario che conta i tipi di mutazioni presenti in un VCF filtrato, per poi generare uno score in base ai risultati presenti nel dictionary e riportare varie voci di output in base al risultato dello score.', 
                    epilog = 'Lo "score" è un valore numerico calcolato dalla script e sta a indicare la presenza di specifche condizioni richieste come input. Nella analisi del VCF, se una specifica condizione è rispettata, lo score aumenterà di un punto e il numero finale sarà dato dalla somma di valori rispettati',
                    add_help = True, )
        # arguments
    parser.add_argument('-f', '--folder', required=True, help='Pathway per la cartella con il file VCF')
    #parser.add_argument('-f2', '--folder2', required=False, help='Pathway per la cartella con un possibile secondo file VCF. utilizzabile per comparazione tramite funzione Dictoplot')
    parser.add_argument('-t', '--TMB', required=True, help='Tumour Mutational Burden (TMB)')
    #parser.add_argument('-m', '--MSI', required=True, help='Microsatellite instability (MSI)')
    #parser.add_argument('-r', '--mutazioniricorrenti', type= list, required=False, help='lista presentante posizioni in cui sono ricorrenti precisi eventi di ricombinazione')
    #parser.add_argument('-i', '--indels', type= list, required=False, help='lista di inserzioni e delezioni osservate')
    #parser.add_argument('-s', '--analizedscore', type= int, required=False, help='score che rappresenta la presenza o assenza della caratterstiche da noi richieste, analizzando i dati del file VCF di input')
    parser.add_argument('-c', '--Category', required=False, help='Category of mutations : hotspot,wt,exo,vus')
    parser.add_argument('-o', '--output', required=False, help='Path of file to store results')
    args = parser.parse_args()

    folder = args.folder
    #folder2 = args.folder2
    TMB = args.TMB
    #MSI = args.MSI
    #recurrent=args.mutazioniricorrenti
    #indels=args.indels
    #analizedscore=args.analizedscore
    category=args.Category
    output=args.output
    
    print("\n")
    print("VCF analizzato: ", folder, "\n")
    # print ("HEADER: ", printheader, "\n")
    # print('--------------------------', "\n")
    
    # print('Numero mutazioni osservate nel VCF =',"\n", dizionario(folder),'\n')
    
    # print('Mutazioni totali osservate = ', eventimutazionetotali(folder),'\n')
    
    # print(valorepercentuale(folder), '\n')

    
    print('--------------------------', "\n")
    
    # print('Eventi di inserzione e delezione trovati = ', listaindels(folder),'\n')
    # print("Numero totale di mutazioni indels = ", len(listaindels(folder)), '\n')
    
    
    # print('##########################', "\n")
    
    # print('Recurrent Mutations = ', recurrentmutations(folder),'\n')
    
    # print('##########################', "\n")    
    # print('##########################', "\n") 
    
    
    score_castillo,ca_perc, tg_perc,cg_perc, ratio_indels=calcoloscore_CASTILLO(folder, TMB)
    
    print('Commento sullo score (CASTILLO) = ',score_castillo ,'\n')
    
    print('##########################', "\n") 
    score_tso=calcoloscore_TS0500(folder, TMB)
    print('Commento sullo score (TS0500) = ', calcoloscore_TS0500(folder, TMB),'\n')
    #print('Commento sullo score = ', calcoloscore(folder, TMB, MSI),'\n')
    
    
    #print(dictoplot(folder), '\n') 
'''
    if os.path.exists(output):
        file=open(output,"a")
    else:
        file=open(output,"w")
        file.write("\t".join(["SampleID","TMB","C>A(%)","C>G(%)","T>G(%)","Indels","Score_Castillo","Score_Prop","Category"])+"\n")
    
    sample_id=folder.split("/")[-1]
    file.write("\t".join([str(sample_id),str(TMB),str(ca_perc), str(cg_perc),str(tg_perc), str(ratio_indels),str(score_castillo),str(score_tso),str(category)])+"\n")
'''