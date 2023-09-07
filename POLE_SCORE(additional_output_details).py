#The following script generates a dictionari depending on the mutation frequency in an input VCF
#Dopo questo genera uno score in base ai risultati presenti nel dictionary e ti riporta varie voci di output in base al risultato dello score.
#le richieste di output presenti in questo script sono basate sull'articolo "Interpretation of somatic POLE mutations in endometrial carcinoma" (Castillo et. all, 2020)


import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#this function generates a dictionary with the type of mutations as keys and the amount of times the specific mutations occurs as values dicto che in cui le mutazioni sono le keys e i valori sono la Frequency in cui queste mutazioni appaiono nel VCF
def dicto(file):
   
#this function generates a dictionary with the type of mutations as keys and the amount of times the specific mutations occurs as values
   
    #generation of an initial empty dictionary
    mutations= ('A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G')
    dictionary=dict.fromkeys(mutations)
    number=0
    dictionary=dict.fromkeys(mutations,number)
    
    # the following dictionary will be {'A>C': 0, 'A>G': 0, 'A>T': 0, 'C>A': 0, 'C>G': 0, 'C>T': 0, 'G>A': 0, 'G>C': 0, 'G>T': 0, 'T>A': 0, 'T>C': 0, 'T>G': 0}
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



#funzione che crea un dicto con posizioni relevant come keys e il tipo di mutazione in quella posizione, dati basati sulla tabella del analisi di POLE nel articolo: "Interpretation of somatic POLE mutations in endometrial carcinoma"
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
    
    #Below the list of "Recurrent" mutations. Le keys sono le posizioni trovate tramite sito Varsome per le rispettive mutazioni, e i values sono gli eventi di mutazione chiave
    #se tali valori coincidono con quelli visti nel VCF filtrato, allora ci sono "current mutations in EC", condizione che aumenta di +1 lo score generale
    listrecurrent_extended = {133253184: 'G>C', 133250289: 'C>A', 133250208: 'C>A', 133253151: 'G>A', 133249349: 'G>A', 133252729: 'G>A', 133249847: 'G>A', 133249766: 'G>A', 133249857: 'C>G', 133252327: 'A>G', 133250250: 'G>T', 133250169: 'G>T', 133253157: 'A>C', 133249355: 'A>C', 133225944: 'A>C', 133252735: 'A>C', 133250213: 'G>C', 133249835: 'G>C', 133250189: 'A>T', 133249811: 'A>T', 133252325: 'C>A', 133248833: 'C>A', 133252027: 'C>A', 133245002: 'G>A', 133218351: 'T>C', 133235946: 'C>A', 133250250: 'G>C', 133250238: 'C>T', 133244183: 'C>T', 133225894: 'G>A', 133253208: 'G>A', 133249829: 'G>A', 133249841: 'G>A', 133256623: 'G>A', 133237646: 'A>C', 133215791: 'C>A', 133252023: 'T>C', 133233976: 'C>T', 133214612: 'T>C', 133242015: 'C>A', 133257828: 'G>A', 133250250: 'G>C', 133253208: 'G>A', 133249829: 'G>A', 133256623: 'G>A', 133252023: 'T>C', 133257828: 'G>A'}

    #print("Mut ricorrenti = ", listrecurrent)
    shared_items = {k: listrecurrent_extended[k] for k in listrecurrent_extended if k in totalmuts and listrecurrent_extended[k] == totalmuts[k]}
    if len(shared_items) == 0:
        print("     No Recurrent Mutations found", "\n") 
        return shared_items   
    else: 
        return shared_items



#funzione che ritorna la lista di eventi di inserzione e delezione
def list_indels(file):
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
def totalmutationevents(file):
    totalemutazioni=0
    #totalemutazioni è il numero totale di elementi nel vcf filtrato, siano essi indels o non indels.    
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
            totalemutazioni+=1

    return totalemutazioni


#Estimation of the mutation frequencies
def mutationsfrequency(file):
    dict = dicto(file)
    totalsos = sum(dict.values())
    #totalsos è la conta totale degli eventi di sostituzione osservati (eventi come gli indels sono scartati da questa conta)
    listpercent = []
    print('PERCENTUALE DELLE MUTAZIONI OSSERVATE RISPETTO AL TOTALE: ')
    for key, value in dict.items():
        percent = round((int(value) / totalsos),2)*100 #Frequency percentuale "classico" 
        #percent2 = percent*1.98/ 38 # Frequency percentuale modificata in relazione alle proporzioni usate da Castillo; 38 : Frequency_castillo = 1,98 : Frequency_TSO500
    #########    if user_input.lower() == 'exo':
        print("Frequency in Percentuale di mutazione:", key, "=", percent, "%")
        listpercent.append(percent)
    #########    elif user_input.lower() == 'panel':
    #########        print("Frequency in Percentuale di mutazione:", key, "=", percent2, "%")
    #########        listpercent.append(percent2)
    #return listpercent
    return ' '



#CALCOLO SCORE IN BASE AI NOSTRI PARAMETRI PRESI DA DICTIONARY (LISTA_INDELS E NUMERO DI EVENTI TOTALI)
#def polescore(file, TMB):
def polescore(file, TMB):
  
    score=0
    
    print ("POLE SCORE", "\n")
    print("A seguire sono mostrate le condizioni rispettate, ognuna di esse aumenta lo Score di +1.")
    print("(Lack of Recurrent Mutations does not raise the Score)", "\n")
    dictionary = dicto(file)
    lista_indels= list_indels(file)
    totalemutazioni = totalmutationevents(file)
    totalsos = sum(dictionary.values())
    #totalsos anche qui è la conta totale degli observed substitution events (indels are excluded)


    print('Percentuali utilizzate per il calcolo dello score modificate rispetto al paper di Castillo et al.', '\n', 'Modifica in base alla proporzione: [38 : frequency_Castillo = 1,98 : frequency_TS0500]', '\n')
    if float(TMB)>= 5.21: #(5.21 = 100*1.98/ 38) #TMB over 100mut/Mb
        print ("+    Tumour Mutational Burden maggiore di 100 mut/Mb. TMB =", TMB, "\n")
        score +=1
    for key, value in dictionary.items():

        if key=="C>A": #C>A over 20%

            if round((int(value) / totalsos),4)*100 >= 1.0421: ## (1.0421 = 20*1.98/ 38)
                print('+    Mutazioni C>A uguali o superiori al valore richiesto (maggiore del 20% in Castillo)')
                score +=1
                
            print("                    Quantità di C>A: ", dictionary["C>A"])
            print("                    Frequency C>A: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>A compared to total mutations: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            print("\n")   
               
        if key=="T>G": #T>G over 4%

            if round((int(value) / totalsos),4)*100 >= 0.208: ## (0.208 = 4*1.98/ 38)
                print('+    Mutazioni T>G uguali o superiori al valore richiesto (maggiore del 4% in Castillo)')
                score +=1
                
            print("                    Quantità di T>G: ", dictionary["T>G"])
            print("                    Frequency T>G: ", round((int(value) / totalsos),3)*100, " %")
            print("                    Frequency T>G compared to total mutations: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            print("\n")  
            
        if key=="C>G": #C>G below 0.6%
            if round((int(value) / totalsos),4)*100 <= 0.0312: ## (0.0312 = 0.6*1.98/ 38)
                print('+    Mutazioni C>G uguali o inferiori al valore richiesto (minore del 0.6% in Castillo)')
                score +=1
                
            print("                    Quantità di C>G: ", dictionary["C>G"])
            print("                    Frequency C>G: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>G compared to total mutations: ", round((int(value) / totalemutazioni), 3)*100, " %", "\n")
            print("\n") 
            
    ratio_indels=round((len(lista_indels)/ totalemutazioni),4) * 100
    #This is the Indels Frequency compared to the total amount of mutations
    
    
    score += 1 
    if ratio_indels <= 0.26:  # ##(0.26 = 5*1.98/ 38)  Indels below 5%
        print('+    Quantità di indels uguale o inferiore al valore richiesto (minore del 5% in Castillo)','\n')
        score-=1
        
    print ("        ###################################################")
    print ("                    Numero indels =", len(lista_indels))    
    print ("                    Ratio indels =", ratio_indels, "%")
    print ("        ###################################################")
    print ("\n")

    
    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
        print(response, "\n")
        score+=1
    #se il dictionary recurrent non è vuoto, sono state identificate Recurrent mutations, perciò lo Score aumenta di +1

    
    print('--------------------------', "\n")
    
    print('Current score =', score,'\n')
    
    #Commento basato sul valore finale dello score
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return output


def polescore_CASTILLO(file, TMB):
  
    score=0
    
    print("CALCOLO DELLO SCORE", "\n")
    print("A seguire sono mostrate le condizioni rispettate, ognuna di esse aumenta lo Score di +1.")
    print("(La mancanza di Recurrent Mutations non porta a un aumento dello score)", "\n")
    dictionary = dicto(file)
    lista_indels= list_indels(file)
    totalemutazioni = totalmutationevents(file)
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
        print ("+    Tumour Mutational Burden maggiore di 100 mut/Mb. TMB =", TMB, "\n")
        score +=1
    for key, value in dictionary.items():

        if key=="C>A": #C>A over 20%

            if round((int(value) / totalsos),4)*100 >= 20: 
                print('+    Mutazioni C>A uguali o superiori al valore richiesto (maggiore del 20% in Castillo)')
                score +=1
                
            print("                    Quantità di C>A: ", dictionary["C>A"])
            print("                    Frequency C>A: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>A compared to total mutations: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            print("\n")   
            ca_perc=round((int(value) / totalsos),4)*100 
        if key=="T>G": #T>G over 4%

            if round((int(value) / totalsos),4)*100 >= 4:
                print('+    Mutazioni T>G uguali o superiori al valore richiesto (maggiore del 4% in Castillo)')
                score +=1
                
            print("                    Quantità di T>G: ", dictionary["T>G"])
            print("                    Frequency T>G: ", round((int(value) / totalsos),3)*100, " %")
            print("                    Frequency T>G compared to total mutations: ", round((int(value) / totalemutazioni),3)*100, " %", "\n")
            print("\n")  
            tg_perc=round((int(value) / totalsos),4)*100 
        if key=="C>G": #C>G below 0.6%
            if round((int(value) / totalsos),4)*100 <= 0.6: 
                print('+    Mutazioni C>G uguali o inferiori al valore richiesto (minore del 0.6% in Castillo)')
                score +=1
                
            print("                    Quantità di C>G: ", dictionary["C>G"])
            print("                    Frequency C>G: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>G compared to total mutations: ", round((int(value) / totalemutazioni), 3)*100, " %", "\n")
            print("\n") 
            cg_perc=round((int(value) / totalsos),4)*100 
            
    ratio_indels=round((len(lista_indels)/ totalemutazioni),4) * 100
    #Qui abbiamo la Frequency delle indels rispetto al totale
    
    

    if ratio_indels <= 5:  # ##(0.26 = 5*1.98/ 38)  Indels below 5%
        print('+    Quantità di indels uguale o inferiore al valore richiesto (minore del 5% in Castillo)','\n')
        score+=1
        
    print ("        ###################################################")
    print ("                    Numero indels =", len(lista_indels))    
    print ("                    Ratio indels =", ratio_indels, "%")
    # print ("        ###################################################")
    print ("\n")
    
    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
        print(response, "\n")
        score+=1
    #se il dictionary recurrent non è vuoto, sono state identificate Recurrent mutations, perciò lo Score aumenta di +1

    
    print('--------------------------', "\n")
    
    print('Current score (CASTILLO) =', score,'\n')
    
    #Commento basato sul valore finale dello score
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return output


#CALCOLO SCORE IN BASE AI NOSTRI PARAMETRI PRESI DA DICTIONARY (LISTA_INDELS E NUMERO DI EVENTI TOTALI)
#def polescore(file, TMB):
def polescore_TS0500(file, TMB):
    score=0
    dictionary = dicto(file)
    lista_indels= list_indels(file)
    totalemutazioni = totalmutationevents(file)
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
    
    print('Current score (TS0500) =', score,'\n')
    
    #Commento basato sul valore finale dello score
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return output


import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
                    prog = 'Parser e score per dati di VCF', 
                    description= 'Codice che genera un dicto che conta i tipi di mutazioni presenti in un VCF filtrato, per poi generare uno score in base ai risultati presenti nel dictionary e riportare varie voci di output in base al risultato dello score.', 
                    epilog = 'Lo "score" è un valore numerico calcolato dalla script e sta a indicare la presenza di specifche condizioni richieste come input. Nella analisi del VCF, se una specifica condizione è rispettata, lo score aumenterà di un punto e il numero finale sarà dato dalla somma di valori rispettati',
                    add_help = True, )
        # arguments
    parser.add_argument('-f', '--folder', required=True, help='Pathway per la cartella con il file VCF')
    parser.add_argument('-f2', '--folder2', required=False, help='Pathway per la cartella con un possibile secondo file VCF. utilizzabile per comparazione tramite funzione Dictoplot')
    parser.add_argument('-t', '--TMB', required=True, help='Tumour Mutational Burden (TMB)')
    #parser.add_argument('-m', '--MSI', required=True, help='Microsatellite instability (MSI)')
    #parser.add_argument('-r', '--mutazioniricorrenti', type= list, required=False, help='lista presentante posizioni in cui sono ricorrenti precisi eventi di ricombinazione')
    parser.add_argument('-i', '--indels', type= list, required=False, help='lista di inserzioni e delezioni osservate')
    parser.add_argument('-s', '--analizedscore', type= int, required=False, help='score che rappresenta la presenza o assenza della caratterstiche da noi richieste, analizzando i dati del file VCF di input')

    args = parser.parse_args()

    folder = args.folder
    folder2 = args.folder2
    TMB = args.TMB
    #MSI = args.MSI
    #recurrent=args.mutazioniricorrenti
    indels=args.indels
    analizedscore=args.analizedscore

    print("\n")
    print("VCF analizzato: ", folder, "\n")
    print ("HEADER: ", printheader, "\n")
    print('--------------------------', "\n")
    
    print('Numero mutazioni osservate nel VCF =',"\n", dicto(folder),'\n')
    
    print('Mutazioni totali osservate = ', totalmutationevents(folder),'\n')
    
    print(mutationsfrequency(folder), '\n')

    
    print('--------------------------', "\n")
    
    print('Eventi di inserzione e delezione trovati = ', list_indels(folder),'\n')
    print("Numero totale di mutazioni indels = ", len(list_indels(folder)), '\n')
    
    
    print('--------------------------', "\n")  
    
    
    print('Recurrent Mutations = ', recurrentmutations(folder),'\n')
    
    
    print('--------------------------', "\n")
    

    #print('Commento sullo score = ', polescore(folder, TMB, MSI),'\n')
    
    print('Commento sullo score (CASTILLO) = ',polescore_CASTILLO(folder, TMB),'\n')
    print('Commento sullo score_TS0500 = ', polescore_TS0500(folder, TMB),'\n')



