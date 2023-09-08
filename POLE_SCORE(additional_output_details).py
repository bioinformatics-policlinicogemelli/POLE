#The following script generates a dictionary depending on the mutation frequency in an input VCF
#From these results, a score is generated to report various outputs depending on the vallue of the score itself, in order to describe the level of pethogenicity of the mutation.
#Based on "Interpretation of somatic POLE mutations in endometrial carcinoma" (Castillo et. all, 2020)


import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#This function generates a dictionary with the type of mutations as keys and the amount of times the specific mutations occurs as values
def dicto(file):
   
#This function generates a dictionary with the type of mutations as keys and the amount of times the specific mutations occurs as values
   
    #generation of an initial empty dictionary
    mutations= ('A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G')
    dictionary=dict.fromkeys(mutations)
    number=0
    dictionary=dict.fromkeys(mutations,number)
    
    # the following dictionary will be {'A>C': 0, 'A>G': 0, 'A>T': 0, 'C>A': 0, 'C>G': 0, 'C>T': 0, 'G>A': 0, 'G>C': 0, 'G>T': 0, 'T>A': 0, 'T>C': 0, 'T>G': 0}
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
           
        # (Mutations with both mutated alleles do not add more point to the Score)
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



#the function generates a dictionary with recurrent mutations as keys, along with the type of mutation in that position
def recurrentmutations(file):
    #positions=[]
    #mutations=[]
    vcf_reader = vcf.Reader(open(file, 'r'))
    totalmuts={} #dictionary with POS as keys and REF>ALT as values
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if record.CHROM == 'chr12':
                x=record.POS
                y= str(record.REF)+">"+str(record.ALT[alt])
                totalmuts[x] = y
    
    #Below the list of "Recurrent" mutations. Keys are the positions on the chromosome obtained by Varsome, values are the specific mutation events defined as recurrent in that position
    listrecurrent_extended = {133253184: 'G>C', 133250289: 'C>A', 133250208: 'C>A', 133253151: 'G>A', 133249349: 'G>A', 133252729: 'G>A', 133249847: 'G>A', 133249766: 'G>A', 133249857: 'C>G', 133252327: 'A>G', 133250250: 'G>T', 133250169: 'G>T', 133253157: 'A>C', 133249355: 'A>C', 133225944: 'A>C', 133252735: 'A>C', 133250213: 'G>C', 133249835: 'G>C', 133250189: 'A>T', 133249811: 'A>T', 133252325: 'C>A', 133248833: 'C>A', 133252027: 'C>A', 133245002: 'G>A', 133218351: 'T>C', 133235946: 'C>A', 133250250: 'G>C', 133250238: 'C>T', 133244183: 'C>T', 133225894: 'G>A', 133253208: 'G>A', 133249829: 'G>A', 133249841: 'G>A', 133256623: 'G>A', 133237646: 'A>C', 133215791: 'C>A', 133252023: 'T>C', 133233976: 'C>T', 133214612: 'T>C', 133242015: 'C>A', 133257828: 'G>A', 133250250: 'G>C', 133253208: 'G>A', 133249829: 'G>A', 133256623: 'G>A', 133252023: 'T>C', 133257828: 'G>A'}

    shared_items = {k: listrecurrent_extended[k] for k in listrecurrent_extended if k in totalmuts and listrecurrent_extended[k] == totalmuts[k]}
    if len(shared_items) == 0:
        print("     No Recurrent Mutations found", "\n") 
        return shared_items   
    else: 
        return shared_items



#Function which returns the events of insertions and deletions
def list_indels(file):
    indels_list=[]
    ins_list=[]
    del_list=[]
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if len(record.ALT[alt])!=len(record.REF):
                indels_list.append(record.REF)
                if len(record.ALT[alt])>len(record.REF):
                    ins_list.append(record.ALT[alt])
                if len(record.ALT[alt])<len(record.REF):
                    del_list.append(record.REF)
    
    return indels_list



#total number of mutation events observed in the filtered VCF, both indels and non-indels.
def totalmutationevents(file):
    mutation_totals=0  
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
            mutation_totals+=1

    return mutation_totals


#Estimation of the mutation frequencies
def mutationsfrequency(file):
    dict = dicto(file)
    totalsos = sum(dict.values())
    listpercent = []
    print('Mutations observed compaared to the toal amount: ')
    for key, value in dict.items():
        percent = round((int(value) / totalsos),2)*100
        print("Mutation frequency:", key, "=", percent, "%")
        listpercent.append(percent)
    return ' '



#Score calculated based on the dictionary, list of indels and the number of total mutation events
def polescore(file, TMB):
  
    score=0
    
    print("POLE Score", "\n")
    print("Each one of the parameters which fall into the pathogenic threshold adds +1 to the Score.")
    print("(Lack of Recurrent Mutations does not raise the Score)", "\n")
    dictionary = dicto(file)
    indels_list= list_indels(file)
    mutation_totals = totalmutationevents(file)
    totalsos = sum(dictionary.values())
    #totalsos is the total count of observed substitution events (indels are excluded)

    if float(TMB)>= 100:
        print ("+    Tumour Mutational Burden higher than 100 mut/Mb. TMB =", TMB, "\n")
        score +=1
    for key, value in dictionary.items():

        if key=="C>A": #C>A over 6%

            if round((int(value) / totalsos),4)*100 >= 6: 
                print('+    C>A Mutations equal or higher than 6%)')
                score +=1
                
            print("                    Amount of C>A: ", dictionary["C>A"])
            print("                    Frequency C>A compared to all SNV: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>A compared to total mutations: ", round((int(value) / mutation_totals),3)*100, " %", "\n")
            print("\n")   
               
        if key=="T>G":

            if round((int(value) / totalsos),4)*100 >= 4:
                print('+    T>G mutations equal or higher than 4%')
                score +=1
                
            print("                    Amount of T>G: ", dictionary["T>G"])
            print("                    Frequency T>G compared to all SNV: ", round((int(value) / totalsos),3)*100, " %")
            print("                    Frequency T>G compared to total mutations: ", round((int(value) / mutation_totals),3)*100, " %", "\n")
            print("\n")  
            
        if key=="C>G": #C>G below 5%
            if round((int(value) / totalsos),4)*100 <= 5:
                print('+    C>G mutations equal or smaller than 5%')
                score +=1
                
            print("                    Amount of C>G: ", dictionary["C>G"])
            print("                    Frequency C>G compared to all SNV: ", round((int(value) / totalsos),4)*100, " %")
            print("                    Frequency C>G compared to total mutations: ", round((int(value) / mutation_totals), 3)*100, " %", "\n")
            print("\n") 
            
    ratio_indels=round((len(indels_list)/ mutation_totals),4) * 100
    #This is the Indels Frequency compared to the total amount of mutations
    
    
    score += 1 
    if ratio_indels <= 4:
        print('+    Amount of indels below 4%','\n')
        score-=1
        
    print ("        ###################################################")
    print ("                    Number of indels =", len(indels_list))    
    print ("                    Ratio indels =", ratio_indels, "%")
    print ("        ###################################################")
    print ("\n")

    
    response= "Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
        print(response, "\n")
        score+=1
    #Recurrent mutations are identified if this dictionary is not empty, leading to a +1 to the Score

    
    print('--------------------------', "\n")
    
    print('Current score =', score,'\n')
    
    #Output string forthe final value of the score
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
                    prog = 'VCF scoring system for POLE pathogenicity', 
                    description= 'Script which generates a set of objects that takes in account the mutation events in an input filtered VCF, to thene generate a score which estimates the patogenic nature of the POLE mutation.', 
                    epilog = 'The generated score is a numeric integer which is derived by the count of specific conditions correlated to the pathogenicity of POLE mutations: during a vcf analysis, if a pathogenicity threshold is met, the score will go up by +1. The final score is given by the sum of all the respected thresholds',
                    add_help = True, )
        # arguments
    parser.add_argument('-f', '--folder', required=True, help='Path to the VCF folder')
    parser.add_argument('-t', '--TMB', required=True, help='Tumour Mutational Burden (TMB)')
    parser.add_argument('-i', '--indels', type= list, required=False, help='List of Indels')
    parser.add_argument('-s', '--analizedscore', type= int, required=False, help='score che rappresenta la presenza o assenza della caratterstiche da noi richieste, analizzando i dati del file VCF di input')
    parser.add_argument('-c', '--Category', required=False, help='Category of mutations : hotspot,wt,exo,vus')
    parser.add_argument('-o', '--output', required=False, help='Path of file to store results')

    args = parser.parse_args()

    folder = args.folder
    TMB = args.TMB
    indels=args.indels
    analizedscore=args.analizedscore
    category=args.Category
    output=args.output

    print("\n")
    print("VCF: ", folder, "\n")
    print('--------------------------', "\n")
    
    print('VCF observed Mutations =',"\n", dicto(folder),'\n')
    
    print('Total number of observed mutations = ', totalmutationevents(folder),'\n')
    
    print(mutationsfrequency(folder), '\n')

    
    print('--------------------------', "\n")
    
    print('Indels = ', list_indels(folder),'\n')
    print("Total number of Indels = ", len(list_indels(folder)), '\n')
    
    
    print('--------------------------', "\n")  
    
    
    print('Recurrent Mutations = ', recurrentmutations(folder),'\n')
    
    
    print('--------------------------', "\n")
    

    
    print('Score Pathogenicity = ', polescore(folder, TMB),'\n')
