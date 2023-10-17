import vcf
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 

def dicto(file):

    mutations= ('A>C', 'A>G', 'A>T', 'C>A', 'C>G', 'C>T', 'G>A', 'G>C', 'G>T', 'T>A', 'T>C', 'T>G')
    dictionary=dict.fromkeys(mutations)
    number=0
    dictionary=dict.fromkeys(mutations,number)
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
        
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



def recurrentmutations(file):
    vcf_reader = vcf.Reader(open(file, 'r'))
    totalmuts={}
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if record.CHROM == 'chr12':
                x=record.POS
                y= str(record.REF)+">"+str(record.ALT[alt])
                totalmuts[x] = y

    listrecurrent = {133253184: 'G>C', 133250289: 'C>A', 133253151: 'G>A', 133249847: 'G>A', 133249857: 'C>G', 133252327: 'A>G', 133250250: 'G>T', 133253157: 'A>C', 133250213: 'G>C', 133250189: 'A>T', 133252325: 'C>A', 133250250: 'G>C', 133250238: 'C>T', 133250250: 'G>C', 133253208: 'G>A', 133249829: 'G>A', 133256623: 'G>A', 133252023: 'T>G', 133257828: 'G>A'}

    
    shared_items = {k: listrecurrent[k] for k in listrecurrent if k in totalmuts and listrecurrent[k] == totalmuts[k]}
    if len(shared_items) == 0:
        return shared_items   
    else:
        return shared_items


def list_indels(file):
    indels_list=[]
    lista_inserzioni=[]
    lista_delezioni=[]
    vcf_reader = vcf.Reader(open(file, 'r')) #leggi il vcf di input
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if len(record.ALT[alt])!=len(record.REF):
                indels_list.append(record.REF)
                if len(record.ALT[alt])>len(record.REF):
                    lista_inserzioni.append(record.ALT[alt])
                if len(record.ALT[alt])<len(record.REF):
                    lista_delezioni.append(record.REF)
    
    return indels_list


def totalmutationevents(file):
    mutation_totals=0   
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
            mutation_totals+=1

    return mutation_totals



def mutationsfrequency(file):
    dict = dicto(file)
    totalsos = sum(dict.values())
    listpercent = []
    for key, value in dict.items():
        percent = round((int(value) / totalsos),2)*100
        listpercent.append(percent)
    return ' '




def polescore(file, TMB):
    score=0
    dictionary = dicto(file)
    indels_list= list_indels(file)
    mutation_totals = totalmutationevents(file)
    totalsos = sum(dictionary.values())

    if float(TMB)>= 100: 
        score +=1
    for key, value in dictionary.items():

        if key=="C>A":
            if round((int(value) / totalsos),4)*100 >= 6:
                score +=1
               
        if key=="T>G": 
            if round((int(value) / totalsos),4)*100 >= 4:
                score +=1
            
        if key=="C>G":
            if round((int(value) / totalsos),4)*100 <= 5:
                score +=1

    ratio_indels=round((len(indels_list)/ mutation_totals),4) * 100
    if ratio_indels <= 4:
        score+=1

    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if len(recurrent) != 0:
        response = "+    Recurrent Mutations in EC"
        score+=1
    
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return score


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
    parser.add_argument('-o', '--output', required=False, help='Path of file to store results')
    
    args = parser.parse_args()

    folder = args.folder
    TMB = args.TMB
    output=args.output
    
    print("\n")
    print("Analyzed VCF: ", folder, "\n")

    
    print('--------------------------', "\n")
    
    
    print('##########################', "\n") 
    score_tso=polescore(folder, TMB)
    print('Score Pathogenicity = ', polescore(folder, TMB),'\n')







