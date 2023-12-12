######################################
# NAME: pathopoleanalyzer.py
# version ="1.0"
# ====================================

import os
import sys
import vcf
import numpy as np
import pandas as pd
import argparse
from loguru import logger
from configparser import ConfigParser
from filterVCF import main_filter

config = ConfigParser()
configFile = config.read("conf.ini")



def dicto(file):
    mutations= ('C>A', 'C>G','T>G')
    dictionary=dict.fromkeys(mutations)
    number=0
    dictionary=dict.fromkeys(mutations,number)
    vcf_reader = vcf.Reader(open(file, 'r'))
    for record in vcf_reader:
        for alt in range(0,len(record.ALT)):
            if record.REF == 'C' and record.ALT[alt]=='A':
                dictionary['C>A']+=1
            if record.REF == 'C' and record.ALT[alt]=='G':
                dictionary['C>G']+=1
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

    listrecurrent = {133253184: 'G>C',
                     133250289: 'C>A', 
                     133253151: 'G>A', 
                     133249847: 'G>A',
                     133249857: 'C>G',
                     133252327: 'A>G', 
                     133250250: 'G>T',
                     133253157: 'A>C',
                     133250213: 'G>C',
                     133250189: 'A>T', 
                     133252325: 'C>A', 
                     133250250: 'G>C', 
                     133250238: 'C>T', 
                     133250250: 'G>C',
                     133253208: 'G>A',
                     133249829: 'G>A', 
                     133256623: 'G>A',
                     133252023: 'T>G',
                     133257828: 'G>A'}
    
    shared_items = {k: listrecurrent[k] for k in listrecurrent if k in totalmuts and listrecurrent[k] == totalmuts[k]}
    if len(shared_items) == 0:
        return "None"   
    else:
        return True


def list_indels(file):
    indels_list=[]
    lista_inserzioni=[]
    lista_delezioni=[]
    vcf_reader = vcf.Reader(open(file, 'r')) 
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
    if float(TMB)>= float(config.get('Threshold', 'TMB')): 
        score +=1
    for key, value in dictionary.items():

        if key=="C>A":
            if round((int(value) / totalsos),4)*100 >= float(config.get('Threshold', 'C_A')):
                score +=1
               
        if key=="T>G": 
            if round((int(value) / totalsos),4)*100 >= float(config.get('Threshold', 'T_G')):
                score +=1
            
        if key=="C>G":
            if round((int(value) / totalsos),4)*100 <= float(config.get('Threshold', 'C_G')):
                score +=1

    
    c_a=round(dictionary["C>A"]/totalsos,2)*100
    c_g=round(dictionary["C>G"]/totalsos,2)*100
    t_g=round(dictionary["T>G"]/totalsos,2)*100

    ratio_indels=round((len(indels_list)/ mutation_totals),4) * 100
    if ratio_indels <= float(config.get('Threshold', 'Indels')):
        score+=1

    response="Recurrent Mutations not found"
    recurrent = recurrentmutations(file)    
    if recurrent != "None":
        response = "+    Recurrent Mutations in EC"
        score+=1
    if score >= 4:
        output = 'Pathogenic POLE mutation'
    if score == 3:
        output = 'Variant of unknown significance'
    if score < 3:
        output= 'Non-Pathogenic POLE mutation'
    return score,TMB,recurrent,round(ratio_indels,2),round(c_a,2),round(c_g,2),round(t_g,2)



if __name__ == '__main__':
    
    logfile="PathoPOLE_{time:YYYY-MM-DD_HH-mm-ss.SS}.log"
    logger.level("INFO", color="<green>")
   # logger.add(sys.stderr, format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",colorize=True, catch=True, backtrace=True, diagnose=True)
    logger.add(os.path.join('Logs',logfile),format="{time:YYYY-MM-DD_HH-mm-ss.SS} | <lvl>{level} </lvl>| {message}",mode="w", backtrace=True, diagnose=True)
    logger.info("Welcome to Patho POLE Analyzer")

    parser = argparse.ArgumentParser(
                    prog = 'VCF scoring system for POLE pathogenicity', 
                    description= 'Script which generates a set of objects that takes in account the mutation events in an input filtered VCF, to thene generate a score which estimates the patogenic nature of the POLE mutation.', 
                    epilog = 'The generated score is a numeric integer which is derived by the count of specific conditions correlated to the pathogenicity of POLE mutations: during a vcf analysis, if a pathogenicity threshold is met, the score will go up by +1. The final score is given by the sum of all the respected thresholds',
                    add_help = True, )
    parser.add_argument('-i', '--Input', required=True, help='Path to the VCF folder')
    parser.add_argument('-t', '--TMB', required=True, help='Tumour Mutational Burden (TMB) or file')
    parser.add_argument('-o', '--output', required=False, help='Path of file to store results')
    parser.add_argument('-w', '--overWrite', required=False, action='store_true',
                        help='Overwrite result file')
    
    try:
        args = parser.parse_args()
    except ValueError:
        logger.critical("Error Argument: missing required input")
        exit(1)
    

    filter_output=config.get("Filter","path_filtered_vcf")
    
    input=args.Input

    if args.output:
            output=args.output
    else:
        logger.warning("Using default output file 'Results_pathopoleanalyzer.txt'")
        output=config.get("Results","results_file")
        


    if os.path.isfile(input):

        if not input.endswith(".vcf"):
            logger.critical("Input file is not a VCF")
            exit(1)

        main_filter(input,filter_output)
        TMB = args.TMB
        
      
        logger.info("Analyzing VCF "+input)
        score_tso,TMB,recurrent,ratio_indels,c_a,c_g,t_g=polescore(os.path.join(filter_output,input.split("/")[-1]+"_filtered.vcf").replace(".vcf_","_"),TMB)
        logger.info('Score Pathogenicity = +'+str(score_tso))
    
        if config.get("Results","store_results"): 
            if os.path.exists(output):
                file=open(output,"a")
            else:
                file=open(output,"w")
                file.write("\t".join(["SampleID","TMB","RecurrentMutation","Indels","C>A","C>G","T>G","Score"])+"\n")
            
            sample_id=input.split("/")[-1]

            file.write("\t".join([str(sample_id),str(TMB).strip(),str(recurrent), str(ratio_indels),str(c_a),str(c_g),str(t_g),str(score_tso)])+"\n")

        logger.success("VCF analysis completed successfully!")

    elif os.path.isdir(args.Input):
        logger.warning("Input is a folder. Reading multiple files")
        
        files=[file for file in os.listdir(args.Input) if ".vcf" in file]
        if len(files)==0:
            logger.critical("No vcf found in input folder")
            exit(1)

        logger.info("Found "+str(len(files))+" vcf files")

        for file in files:
            main_filter(os.path.join(input,file),filter_output)
             

        if not os.path.isfile(args.TMB):
            logger.critical("Error: TMB must be a file")
            exit()
        else:
            tmb=pd.read_csv(args.TMB,sep="\t")
            
            results=pd.DataFrame(columns=["SampleID","TMB","RecurrentMutation","Indels","C>A","C>G","T>G","Score"])
            for file in files:

                tmb_sample=tmb.loc[tmb["Sample"]==file,"TMB"].values[0]
           
               
                logger.info("Analyzing VCF "+file)
                score_tso,tmb_sample,recurrent,ratio_indels,c_a,c_g,t_g=polescore(os.path.join(filter_output,file.split("/")[-1]+"_filtered.vcf").replace(".vcf_","_"),tmb_sample)
                logger.info('Score Pathogenicity of '+file +': '+str(score_tso))

                sample_id=file.split("/")[-1]
                temp=pd.DataFrame({"SampleID":str(sample_id),"TMB":tmb_sample,"RecurrentMutation":recurrent,"Indels":ratio_indels,
                                   "C>A":c_a,"C>G":c_g,"T>G":t_g,"Score":score_tso},index=[0])


                results=pd.concat([results,temp],axis=0).reset_index(drop=True)

            if config.get("Results","store_results"): 
                logger.info("Writing results in output file")
                if os.path.exists(output) and not args.overWrite:
                   results.to_csv(output, sep="\t", mode="a",header=False,index=False)
                else:
                    results.to_csv(output, sep="\t", mode="w",index=False)
                    
            logger.success("VCF analysis completed successfully!")