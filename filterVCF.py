import argparse
import pandas as pd
import vcf
import os
from configparser import ConfigParser
from loguru import logger


def parsing_vcf(file_input,file_output):
    #lista contenente i dati da #CHROM in poi (no meta-information)
    with open (file_input) as f:
        correct_data=[line.strip().split("\t") for line in f if not line.startswith("##")]

    #creazione dataframe contenente i dati
    data=pd.DataFrame(correct_data[1:],columns=correct_data[0])
    
    #filtro per dati in cui FILTER == "PASS"
    data=data[data["FILTER"]=="PASS"]
    
    #filtro per dati in cui ALT!="."
    data=data[data["ALT"]!="."]
    
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
    
    

def main_filter(input, output):

    
    if not os.path.exists(output):
        logger.warning("Creating folder to store filtered vcf")
        os.makedirs(output)
    elif not os.path.isdir(output):
        print("Output is not a folder")
        exit(1)
    
    metadata = read_vcf_meta(input)
  
    filtered_path=os.path.join(output,input.split("/")[-1]+"_filtered.vcf").replace(".vcf_","_")    
    
    write_vcf_meta(filtered_path, metadata)
    parsing_vcf(input,filtered_path)






