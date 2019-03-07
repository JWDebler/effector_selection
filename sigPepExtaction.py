#!/usr/bin/env python3

import os
import re
import csv
import argparse
import glob
from pathlib import Path

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a folder containing another folder with input files (default is a folder called "input" inside the directory where the script is)')
parser.add_argument('-o', '--output', help='folder for output files (default "output" inside the directory where the script is)')
args = parser.parse_args()

if args.input:
    input_path = Path(args.input)
else:
    input_path = Path(os.path.join(os.getcwd(),'input'))

if args.output:
    output_path = Path(args.output)
else:
    output_path = Path(os.path.join(os.getcwd(), 'output'))


print("Reading content of input directory")
print('============================================')
samples = [dI for dI in os.listdir(input_path) if os.path.isdir(os.path.join(input_path,dI))]
print("found following isolates:")
for sample in samples:
    print(sample)
for sample in samples:

     #location of input files
    print('============================================')
    print(sample,"- Checking if all required files are present")

    input_dir = input_path / sample
    
    if os.listdir(input_dir):
        
        protein_file = '*proteins.fasta'
        protein_file_name = os.path.join(input_dir, protein_file)
        if glob.glob(protein_file_name) == []:
            print('ERROR - missin proteins.fasta')
        else:
            for element in glob.glob(protein_file_name):
                protein_file_path = Path(element)  
                protein_file_path = Path(element)

        interproscan_file = '*interproscan.tsv'
        interproscan_file_name = os.path.join(input_dir, interproscan_file)
        if glob.glob(interproscan_file_name) == []:
            print('ERROR - missing interproscan.tsv')
            break
        else:
            for element in glob.glob(interproscan_file_name):
                interproscan_file_path = Path(element)
                interproscan_file_path = Path(element)

        deepsig_file = '*deepsig.out'
        deepsig_file_name = os.path.join(input_dir, deepsig_file)

        if glob.glob(deepsig_file_name) == []:
            print('ERROR - missing deepsig.out')
            break
        else:
            for element in glob.glob(deepsig_file_name):
                deepsig_file_path = Path(element)
                deepsig_file_path = Path(element)

    else:
        print("ERROR - no input files for", sample)
        break

    candidate_output_dir = output_path / sample 

    try:
        os.makedirs(candidate_output_dir)
    except Exception:
        pass

    #delete content of folder to avoid appending files
    sigPeptideFasta = sample + ".proteinsWithSig.fasta"
    for file in os.listdir(candidate_output_dir):
        if file == sigPeptideFasta:
            print('INFO: outputfile already present, overwriting it')
            pathToSigPeptideFasta = os.path.join(candidate_output_dir, file)
            try:
                os.unlink(pathToSigPeptideFasta)
            except Exception as e:
                print(e)
        
    proteins = {}
    deepsig = {}
    interproScan = {}
    signalP = {}
    signalPandDeepsig = {}

    #Read fasta file with protein sequences and IDs
    #And fill lists with default values
    with open(protein_file_path) as file:
        print(sample,'- reading proteins file')
        input = file.read().splitlines()
        for line in input:
            if not line.strip(): continue
            if line[0] == '>':
                name = line[1:]
                proteins[name] = ''
                signalP[name] = 0
                deepsig[name] = 0
            else:
                proteins[name] += line

    print(sample,'- processed proteins data')

    #Read deepsig output
    with open(deepsig_file_path) as file:
        print(sample,'- reading deepsig file')
        input = csv.reader(file, delimiter='\t')
        for line in input:
            if line[1] == 'SignalPeptide':
                deepsig[line[0]] = int(line[3])

    print(sample,'- processed deepsig data')

    #Read interproScan output
    with open(interproscan_file_path) as file:
        print(sample,'- reading interproScan file')
        input = csv.reader(file, delimiter='\t')
        for line in input:
            id = line[0]
            if id in interproScan:
                interproScan[id].append(line[3:])
                if line[3] == 'SignalP_EUK':
                    signalP[id] = int(line[7])
            else:
                interproScan[id] = [line[3:]]
                if line[3] == 'SignalP_EUK':
                    signalP[id] = int(line[7])

    print(sample, '- processed interproScan data')

    # Combine predictions where both SignalP and deepsig predict a signal peptide,
    # but if they can't agree on the cleavage site prefer SignalP over deepsig
    for key, value in signalP.items():
        if value > 0 and deepsig[key] > 0:
            if value == deepsig[key]:
                signalPandDeepsig[key] = value
            else:
                signalPandDeepsig[key] = value

    print('============================================')
    print("Extracting proteins with signal Peptide")
    print("There are", len(signalPandDeepsig),"proteins with a SP")

    output_name_fasta = sample+'.proteinsWithSig.fasta'
    output_file_fasta = candidate_output_dir / output_name_fasta

    for element in signalPandDeepsig:
        print('>'+element, '\n'+proteins[element], file=open(os.path.join(candidate_output_dir, output_name_fasta),'a'))

    print("All Done, check:", output_file_fasta)

