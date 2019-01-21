import os
import re
import csv
import argparse
from pathlib import Path

# This script expects a folder with the name of the sample containing following input files:

#  ├── 206-15
#      ├── 206-15.deepsig.out
#      ├── 206-15.effectorP.tsv
#      ├── 206-15.genemark.proteins.fasta
#      ├── 206-15.interproscan.tsv
#      ├── 206-15.cazymes.txt <-- this one is optional

################ DEFAULT VALUES ################
#Conditions for selecting effector candidates
effectorPscore = 0.8 # min EffectorP score
MWcutoff = 25  # molecular weight of mature (cleaved signal peptide) candidate
Cyscutoff = 2 # min number of cysteins
################################################

#parse commandline arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='path to a folder containing another folder with input files (default is a folder called "input" inside the directory where the script is)')
parser.add_argument('-o', '--output', help='folder for output files (default "output" inside the directory where the script is)')
parser.add_argument('-e','--effectorPscore', help='set minimum effectoP score (default = 0.8)', type=float)
parser.add_argument('-m', '--MWcutoff', help='set maximum molecular weight of mature candidate in kDa (default = 25) ', type=int)
parser.add_argument('-c', '--cysteins', help='set minimum number of cysteins in mature candidate (default = 2)', type=int)
args = parser.parse_args()

if args.input:
    input_path = Path(args.input)
else:
    input_path = Path(os.path.join(os.getcwd(),'input'))

if args.output:
    output_path = Path(args.output)
else:
    output_path = Path(os.path.join(os.getcwd(), 'output'))

if args.effectorPscore == 0:
    effectorPscore = 0.0
elif args.effectorPscore:
    effectorPscore = args.effectorPscore

if args.MWcutoff == 0:
    MWcutoff = 0
elif args.MWcutoff:
    MWcutoff = args.MWcutoff

if args.cysteins == 0:
    Cyscutoff = 0
elif args.cysteins:
    Cyscutoff = args.cysteins

print("Reading content of input directory")
print('============================================')
samples = [dI for dI in os.listdir(input_path) if os.path.isdir(os.path.join(input_path,dI))]
print("found following isolates:")
for sample in samples:
    print(sample)

print('============================================')
print('Candidate Selection criteria:')
print('EffectoP score: >=', effectorPscore)
print('Molecular weight: <=', MWcutoff)
print('Number of Cysteins: >=', Cyscutoff)
print('============================================')

#begin of actual analysis
for sample in samples:
    #folder where output files will be saved
    #folder name will contain selection criteria
    output = 'candidates.'+'effectorP='+str(effectorPscore)+'.MW='+str(MWcutoff)+'.Cysteins='+str(Cyscutoff)

    #location of input files
    print('============================================')
    print(sample,"- Checking if all required files are present")

    input_dir = input_path / sample
    protein_file = sample+'.genemark.proteins.fasta'
    protein_file_path = os.path.join(input_dir, protein_file)
    interproscan_file = sample+'.interproscan.tsv'
    interproscan_file_path = os.path.join(input_dir, interproscan_file)
    deepsig_file = sample+'.deepsig.out'
    deepsig_file_path = os.path.join(input_dir, deepsig_file)
    effectorP_file = sample+'.effectorP.tsv'
    effectorP_file_path = os.path.join(input_dir, effectorP_file)

    cazyme_file = sample+'.cazymes.txt'
    cazyme_file_path = os.path.join(input_dir, cazyme_file)
    cazyme_file_present = Path(cazyme_file_path).is_file()

    if not os.path.isfile(protein_file_path):
        print('ERROR - missing:', protein_file)
    if not os.path.isfile(interproscan_file_path):
        print('ERROR - missing:', interproscan_file)
    if not os.path.isfile(deepsig_file_path):
        print('ERROR - missing:', deepsig_file)
    if not os.path.isfile(effectorP_file_path):
        print('ERROR - missing:', effectorP_file)
    print(sample,"- All good, ready to roll")
    #Change this to where outputfiles shall go
    candidate_output_dir = output_path / sample / output

    try:
        os.makedirs(candidate_output_dir)
    except Exception:
        pass

    #delete content of folder to avoid appending files
    for file in os.listdir(candidate_output_dir):
        file_path = os.path.join(candidate_output_dir, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            print(e)

    proteins = {}
    deepsig = {}
    signalP = {}
    signalPandDeepsig = {}

    effectorP = {}
    interproScan = {}
    cysteins = {}
    proteinMW = {}
    matureProtein = {}

    everything = {}
    candidates = {}
    cazymes = {}

    #Preprocessing
    print(sample,'- Processing all the raw data')

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
                effectorP[name] = 0
                if cazyme_file_present:
                    cazymes[name] = 'No'
                else:
                    cazymes[name] = 'NA'

            else:
                proteins[name] += line

    print(sample,'- processed proteins data')

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

    #Read dbCAN2 CAZyme output
    if Path(cazyme_file_path).is_file():
        with open(cazyme_file_path) as file:
            print(sample,'- reading CAZyme file')
            input = csv.reader(file, delimiter='\t')
            next(input, None)
            for line in input:
                if int(line[5]) >= 2:
                    cazymes[line[0]] = 'Yes'
                else:
                    cazymes[line[0]] = 'No'
            print(sample,'- processed dbCAN cazyme data')

    else:
        print('no cazymes file')


    #Read deepsig output
    with open(deepsig_file_path) as file:
        print(sample,'- reading deepsig file')
        input = csv.reader(file, delimiter='\t')
        for line in input:
            if line[1] == 'SignalPeptide':
                deepsig[line[0]] = int(line[3])


    print(sample,'- processed deepsig data')

    #Read effectorP output
    with open(effectorP_file_path) as file:
        print(sample,'- reading effectorP file')
        input = csv.reader(file, delimiter='\t')
        for line in input:
            if line[1] == 'Effector':
                effectorP[line[0]] = float(line[2])

    print(sample, '- processed effectorP data')

    #calculate cystein content and MW in kDA
    def cysteineContent(aaSequence):
        cystein = 0
        proteinMass = 0.0
        weights = dict(
            A = 71.03711,
            C = 103.00919,
            D = 115.02694,
            E = 129.04259,
            F = 147.06841,
            G = 57.02146,
            H = 137.05891,
            I = 113.08406,
            K = 128.09496,
            L = 113.08406,
            M = 131.04049,
            N = 114.04293,
            P = 97.05276,
            Q = 128.05858,
            R = 156.10111,
            S = 87.03203,
            T = 101.04768,
            V = 99.06841,
            W = 186.07931,
            Y = 163.06333)
        for aminoAcid in aaSequence:
            proteinMass += weights[aminoAcid]
            if aminoAcid == 'C':
                cystein += 1
        return cystein, round(proteinMass/1000,3)

    #Filling up 'everything' dictionary
    print(sample, '- Putting everything in one place')
    for e in proteins:
        everything[e] = ('deepsig', deepsig[e]), ('SignalP', signalP[e]), ('EffectorP', effectorP[e]),('Cazyme', cazymes[e]), ('Sequence', proteins[e])

    #Effector Prediction
    # Combine predictions where both SignalP and deepsig predict a signal peptide,
    # but if they can't agree on the cleavage site prefer SignalP over deepsig
    for key, value in signalP.items():
        if value > 0 and deepsig[key] > 0:
            if value == deepsig[key]:
                signalPandDeepsig[key] = value
            else:
                signalPandDeepsig[key] = value


    #Process cystein content and protein mass for candidates with signal peptides
    print(sample,'- calculating cystein content and MW for candidates')
    for key, value in signalPandDeepsig.items():
        mature = proteins[key][value:]
        temp = cysteineContent(mature)
        cysteins[key], proteinMW[key] = temp[0], temp[1]
        matureProtein[key] = mature

    print(sample,'- processed cystein content and MW for candidates')

    #defining ouput files
    output_name_candidates = sample+'.candidates.txt'
    output_file_candidates = candidate_output_dir / output_name_candidates

    output_name_candidates_mature_fasta = sample+'.candidates.mature.fasta'
    output_file_candidates_mature_fasta = candidate_output_dir / output_name_candidates_mature_fasta

    output_name_candidates_fasta = sample+'.candidates.fasta'
    output_file_candidates_fasta = candidate_output_dir / output_name_candidates_fasta

    #everything: SignalP[0/1], EffectorP[float], Cysteins[int], MW[float], aaSequence[string]
    for element in matureProtein:
        #If effectorP score > 0:
        if effectorP[element] >= effectorPscore:
            if proteinMW[element] <= MWcutoff:
                if cysteins[element] >= Cyscutoff:
                    if element in interproScan:
                        candidates[element] = effectorP[element], cysteins[element], proteinMW[element], cazymes[element], interproScan[element], matureProtein[element]
                        print(element,'\t','SignalP=SP', '\t', 'EffectorP_Score:', effectorP[element],'\t', 'Cysteins:', cysteins[element], '\t', 'MW:', proteinMW[element], '\t', 'CAZyme:', cazymes[element], '\t', 'InterProScan_Results:', interproScan[element], file=open(os.path.join(candidate_output_dir, output_name_candidates),"a"))
                        print('>'+element, '\n'+matureProtein[element], file=open(os.path.join(candidate_output_dir, output_name_candidates_mature_fasta),'a'))
                        print('>'+element, '\n'+proteins[element], file=open(os.path.join(candidate_output_dir, output_name_candidates_fasta),'a'))
                    else:
                        candidates[element] = effectorP[element], cysteins[element], proteinMW[element], cazymes[element], 'No InterProScan Hits', matureProtein[element]
                        print(element,'\t','SignalP=SP', '\t', 'EffectorP_Score:', effectorP[element],'\t', 'Cysteins:', cysteins[element], '\t', 'MW:', proteinMW[element],  '\t', 'CAZyme:', cazymes[element], '\t', 'InterProScan_Results:', 'NONE', file=open(os.path.join(candidate_output_dir, output_name_candidates),"a"))
                        print('>'+element, '\n'+matureProtein[element], file=open(os.path.join(candidate_output_dir, output_name_candidates_mature_fasta),'a'))
                        print('>'+element, '\n'+proteins[element], file=open(os.path.join(candidate_output_dir, output_name_candidates_fasta),'a'))

    print('============================================')
    print(sample,'============= Summary ============')
    print(sample,'- Selection of effector candidates:')
    print(sample,'- has Signal peptide (SignalP & deepsig)')
    print(sample,'- EffectorP score >=', effectorPscore)
    print(sample,'- Cysteins >=', Cyscutoff)
    print(sample,'- MW of mature protein <=', MWcutoff, 'kDa')
    print(sample,'- Effector candidates:',len(candidates))
    print('============================================')
