# KFERQminer.

from Bio import SeqIO  # reading fasta files
import re  # counting aminoacids for motif approval
import csv # reading motif database


# is_bonafide function takes 5aa subsequence from a protein as an argument and checks if it meets the rules of canonical bonafide KFERQ motif, returns True or False
def is_bonafide(motif):
    is_motif = False
    side = ["Q", "N"]
    positive = ["K", "R"]
    negative = ["D", "E"]
    hydrophobic = ["F", "L", "I", "V"]

    if any(aminoacid not in side+positive+negative+hydrophobic for aminoacid in motif):
        return False
    if motif[0] in side or motif[-1] in side:
        if len(re.findall('[' + ''.join(side) + ']', motif)) == 1:
            if len(re.findall('['+''.join(positive)+']', motif)) in [1, 2]:
                if len(re.findall('['+''.join(negative)+']', motif)) == 1:
                    if len(re.findall('['+''.join(hydrophobic)+']', motif)) in [1, 2]:
                        is_motif = True
    return is_motif


# is_phosphorylated function takes 5aa subsequence from a protein as an arguemnt and checks if it meets the rules of phosphorylated KFERQ-like motif, returns True or False
def is_phosphorylated(motif):
    is_motif = False
    side = ["Q", "N"]
    positive = ["K", "R"]
    negative = ["T", "Y", "S"]
    hydrophobic = ["F", "L", "I", "V"]

    if any(aminoacid not in side+positive+negative+hydrophobic for aminoacid in motif):
        return False
    if motif[0] in side or motif[-1] in side:
        if len(re.findall('[' + ''.join(side) + ']', motif)) == 1:
            if len(re.findall('['+''.join(positive)+']', motif)) in [1, 2]:
                if len(re.findall('['+''.join(negative)+']', motif)) == 1:
                    if len(re.findall('['+''.join(hydrophobic)+']', motif)) in [1, 2]:
                        is_motif = True
    return is_motif


# is_acetylated function takes 5aa subsequence from a protein as an arguemnt and checks if it meets the rules of acetylated KFERQ-like motif, returns True or False
def is_acetylated(motif):
    is_motif = False
    side = ["K"]
    positive = ["K", "R"]
    negative = ["D", "E"]
    hydrophobic = ["F", "L", "I", "V"]

    if any(aminoacid not in side+positive+negative+hydrophobic for aminoacid in motif):
        return False
    if motif[0] in side or motif[-1] in side:
        if len(re.findall('['+''.join(positive)+']', motif)) in [1, 2, 3]:
            if len(re.findall('['+''.join(negative)+']', motif)) == 1:
                if len(re.findall('['+''.join(hydrophobic)+']', motif)) in [1, 2]:
                    is_motif = True
    return is_motif

# charge_transformation transforms any amino acid sequence into charge string
def charge_transformation(seq):
    side = ["Q", "N"]
    positive = ["K", "R"]
    negative = ["D", "E", "T", "Y", "S"]
    hydrophobic = ["F", "L", "I", "V"]

    seq_charge = ''
    for aa in seq:
      if aa in side:
        seq_charge = seq_charge + '!'
      elif aa in positive:
        seq_charge = seq_charge + '+'
      elif aa in negative:
        seq_charge = seq_charge + '-'
      elif aa in hydrophobic:
        seq_charge = seq_charge + '^'
    return seq_charge


# The program starts
output = open("output.txt", "w")
input_f = 'input.fasta'

for record in SeqIO.parse(input_f, "fasta"):
    sequence = str(record.seq)
    
    # get 5_mers from the sequence
    n_5_mers = len(sequence) - 4
    for i in range(n_5_mers):
        mer = sequence[i:i+5]
        # check mer
        if is_bonafide(mer):
            result = [str(record.description), mer, 'canonical', str(i+1) + '\n']
            output.write(','.join(result))
        if is_phosphorylated(mer):
            result = [str(record.description), mer, 'phosphorylated', str(i+1) + '\n']
            output.write(','.join(result))
        if is_acetylated(mer):
            result = [str(record.description), mer, 'phosphorylated', str(i+1) + '\n']
            output.write(','.join(result))

output.close()

# checking validated motifs from db
output = open('output_validated.txt', 'w')
head = ['record_id', 'motif', 'position', 'validation', 'protein', 'author',\
        'year', 'doi']

output.write(','.join(head))
for record in SeqIO.parse('input.fasta', 'fasta'):
    sequence = str(record.seq)

    with open('motif_db.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')

        # identical
        for row in csv_reader:
          motif = row[0]
          n_mers = len(sequence) - len(motif) + 1
          for i in range(n_mers):
            mer = sequence[i:i+5]
            if mer == motif:
              result = [str(record.description), mer, str(i+1), 'validated',\
                        row[2], row[3], str(row[4]), row[5] + '\n']
              output.write(','.join(result))

output.close()

# checking motifs similar to validated from db
output = open('output_similar.txt', 'w')
head = ['record_id', 'motif', 'position', 'validation', 'protein', 'author',\
        'year', 'doi']
output.write(','.join(head))
for record in SeqIO.parse('input.fasta', 'fasta'):
    sequence = str(record.seq)

    with open('motif_db.txt') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=';')

        # similar
        for row in csv_reader:
          motif = row[1]
          n_mers = len(sequence) - len(motif) + 1
          for i in range(n_mers):
            mer = sequence[i:i+5]
            if charge_transformation(mer) == motif:
              result = [str(record.description), mer, str(i+1), 'similar',\
                        row[2], row[3], str(row[4]), row[5] + '\n']
              output.write(','.join(result))

output.close()
