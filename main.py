# KFERQminer.

from Bio import SeqIO

# motif_counts is dictionary storaging information about how many motifs of different categories were found in input file
# index 0 refers to every potential motif
# index 1 refers to potential motifs which is similar based on aminoacid charge order to motif already validated and storaged in database
# index 2 refers to motifs already validated and storaged in database

global motif_counts

motif_counts = {
    "canonical": [0, 0, 0],
    "phosphorylation": [0, 0, 0],
    "acetylation": [0, 0, 0]
}


# is_canonical function takes 5aa subsequence from a protein as an arguemnt and checks if it meets the rules of canonical KFERQ motif

def is_canonical(sequence):
    output.write("CANONICAL:\n")
    side = ["Q", "N"]
    positive = ["K", "R"]
    negative = ["D", "E"]
    hydrophobic = ["F", "L", "I", "V"]
    a = 0
    b = 5
    while b <= len(sequence):
        motif = sequence[a:b]
        motif = str(motif)
        if motif[0] in side or motif[-1] in side:
            p = 0
            n = 0
            h = 0
            s = 0
            other = 0
            for aminoacid in motif:
                if aminoacid in positive:
                    p += 1
                elif aminoacid in negative:
                    n += 1
                elif aminoacid in hydrophobic:
                    h += 1
                elif aminoacid in side:
                    s += 1
                else:
                    other += 1
                    break
            if s != 1 or p not in range(1, 3) or n != 1 or h not in range(1, 3) or other != 0:
                # the sequence does not meet the canonical motif rules
                pass
            else:
                # we have found a potential canonical KFERQ motif
                motif_counts["canonical"][0] += 1
                # translation of motif into charge order (quasi_motif) to check if similar motifs were already validated
                quasi_motif = motif
                for aa in positive:
                    quasi_motif = quasi_motif.replace(aa, "+")
                for aa in negative:
                    quasi_motif = quasi_motif.replace(aa, "-")
                for aa in hydrophobic:
                    quasi_motif = quasi_motif.replace(aa, "^")
                for aa in side:
                    quasi_motif = quasi_motif.replace(aa, "!")
                is_validated(motif, quasi_motif, a, "canonical")
        a += 1
        b += 1


# is_phosphorylated function takes 5aa subsequence from a protein as an arguemnt and checks if it meets the rules of phosphorylated KFERQ-like motif

def is_phosphorylated(sequence):
    output.write("PHOSPHORYLATION (T/S/Y):\n")
    side = ["Q", "N"]
    positive = ["K", "R"]
    negative = ["Y", "T", "S"]
    hydrophobic = ["F", "L", "I", "V"]
    a = 0
    b = 5
    while b <= len(sequence):
        motif = sequence[a:b]
        motif = str(motif)
        if motif[0] in side or motif[-1] in side:
            p = 0
            n = 0
            h = 0
            s = 0
            other = 0
            for aminoacid in motif:
                if aminoacid in positive:
                    p += 1
                elif aminoacid in negative:
                    n += 1
                elif aminoacid in hydrophobic:
                    h += 1
                elif aminoacid in side:
                    s += 1
                    if s > 1:
                        break
                else:
                    other += 1
                    break
            if s != 1 or p not in range(1, 3) or n != 1 or h not in range(1, 3) or other != 0:
                pass
                # the sequence does not meet the phosphorylated KFERQ-like motif rules
            else:
                # we have found a potential phosphorylated KFERQ-like motif
                motif_counts["phosphorylation"][0] += 1
                # translation of motif into charge order (quasi_motif) to check if similar motifs were already validated
                quasi_motif = motif
                for aa in positive:
                    quasi_motif = quasi_motif.replace(aa, "+")
                for aa in negative:
                    quasi_motif = quasi_motif.replace(aa, "-")
                for aa in hydrophobic:
                    quasi_motif = quasi_motif.replace(aa, "^")
                for aa in side:
                    quasi_motif = quasi_motif.replace(aa, "!")
                is_validated(motif, quasi_motif, a, "phosphorylation")
        a += 1
        b += 1


# is_acetylated function takes 5aa subsequence from a protein as an arguemnt and checks if it meets the rules of acetylated KFERQ-like motif

def is_acetylated(sequence):
    output.write("ACETYLATION (side K):\n")
    side = ["K"]
    positive = ["K", "R"]
    negative = ["D", "E"]
    hydrophobic = ["F", "L", "I", "V"]
    a = 0
    b = 5
    while b <= len(sequence):
        motif = sequence[a:b]
        motif = str(motif)
        if motif[0] in side or motif[-1] in side:
            p = 0
            n = 0
            h = 0
            s = 0
            other = 0
            for aminoacid in motif:
                if aminoacid in positive:
                    p += 1
                elif aminoacid in negative:
                    n += 1
                elif aminoacid in hydrophobic:
                    h += 1
                else:
                    other += 1
                    break
            if p not in range(2, 4) or n != 1 or h not in range(1, 3) or other != 0:
                pass
                # the sequence does not meet the acetylated KFERQ-like motif rules
            else:
                # we have found a potential acetylated KFERQ-like motif
                motif_counts["acetylation"][0] += 1
                # translation of motif into charge order (quasi_motif) to check if similar motifs were already validated
                quasi_motif = motif
                for aa in positive:
                    quasi_motif = quasi_motif.replace(aa, "+")
                for aa in negative:
                    quasi_motif = quasi_motif.replace(aa, "-")
                for aa in hydrophobic:
                    quasi_motif = quasi_motif.replace(aa, "^")
                for aa in side:
                    quasi_motif = quasi_motif.replace(aa, "!")
                is_validated(motif, quasi_motif, a, "acetylation")
        a += 1
        b += 1


# is_other function takes 5aa subsequence from a protein as an arguemnt and checks if it was already validated as a functional motif frome those which does not meet the rules of KFERQ-like motifs

def is_other(sequence):
    validated_motifs = open("validated_motifs_other.txt", "r")
    validated_motifs = validated_motifs.readlines()[1:]
    a = 0
    b = 5
    while b <= len(sequence):
        motif = sequence[a:b]
        motif = str(motif)
        references = []
        for record in validated_motifs:
            line = record.replace("\n", "").split("\t")
            if motif in line[2:] or motif[::-1] in line[2:]:
                # we have found a motif which was validated
                references.append(str(line[0]) + ", " + str(line[1]))
        if len(references) > 0:
            output.write("Other:\n")
            output.write(str(motif) + " (" + str(a + 1) + ") *\n")
            output.write("This motif was validated as functional in:\n")
            for reference in references:
                output.write(reference + "\n")
        a += 1
        b += 1


# is_validated function takes 5aa already checked motif (motif), its charge order (quasi_motif), localization in protein sequence (start), and category of protein eg. canonical or phosphorylated (category) as arguemnts and checks if the motif or its charge order are already storaged in the KFERQminer database

def is_validated(motif, quasi_motif, start, category):
    motif = str(motif)
    validated_motifs = open("validated_motifs_" + category + ".txt", "r")
    validated_motifs = validated_motifs.readlines()[1:]
    references = []
    for record in validated_motifs:
        line = record.replace("\n", "").split("\t")
        if motif in line[2:] or motif[::-1] in line[2:]:
            # the program found motif in KFERQminer database; it is validated
            references.append(str(line[0]) + ", " + str(line[1]))

    validated_quasi_motifs = open("validated_quasi_motifs_" + category + ".txt", "r")
    validated_quasi_motifs = validated_quasi_motifs.readlines()[1:]
    references_charge = []
    quasi_motif = str(quasi_motif)
    for record in validated_quasi_motifs:
        line = record.replace("\n", "").split("\t")
        if quasi_motif in line[2:] or quasi_motif[::-1] in line[2:]:
            # the program found motif charge order in KFERQminer database; the motif is similar to another validated motif
            references_charge.append(str(line[0]) + ", " + str(line[1]))

    if len(references) > 0 and len(references_charge) > 0:
        # both validated motifs and similar motif charge orders were found in database
        motif_counts[category][2] += 1
        output.write(motif + " (" + str(start + 1) + ") *\n")
        output.write("This motif was validated as functional in:\n")
        for reference in references:
            output.write(reference + "\n")
        if len([reference for reference in references_charge if reference not in references]) > 0:
            output.write("This motif has the same charge order as validated motif in:\n")
            for reference in references_charge:
                if reference in references:
                    pass
                else:
                    output.write(reference + "\n")
    elif len(references_charge) > 0:
        # only the motif charge order was found in the KFERQminer database
        motif_counts[category][1] += 1
        output.write(motif + " (" + str(start + 1) + ") #\n")
        output.write("This motif has the same charge order as validated motif in:\n")
        for reference in references_charge:
            output.write(reference + "\n")
    else:
        output.write(motif + " (" + str(start + 1) + ")\n")


# The program starts

output = open("output.txt", "w")

with open("input.fasta", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        output.write(str(record.id) + "\n")
        sequence = record.seq
        is_canonical(sequence)
        print("canonical done")
        is_phosphorylated(sequence)
        print("phosphorylation done")
        is_acetylated(sequence)
        print("acetylation done")
        is_other(sequence)
        print("other done")
        output.write("\n")
        print(str(record.id) + " done")

output.close()

for category in motif_counts:
    print(category, motif_counts[category])

