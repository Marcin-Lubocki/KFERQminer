# KFERQminer
Searching for KFERQ-like motifs in whole proteomes

The KFERQminer is a database of described in PubMed literature KFERQ-related motifs and bioinformatic tool that accelerate discovery of new KFERQ-bearing proteins.
Web KFERQminer is avaible at this link:

Motifs in the KFERQminer database
KFERQminer database contains KFERQ-related motifs described in the literature as functional. All motifs currently in it were collected from papers available in PubMed. For a motif to be included in the database, it had to be confirmed by at least one laboratory experiment. It means that in silico demonstration of the mere presence of a potential motif in the substrate is not sufficient.

To increase user’s comfort, the source paper is attached to every record in the KFERQminer database. If you find the KFERQminer useful in your research, please quote us.
If you find that a particular motif is missing from the database, please let us know.

KFERQminer tool structure
KFERQminer is a Python-based tool that uses a protein sequence in fasta format as an input from the user and subsequently searches for every possible KFERQ-like motif (based on canonical literature descriptions). Additionally, to narrow the results to the most probable true positives KFERQminer checks the input for the presence of motifs from the KFERQminer database. This way, KFERQminer is able to give the user 3 level scoring:
1. The user’s protein contains a putative motif based on the classical motif descriptions.
2. The user’s protein contains a motif that is similar to one described in literature. The claim is based on the amino acid charge order in a sequence, eg. hydrophobic amino acids are replaced with other hydrophobic amino acids.
3. The user’s protein contains a motif that is identical (the exact match) to a motif that has been already described in another study.
In the case of levels 2 and 3, due to the fact the source studies used different methodologies for verifying a motif, we decided to give the user information about the source paper and let him decide on his own about its utility in his research.
