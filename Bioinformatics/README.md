# Structural Genomics


### Local Alignment:

A python script to perform local alignment on two DNA or protein sequences.
User chooses the type of molecule to align with the DNA scoring:

 - Match: +1
 - Mismatch: -2
 - Gap:-1

and protein using BLOSUM62 similarity matrix


### Global Alignment:

A python script to perform global alignment on two DNA or protein sequences.
User chooses the type of molecule to align with the DNA scoring:

 - Match: +1
 - Mismatch: -2
 - Gap:-1

and protein using BLOSUM62 similarity matrix


### Phylogenetic Tree:

A python script to construct a phylogenetic tree given an input distance matrix using the UPGMA method.
User inserts the distance matrix and the output clarifies all the intermediate calculations, merging steps, and matrices alongside the final result.


### Blast Project:

A python script for the Blast technique to search for a protein query in a protein database.
User inputs a protein sequence query, word threshold, word length, and HSP threshold.
The output contains:
1. HSP(s) with their score if they exist. Along with the ID of the sequence in the database that the hsp was found in (For example, the first sequence would have ID of zero, the second sequence would have ID of 1...etc)
2. If not, then output a message indicating that the query wasn’t found.
