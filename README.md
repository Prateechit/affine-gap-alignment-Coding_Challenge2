Affine Gap Global Sequence Alignment Tool
-
This project implements a global nucleotide sequence alignment tool using the Needleman–Wunsch algorithm with affine gap penalties. The program aligns two DNA sequences provided in FASTA format and computes the optimal global alignment using user-defined scoring parameters.

Overview:
-
The tool performs global alignment between two nucleotide sequences by applying an affine gap penalty model. Unlike simple linear gap penalties, affine gap scoring differentiates between opening a gap and extending an existing gap. This approach more accurately represents biological insertions and deletions.

Features:

1. Global sequence alignment

2. Affine gap penalty implementation

Three dynamic programming matrices:

1. M → Match/Mismatch

2. Ix → Gap in sequence 1

3. Iy → Gap in sequence 2

Custom scoring parameters provided by the user

Outputs aligned sequences and final alignment score

Works for sequences of any length

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Input
-
The program requires:

1. Two FASTA files containing nucleotide sequences

2. Scoring parameters provided as command-line inputs:

-- Match score

-- Mismatch penalty

-- Gap opening penalty

-- Gap extension penalty

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
How to Run

Run the program using:

python3 affine_align.py seq1.fasta seq2.fasta match mismatch gap_open gap_extend

Example:

python3 affine_align.py seq1_100bp.fasta seq2_100bp.fasta 2 -1 -5 -1

In this example:
Match score = 2
Mismatch penalty = -1
Gap opening penalty = -5
Gap extension penalty = -1

-------------------------------------------------------------------------------------------------------------------------------------------------------
Output
-
The program prints:

The best global alignment

The final alignment score

Example output structure:

===== BEST GLOBAL ALIGNMENT (Affine Gap) =====

Aligned Sequence 1
Aligned Sequence 2

Final Alignment Score: 176

-----------------------------------------------------------------------------------------------------------------
Algorithm Description

This tool implements the Needleman–Wunsch global alignment algorithm using an affine gap penalty model.

The affine gap penalty is calculated as:

Gap penalty = gap_open + (k - 1) × gap_extend

Where:
gap_open = penalty for starting a gap
gap_extend = penalty for extending an existing gap
k = length of the gap

Three dynamic programming matrices are used:

M stores scores for match/mismatch
Ix stores scores when there is a gap in sequence 1
Iy stores scores when there is a gap in sequence 2

This ensures that long gaps are penalized realistically by charging a higher cost for opening a gap and a smaller cost for extending it.

------------------------------------------------------------------------------------------------------------------------------------------
Time and Space Complexity

Time Complexity: O(nm)
Space Complexity: O(nm)

Where n and m are the lengths of the two sequences.

------------------------------------------------------------------------------------------------------------
Project Structure
-
affine_align.py – main alignment script

seq1.fasta – first input sequence

seq2.fasta – second input sequence

output.txt – saved alignment result

README – project documentation

-------------------------------------------------------------------------------------
Applications

Comparative genomics, Mutation analysis, Sequence similarity studies, Bioinformatics algorithm implementation
