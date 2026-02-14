import sys

def read_fasta(filename):
    sequence = ""
    with open(filename, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence


def affine_gap_alignment(seq1, seq2, match, mismatch, gap_open, gap_extend):
    n = len(seq1)
    m = len(seq2)

    NEG_INF = float('-inf')

    # Initialize matrices
    M = [[NEG_INF]*(m+1) for _ in range(n+1)]
    Ix = [[NEG_INF]*(m+1) for _ in range(n+1)]
    Iy = [[NEG_INF]*(m+1) for _ in range(n+1)]

    M[0][0] = 0

    # Initialize first column
    for i in range(1, n+1):
        Ix[i][0] = gap_open + (i-1)*gap_extend

    # Initialize first row
    for j in range(1, m+1):
        Iy[0][j] = gap_open + (j-1)*gap_extend

    # Fill matrices
    for i in range(1, n+1):
        for j in range(1, m+1):

            score = match if seq1[i-1] == seq2[j-1] else mismatch

            M[i][j] = max(
                M[i-1][j-1],
                Ix[i-1][j-1],
                Iy[i-1][j-1]
            ) + score

            Ix[i][j] = max(
                M[i-1][j] + gap_open,
                Ix[i-1][j] + gap_extend
            )

            Iy[i][j] = max(
                M[i][j-1] + gap_open,
                Iy[i][j-1] + gap_extend
            )

    # Traceback
    align1 = ""
    align2 = ""

    i, j = n, m
    final_score = max(M[n][m], Ix[n][m], Iy[n][m])

    if final_score == M[n][m]:
        current = "M"
    elif final_score == Ix[n][m]:
        current = "Ix"
    else:
        current = "Iy"

    while i > 0 or j > 0:

        if current == "M":
            score = match if seq1[i-1] == seq2[j-1] else mismatch

            if M[i][j] == M[i-1][j-1] + score:
                current = "M"
            elif M[i][j] == Ix[i-1][j-1] + score:
                current = "Ix"
            else:
                current = "Iy"

            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1

        elif current == "Ix":
            if Ix[i][j] == M[i-1][j] + gap_open:
                current = "M"
            else:
                current = "Ix"

            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        elif current == "Iy":
            if Iy[i][j] == M[i][j-1] + gap_open:
                current = "M"
            else:
                current = "Iy"

            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1

    return align1, align2, final_score


if __name__ == "__main__":

    if len(sys.argv) != 7:
        print("Usage: python3 affine_align.py seq1.fasta seq2.fasta match mismatch gap_open gap_extend")
        sys.exit(1)

    seq1 = read_fasta(sys.argv[1])
    seq2 = read_fasta(sys.argv[2])

    match = int(sys.argv[3])
    mismatch = int(sys.argv[4])
    gap_open = int(sys.argv[5])
    gap_extend = int(sys.argv[6])

    alignment1, alignment2, score = affine_gap_alignment(
        seq1, seq2, match, mismatch, gap_open, gap_extend
    )

    print("\n===== AFFINE GAP GLOBAL ALIGNMENT =====\n")
    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nBest Alignment:\n")
    print(alignment1)
    print(alignment2)
    print("\nFinal Alignment Score:", score)

