#             A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
blosum62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
            [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
            [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
            [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
            [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
            [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
            [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
            [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
            [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
            [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
            [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
            [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
            [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
            [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
            [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4]]

amino_acids = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
               'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}


# DNA

def local_alignment_DNA(seqA, seqB):
    match = +1
    mismatch = -1
    gap = -2

    row = len(seqA)
    column = len(seqB)

    bond = ""

    max_list = []

    L = [[0 for x in range(row + 1)] for y in range(column + 1)]

    maximum = -1

    for i in range(1, column + 1, 1):
        for j in range(1, row + 1, 1):
            # the local alignment recurrence rule:
            L[i][j] = max(
                L[i][j - 1] + gap,
                L[i - 1][j] + gap,
                L[i - 1][j - 1] + (match if seqB[i - 1] == seqA[j - 1] else mismatch),
                0
            )
            if L[i][j] >= maximum:
                maximum = L[i][j]

    for i in range(1, column + 1, 1):
        for j in range(1, row + 1, 1):
            if L[i][j] == maximum:
                startTraceback = (i, j)
                max_list.append(startTraceback)
    # print(max_list)
    # print(L)

    for maxIndex in max_list:

        i, j = maxIndex
        localAlignmentSeqA = ""
        localAlignmentSeqB = ""

        while i > 0 and j > 0 and L[i][j] > 0:
            up = L[i - 1][j] + gap
            left = L[i][j - 1] + gap
            diagonal = L[i - 1][j - 1] + (match if seqB[i - 1] == seqA[j - 1] else mismatch)

            if L[i][j] == diagonal:
                localAlignmentSeqA += seqA[j - 1]
                localAlignmentSeqB += seqB[i - 1]
                i -= 1
                j -= 1

            elif L[i][j] == up:
                localAlignmentSeqA += "-"
                localAlignmentSeqB += seqB[i - 1]
                i -= 1

            else:
                localAlignmentSeqA += seqA[j - 1]
                localAlignmentSeqB += "-"
                j -= 1

        while i > 0:
            localAlignmentSeqA += "-"
            localAlignmentSeqB += seqB[i - 1]
            i -= 1

        while j > 0:
            localAlignmentSeqA += seqA[j - 1]
            localAlignmentSeqB += "-"
            j -= 1

        localAlignmentSeqA = localAlignmentSeqA[::-1]
        localAlignmentSeqB = localAlignmentSeqB[::-1]

        print(localAlignmentSeqA)
        print(localAlignmentSeqB)


# protein
def local_alignment_protein(seqA, seqB):
    gap = -1

    row = len(seqA)
    column = len(seqB)

    bond = ""

    max_list = []

    L = [[0 for x in range(row + 1)] for y in range(column + 1)]

    maximum = -1

    for i in range(1, column + 1, 1):
        for j in range(1, row + 1, 1):
            # the local alignment recurrence rule:
            L[i][j] = max(
                L[i][j - 1] + gap,
                L[i - 1][j] + gap,
                L[i - 1][j - 1] + (blosum62[amino_acids[seqB[i - 1]]][amino_acids[seqA[j - 1]]]),
                0
            )
            if L[i][j] >= maximum:
                maximum = L[i][j]

    for i in range(1, column + 1, 1):
        for j in range(1, row + 1, 1):
            if L[i][j] >= maximum:
                startTraceback = (i, j)
                max_list.append(startTraceback)
    # print(max_list)
    # print(L)

    for maxIndex in max_list:

        i, j = maxIndex
        localAlignmentSeqA = ""
        localAlignmentSeqB = ""

        while i > 0 and j > 0 and L[i][j] > 0:
            up = L[i - 1][j] + gap
            left = L[i][j - 1] + gap
            diagonal = L[i - 1][j - 1] + (blosum62[amino_acids[seqA[j - 1]]][amino_acids[seqB[i - 1]]])

            if L[i][j] == diagonal:
                localAlignmentSeqA += seqA[j - 1]
                localAlignmentSeqB += seqB[i - 1]
                i -= 1
                j -= 1

            elif L[i][j] == up:
                localAlignmentSeqA += "-"
                localAlignmentSeqB += seqB[i - 1]
                i -= 1

            else:
                localAlignmentSeqA += seqA[j - 1]
                localAlignmentSeqB += "-"
                j -= 1

        while i > 0:
            localAlignmentSeqA += "-"
            localAlignmentSeqB += seqB[i - 1]
            i -= 1

        while j > 0:
            localAlignmentSeqA += seqA[j - 1]
            localAlignmentSeqB += "-"
            j -= 1

        localAlignmentSeqA = localAlignmentSeqA[::-1]
        localAlignmentSeqB = localAlignmentSeqB[::-1]

        print(localAlignmentSeqA)
        print(localAlignmentSeqB)


option = input("Enter your option: \n 1: For DNA. \n 2: For Protein. \n")

if option == "1":
    seqA = input("Enter the first sequence: \n")
    seqB = input("Enter the second sequence: \n")
    local_alignment_DNA(seqA, seqB)

elif option == "2":
    seqA = input("Enter the first sequence: \n")
    seqB = input("Enter the second sequence: \n")
    local_alignment_protein(seqA, seqB)
else:
    print("Invalid Option\n")
