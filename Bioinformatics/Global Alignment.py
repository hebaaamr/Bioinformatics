import numpy as np

option = input ("Enter your option: \n 1: For DNA. \n 2: For Protein. \n")

if option == "1":
    frist_sequence = input ("Enter the first sequence: \n")
    second_sequence = input ("Enter the second sequence: \n")
    input_matrix = np.zeros((len(frist_sequence)+1,len(second_sequence)+1))
    result_matrix = np.zeros((len(frist_sequence),len(second_sequence)))

    for i in range(len(frist_sequence)):
        for j in range(len(second_sequence)):
            if frist_sequence[i] == second_sequence[j]:
                result_matrix[i][j] = 1
            else:
                result_matrix[i][j] = -2

    for i in range(len(frist_sequence)+1):
        input_matrix[i][0] = i*-2
    for j in range(len(second_sequence)+1):
        input_matrix[0][j] = j*-2

    for i in range(1,len(frist_sequence)+1):
        for j in range(1,len(second_sequence)+1):
            input_matrix[i][j] = max(input_matrix[i-1][j-1] + result_matrix[i-1][j-1],
                                    input_matrix[i-1][j] + -2,
                                    input_matrix[i][j-1] + -2)

    frist_sequence_aligned = ""
    second_sequence_aligned = ""
    i = len(frist_sequence)
    j = len(second_sequence)

    while(i>0 and j>0):
        if (input_matrix[i][j] == input_matrix[i-1][j-1] + result_matrix[i-1][j-1] and i>0 and j>0):
            frist_sequence_aligned = frist_sequence[i-1] + frist_sequence_aligned
            second_sequence_aligned = second_sequence[j-1] + second_sequence_aligned
            i = i - 1
            j = j - 1    
        elif(input_matrix[i][j] == input_matrix[i-1][j] + -2 and i>0):
            frist_sequence_aligned = frist_sequence[i-1] + frist_sequence_aligned
            second_sequence_aligned = "-" + second_sequence_aligned
            i = i - 1
        else:
            frist_sequence_aligned = "-" + frist_sequence_aligned
            second_sequence_aligned = second_sequence[j-1] + second_sequence_aligned
            j = j - 1

    print("Frist sequence:  " + frist_sequence_aligned)
    print("Second sequence: " + second_sequence_aligned)

if option == "2":
    blosum62 = {
    '*':{'*':1,'A':-4,'C':-4,'B':-4,'E':-4,'D':-4,'G':-4,'F':-4,'I':-4,'H':-4,'K':-4,'M':-4,'L':-4,'N':-4,'Q':-4,'P':-4,'S':-4,'R':-4,'T':-4,'W':-4,'V':-4,'Y':-4,'X':-4,'Z':-4},
    'A':{'*':-4,'A':4,'C':0,'B':-2,'E':-1,'D':-2,'G':0,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-1,'S':1,'R':-1,'T':0,'W':-3,'V':0,'Y':-2,'X':-1,'Z':-1},
    'C':{'*':-4,'A':0,'C':9,'B':-3,'E':-4,'D':-3,'G':-3,'F':-2,'I':-1,'H':-3,'K':-3,'M':-1,'L':-1,'N':-3,'Q':-3,'P':-3,'S':-1,'R':-3,'T':-1,'W':-2,'V':-1,'Y':-2,'X':-1,'Z':-3},
    'B':{'*':-4,'A':-2,'C':-3,'B':4,'E':1,'D':4,'G':-1,'F':-3,'I':-3,'H':0,'K':0,'M':-3,'L':-4,'N':3,'Q':0,'P':-2,'S':0,'R':-1,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'E':{'*':-4,'A':-1,'C':-4,'B':1,'E':5,'D':2,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-2,'L':-3,'N':0,'Q':2,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4},
    'D':{'*':-4,'A':-2,'C':-3,'B':4,'E':2,'D':6,'G':-1,'F':-3,'I':-3,'H':-1,'K':-1,'M':-3,'L':-4,'N':1,'Q':0,'P':-1,'S':0,'R':-2,'T':-1,'W':-4,'V':-3,'Y':-3,'X':-1,'Z':1},
    'G':{'*':-4,'A':0,'C':-3,'B':-1,'E':-2,'D':-1,'G':6,'F':-3,'I':-4,'H':-2,'K':-2,'M':-3,'L':-4,'N':0,'Q':-2,'P':-2,'S':0,'R':-2,'T':-2,'W':-2,'V':-3,'Y':-3,'X':-1,'Z':-2},
    'F':{'*':-4,'A':-2,'C':-2,'B':-3,'E':-3,'D':-3,'G':-3,'F':6,'I':0,'H':-1,'K':-3,'M':0,'L':0,'N':-3,'Q':-3,'P':-4,'S':-2,'R':-3,'T':-2,'W':1,'V':-1,'Y':3,'X':-1,'Z':-3},
    'I':{'*':-4,'A':-1,'C':-1,'B':-3,'E':-3,'D':-3,'G':-4,'F':0,'I':4,'H':-3,'K':-3,'M':1,'L':2,'N':-3,'Q':-3,'P':-3,'S':-2,'R':-3,'T':-1,'W':-3,'V':3,'Y':-1,'X':-1,'Z':-3},
    'H':{'*':-4,'A':-2,'C':-3,'B':0,'E':0,'D':-1,'G':-2,'F':-1,'I':-3,'H':8,'K':-1,'M':-2,'L':-3,'N':1,'Q':0,'P':-2,'S':-1,'R':0,'T':-2,'W':-2,'V':-3,'Y':2,'X':-1,'Z':0},
    'K':{'*':-4,'A':-1,'C':-3,'B':0,'E':1,'D':-1,'G':-2,'F':-3,'I':-3,'H':-1,'K':5,'M':-1,'L':-2,'N':0,'Q':1,'P':-1,'S':0,'R':2,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':1},
    'M':{'*':-4,'A':-1,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':0,'I':1,'H':-2,'K':-1,'M':5,'L':2,'N':-2,'Q':0,'P':-2,'S':-1,'R':-1,'T':-1,'W':-1,'V':1,'Y':-1,'X':-1,'Z':-1},
    'L':{'*':-4,'A':-1,'C':-1,'B':-4,'E':-3,'D':-4,'G':-4,'F':0,'I':2,'H':-3,'K':-2,'M':2,'L':4,'N':-3,'Q':-2,'P':-3,'S':-2,'R':-2,'T':-1,'W':-2,'V':1,'Y':-1,'X':-1,'Z':-3},
    'N':{'*':-4,'A':-2,'C':-3,'B':3,'E':0,'D':1,'G':0,'F':-3,'I':-3,'H':1,'K':0,'M':-2,'L':-3,'N':6,'Q':0,'P':-2,'S':1,'R':0,'T':0,'W':-4,'V':-3,'Y':-2,'X':-1,'Z':0},
    'Q':{'*':-4,'A':-1,'C':-3,'B':0,'E':2,'D':0,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':0,'L':-2,'N':0,'Q':5,'P':-1,'S':0,'R':1,'T':-1,'W':-2,'V':-2,'Y':-1,'X':-1,'Z':3},
    'P':{'*':-4,'A':-1,'C':-3,'B':-2,'E':-1,'D':-1,'G':-2,'F':-4,'I':-3,'H':-2,'K':-1,'M':-2,'L':-3,'N':-2,'Q':-1,'P':7,'S':-1,'R':-2,'T':-1,'W':-4,'V':-2,'Y':-3,'X':-1,'Z':-1},
    'S':{'*':-4,'A':1,'C':-1,'B':0,'E':0,'D':0,'G':0,'F':-2,'I':-2,'H':-1,'K':0,'M':-1,'L':-2,'N':1,'Q':0,'P':-1,'S':4,'R':-1,'T':1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':0},
    'R':{'*':-4,'A':-1,'C':-3,'B':-1,'E':0,'D':-2,'G':-2,'F':-3,'I':-3,'H':0,'K':2,'M':-1,'L':-2,'N':0,'Q':1,'P':-2,'S':-1,'R':5,'T':-1,'W':-3,'V':-3,'Y':-2,'X':-1,'Z':0},
    'T':{'*':-4,'A':0,'C':-1,'B':-1,'E':-1,'D':-1,'G':-2,'F':-2,'I':-1,'H':-2,'K':-1,'M':-1,'L':-1,'N':0,'Q':-1,'P':-1,'S':1,'R':-1,'T':5,'W':-2,'V':0,'Y':-2,'X':-1,'Z':-1},
    'W':{'*':-4,'A':-3,'C':-2,'B':-4,'E':-3,'D':-4,'G':-2,'F':1,'I':-3,'H':-2,'K':-3,'M':-1,'L':-2,'N':-4,'Q':-2,'P':-4,'S':-3,'R':-3,'T':-2,'W':11,'V':-3,'Y':2,'X':-1,'Z':-3},
    'V':{'*':-4,'A':0,'C':-1,'B':-3,'E':-2,'D':-3,'G':-3,'F':-1,'I':3,'H':-3,'K':-2,'M':1,'L':1,'N':-3,'Q':-2,'P':-2,'S':-2,'R':-3,'T':0,'W':-3,'V':4,'Y':-1,'X':-1,'Z':-2},
    'Y':{'*':-4,'A':-2,'C':-2,'B':-3,'E':-2,'D':-3,'G':-3,'F':3,'I':-1,'H':2,'K':-2,'M':-1,'L':-1,'N':-2,'Q':-1,'P':-3,'S':-2,'R':-2,'T':-2,'W':2,'V':-1,'Y':7,'X':-1,'Z':-2},
    'X':{'*':-4,'A':-1,'C':-1,'B':-1,'E':-1,'D':-1,'G':-1,'F':-1,'I':-1,'H':-1,'K':-1,'M':-1,'L':-1,'N':-1,'Q':-1,'P':-1,'S':-1,'R':-1,'T':-1,'W':-1,'V':-1,'Y':-1,'X':-1,'Z':-1},
    'Z':{'*':-4,'A':-1,'C':-3,'B':1,'E':4,'D':1,'G':-2,'F':-3,'I':-3,'H':0,'K':1,'M':-1,'L':-3,'N':0,'Q':3,'P':-1,'S':0,'R':0,'T':-1,'W':-3,'V':-2,'Y':-2,'X':-1,'Z':4}}
    
    frist_sequence = input ("Enter the first sequence: \n")
    second_sequence = input ("Enter the second sequence: \n")
    input_matrix = np.zeros((len(frist_sequence)+1,len(second_sequence)+1))
    result_matrix = np.zeros((len(frist_sequence),len(second_sequence)))

    for i in range(len(frist_sequence)):
        for j in range(len(second_sequence)):
            if frist_sequence[i] == second_sequence[j]:
                result_matrix[i][j] = 1
            else:
                result_matrix[i][j] = -2

    for i in range(len(frist_sequence)+1):
        input_matrix[i][0] = -i
    for j in range(len(second_sequence)+1):
        input_matrix[0][j] = -j

    for i in range(1,len(frist_sequence)+1):
        for j in range(1,len(second_sequence)+1):
            value = blosum62[frist_sequence[i-1]][second_sequence[j-1]]
            input_matrix[i][j] = max(input_matrix[i-1][j-1] + value,
                                    input_matrix[i-1][j] + -2,
                                    input_matrix[i][j-1] + -2)

    frist_sequence_aligned = ""
    second_sequence_aligned = ""
    i = len(frist_sequence)
    j = len(second_sequence)

    while(i>0 and j>0):
        if (input_matrix[i][j] == input_matrix[i-1][j-1] + blosum62[frist_sequence[i-1]][second_sequence[j-1]] and i>0 and j>0):
            frist_sequence_aligned = frist_sequence[i-1] + frist_sequence_aligned
            second_sequence_aligned = second_sequence[j-1] + second_sequence_aligned
            i = i - 1
            j = j - 1    
        elif(input_matrix[i][j] == input_matrix[i-1][j] + -2 and i>0):
            frist_sequence_aligned = frist_sequence[i-1] + frist_sequence_aligned
            second_sequence_aligned = "-" + second_sequence_aligned
            i = i - 1
        else:
            frist_sequence_aligned = "-" + frist_sequence_aligned
            second_sequence_aligned = second_sequence[j-1] + second_sequence_aligned
            j = j - 1

    print("Frist sequence:  " + frist_sequence_aligned)
    print("Second sequence: " + second_sequence_aligned)
