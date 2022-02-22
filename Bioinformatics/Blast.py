#             A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   X
blosum62 = [[ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -4],
            [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -4],
            [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3, -4],
            [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3, -4],
            [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -4], 
            [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2, -4],
            [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2, -4],
            [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -4],
            [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3, -4],
            [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -4],
            [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4],
            [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2, -4],
            [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -4],
            [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -4],
            [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -4],
            [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2, -4],
            [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -4],
            [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4],
            [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -4],
            [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -4],
            [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1]]

amino_acids = { 'A': 0, 'R': 1, 'N': 2, 'D' : 3, 'C' : 4, 'Q' : 5, 'E' : 6, 'G' : 7, 'H' : 8, 'I' : 9, 'L' : 10, 'K' : 11, 'M' : 12, 'F' : 13, 'P' : 14, 'S' : 15, 'T' : 16, 'W' : 17, 'Y' : 18, 'V' : 19, 'X' : 20 }

query_seq = input ("Enter protein sequence query: ")
word_threshold = int(input("Enter word threshold: "))
word_length = int(input("Enter word length: "))
HSP_threshold = int(input("Enter HSP threshold: "))

def local_alignment_protein(seqA, seqB, index):
    gap = -1
    row = len(seqA)
    column = len(seqB)
    bond = ""
    max_list = []
    test = -1
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

    for maxIndex in max_list:
        i, j = maxIndex
        localAlignmentSeqA = ""
        localAlignmentSeqB = ""

        while i > 0 and j > 0 and L[i][j] > 0:
            up = L[i - 1][j] + gap
            left = L[i][j - 1] + gap
            diagonal = L[i - 1][j - 1] + (blosum62[amino_acids[seqA[j - 1]]][amino_acids[seqB[i - 1]]])

            if L[i][j] == diagonal:
                if len(localAlignmentSeqB) == len(seqA)-1:
                  test = i-1
                localAlignmentSeqA += seqA[j - 1]
                localAlignmentSeqB += seqB[i - 1]
                i -= 1
                j -= 1

            else:
              break

        localAlignmentSeqA = localAlignmentSeqA[::-1]
        localAlignmentSeqB = localAlignmentSeqB[::-1]

        res = localAlignmentSeqA.find('-')
        if (res <= 0 and (localAlignmentSeqA == localAlignmentSeqB)) and len(localAlignmentSeqA) == len(seqA):
          protein_indices = {}
          protein_indices[seqA] = [test, maximum, index]
          return protein_indices

seq = []
for i in query_seq:
  seq.append(blosum62[amino_acids[i]][amino_acids[i]])

tmp_string = ''

count = 0
count2 = 0
flag = False

for i in range(len(query_seq)):
  if i+3 < len(query_seq):
    if seq[i] == seq[i+1]:
      count += 1
    elif seq[i] == seq[i+2] and seq[i+1] == seq[i+3]:
      count2 +=1
    else:
      if count+1 >= 6:
        tmp_string = tmp_string + 'X'
        count = 0
      elif count2+1 >= 4:
        tmp_string = tmp_string + 'X'
        count2 = 0
        flag = True
      else:
        if flag == True:
          for x in range(2, count2+1):
            tmp_string = tmp_string + query_seq[i]
          count2 = 0
        else:
          for x in range(count2+1):
            tmp_string = tmp_string + query_seq[i]
          count2 = 0

  else:
    tmp_string = tmp_string + query_seq[i]

words_seq = []
for i in range(0,len(tmp_string)):
    if i+word_length    <= len(tmp_string):
        words_seq.append(tmp_string[i:i+word_length])

def listToString(s): 
    
    # initialize an empty string
    str1 = "" 
    
    # return string  
    return (str1.join(s))

seeds = []
seq = []
indices = []
dic = dict()
for x in range(word_length):
  seq.append("A")
index = 0

for word in words_seq:
  for aa in amino_acids:
    for i in range(word_length):
      index
      sum = blosum62[amino_acids[word[i]]][amino_acids[aa]]
      seq[i] = aa
      for j in range(word_length):
        if i != j:
          sum += blosum62[amino_acids[word[i]]][amino_acids[word[i]]]
          seq[j] = word[j]
      if sum >= word_threshold and listToString(seq) not in seeds:
        seeds.append(listToString(seq))
        indices.append(index)
        dic[str(listToString(seq))] = index
        
  index +=1

def getIndices(seed, protein_matches):
  start_query = dic[seed]
  end_query = start_query + len(seed) - 1
  start_db = (protein_matches[seed])[0]
  end_db = start_db + len(seed) - 1
  return start_query-1, end_query+1, start_db-1, end_db+1

def getMatches():
  protein_matches = {}
  dataBase = []
  file = open('database.txt', 'r')
  for line in file:
    line = line[0:(len(line)-1)]
    dataBase.append(line)
  file.close()
  matches = []
  print('\nSeeds that successfully matched with DB sequences:\n')
  for sequence in seeds:
    for i in range(len(dataBase)):
      x = local_alignment_protein(sequence, dataBase[i], i)
      y = 5
      if x is not None:
        print(x)
        protein_matches.update(x)
        start_query, end_query, start_db, end_db = getIndices(sequence, protein_matches)
        flag1, flag2 = True, True
        score = protein_matches[sequence][1]
        newScore = score
        tmpScore = score
        match = sequence
        while (start_query >= 0 and start_db >= 0 and flag1) or (end_query < len(query_seq) and end_db < len(dataBase[i]) and  flag2):
          if(flag1 and start_query >= 0 and start_db >= 0):
            tmpScore += blosum62[amino_acids[query_seq[start_query]]][amino_acids[dataBase[i][start_db]]]
            if (newScore - tmpScore) > y:
              flag1 = False
            else:
              match = dataBase[i][start_db] + match
              newScore = tmpScore
              start_query -= 1
              start_db -= 1
          if(flag2 and end_query < len(query_seq) and end_db < len(dataBase[i])):
            newScore += blosum62[amino_acids[query_seq[end_query]]][amino_acids[dataBase[i][end_db]]]
            if (newScore - tmpScore) > y:
              flag2 = False
            else:
              match = match + dataBase[i][end_db]
              newScore = tmpScore
              end_query += 1
              end_db += 1
        if newScore >= HSP_threshold:
            matches.append([i, newScore, match])
  return matches

matches = getMatches()
matches.sort()
total = 0
flag1 = False
flag2 = False
print('\n**Results**\n')
if len(matches) > 0:
  for i in range(1, len(matches)):
    if matches[i][0] == matches[i-1][0] and matches[i][2] != matches[i-1][2]:
      if total == 0:
        total = matches[i][1] + matches[i-1][1]
        print('HPS: ', matches[i-1][2])
        print('HPS: ', matches[i][2])
      else:
        total += matches[i][1]
        print('HPS: ', matches[i][2])
    else:
      if flag1 == False:
        flag1 = True
      if i == len(matches)-1:
        flag2 = True
      if total > 0:
        print('Score: ', total)
        print('Data base ID: ', matches[i][0])
        total = 0
      else:
        print('HPS: ', matches[i-1][2])
        print('Score: ', matches[i-1][1])
        print('Data base ID: ', matches[i-1][0])
        if(i == len(matches)-1):
          print()
          print('HPS: ', matches[i][2])
          print('Score: ', matches[i][1])
          print('Data base ID: ', matches[i][0])
      print()
  if flag1 == False:
    print('Score: ', total)
    print('Data base ID: ', matches[0][0])
    print()
  elif flag2 == False:
    print('Score: ', total)
    print('Data base ID: ', matches[len(matches)-1][0])
    print()
else:
  print('The query wasn\'t found!')
