# Phylogenetic Tree Task


import string


# First we get the matrix
def enterMatrix(size):
    matrix = []
    print("Enter the distance matrix: ")
    for i in range(size):
        temp = []
        for j in range(size):
            if i > j:
                allInputs = float(input())
                temp.append(allInputs)
        matrix.append(temp)
    return matrix


# Then we get matrix labels
def get_labels(matrix):
  # Create a string of all uppercase letters
  alphabet_string = string.ascii_uppercase

  # Create a list of all uppercase letters
  alphabet_list = list(alphabet_string)

  # Get value
  labels = alphabet_list[0:(len(matrix))]

  return labels


# Then we locate the smallest cell in the matrix
def getSmallestCell(matrix):
  min = 1e9
  minRow = -100
  minCol = -100

  for row in range(len(matrix)):
    for col in range(len(matrix[row])):
      if matrix[row][col] < min:
        min = matrix[row][col]
        minRow = row
        minCol = col

  return  minRow, minCol, min


# Then, since we're working on the lower triangle of the matrix, we join the 2 rows intersecting at the minimum value in the upper row 
def join_matrix(matrix, minRow, minCol):
    # Getting the upper and lower rows
    upperRow = min(minRow, minCol)
    lowerRow = max(minRow, minCol)

    # Reconstructing the upper row
    newRow = []
    for i in range(0, upperRow):
        newRow.append((matrix[upperRow][i] + matrix[lowerRow][i])/2)
    matrix[upperRow] = newRow

    # Updating the cells intersecting with the upper row's column
    # Since this works on the lower triangle, so the rows between the upper and lower rows will have a cell in the lower row containing its intersection value
    # while the rows below won't, so we will then access their cells intersecting with the lower row
    for i in range(upperRow+1, len(matrix)):
        # The rows between the upper and lower rows
        if i < lowerRow:
            matrix[i][upperRow] = (matrix[i][upperRow]+matrix[lowerRow][i])/2

        # The rows below the lower row
        elif i > lowerRow:
            matrix[i][upperRow] = (matrix[i][upperRow]+matrix[i][lowerRow])/2

            # Deleting the value of intersection with lower row
            del matrix[i][lowerRow]

    # Deleting the lower row
    del matrix[lowerRow]

    # Printing updated matrix
    for row in range(len(matrix)):
        for col in range(len(matrix[row])):
            print(matrix[row][col], end=" ")
        print()
    print()


# Joining the closest pair
def join_pair(labels, mini, minRow, minCol):
    # Getting the upper and lower rows
    upperRow = min(minRow, minCol)
    lowerRow = max(minRow, minCol)

    # Join the labels in the upperRow label
    labels[upperRow] = "(" + labels[upperRow] + "," + labels[lowerRow] + ")"

    # Remove the redundant lower row label
    del labels[lowerRow]
    print(labels[upperRow], mini)

def UPGMA(matrix, labels):
    # Run until all pairs are joined
    while len(labels) > 1:
        # Locating the smallest cell in the matrix
        x, y, minimum = getSmallestCell(matrix)

        # Join the matrix on the cell co-ordinates
        join_matrix(matrix, x, y)

        # Update the tree
        join_pair(labels, minimum, x, y)

    # Printing the final tree
    print("Final Tree: ", labels[0])


size = int(input("Enter the matrix size: "))

matrix = enterMatrix(size)

labels = get_labels(matrix)

UPGMA(matrix, labels)


#Example:
#size = 5
#0.1, 0.8, 0.5, 0.7, 0.6, 0.3, 1, 0.9, 0.4, 0.2

