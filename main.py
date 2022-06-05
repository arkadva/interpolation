import math


def create_matrix(x_list):
    i = 0
    pow = 0
    n = len(x_list)
    mat = [([0] * n) for i in range(n)]  # initialize the matrix with zeros
    for i in range(0, n):
        for j in range(0, n):
            mat[i][j] = math.pow(x_list[i], pow)
            pow += 1
        pow = 0
    return mat


def CalcMatrix(matrix, b):
    return Matrix_multiplication(InvertMatrix(matrix), b)


def Matrix_multiplication(mat1, mat2):
    if len(mat1[0]) != len(mat2):
        raise Exception("Illegal multiplication between matrix's ")
    result_mat = [([0] * len(mat2[0])) for i in range(len(mat1))]  # initialize the result matrix with zeros

    # iterate through the first matrix rows
    for row1 in range(0, len(mat1)):
        # iterate through the second matrix columns
        for col2 in range(0, len(mat2[0])):
            # iterate through the second matrix rows
            for row2 in range(0, len(mat2)):
                result_mat[row1][col2] += mat1[row1][row2] * mat2[row2][col2]
    return result_mat
def InvertMatrix(matrix):
    if len(matrix) != len(matrix[0]):
        raise Exception("singular matrix. there is no inverted matrix")
    n = len(matrix)
    inverted = Identity(n)
    for j in range(0, n):
        for i in range(0, n):
            if i == j:
                pivot = matrix[i][j]
                for k in range(i + 1, n):
                    if abs(matrix[k][j]) > abs(pivot):  # pivoting
                        elementary_matrix = ExchangeRows(k, i, n)
                        matrix = Matrix_multiplication(elementary_matrix, matrix)
                        inverted = Matrix_multiplication(elementary_matrix, inverted)
                        pivot = matrix[i][j]

                if matrix[i][j] == 0:
                    raise Exception("singular matrix. there is no inverted matrix")

        for i in range(0, n):
            if i != j:
                if matrix[i][j] != 0:
                    elementary_matrix = ResetOrgan(i, j, n, pivot, matrix[i][j])
                    matrix = Matrix_multiplication(elementary_matrix, matrix)
                    inverted = Matrix_multiplication(elementary_matrix, inverted)

    for i in range(0, n):
        if matrix[i][i] != 1:
            if matrix[i][i] < 0:
                elementary_matrix = MultiplyRow(i, -1, n)
                matrix = Matrix_multiplication(elementary_matrix, matrix)
                inverted = Matrix_multiplication(elementary_matrix, inverted)

            elementary_matrix = MultiplyRow(i, 1 / matrix[i][i], n)
            matrix = Matrix_multiplication(elementary_matrix, matrix)
            inverted = Matrix_multiplication(elementary_matrix, inverted)
    for row in range(n):
        for col in range(n):
            inverted[row][col] = round(inverted[row][col], 2)
    return inverted

def Identity(n):
    mat = [([0] * n) for i in range(n)]  # initialize the matrix with zeros
    for i in range(0, n):
        mat[i][i] = 1  # the identity matrix includes 1 all over its diagonal, starts at [0][0]
    return mat
def ResetOrgan(row, col, n, pivot, a):
    elementary_matrix = Identity(n)
    elementary_matrix[row][col] = -(a / pivot)
    return elementary_matrix

def MultiplyRow(row, a, n):
    elementary_matrix = Identity(n)
    elementary_matrix[row][row] = a
    return elementary_matrix

def ExchangeRows(row1, row2, n):
    elementary_matrix = Identity(n)
    elementary_matrix[row1][row1] = 0
    elementary_matrix[row1][row2] = 1
    elementary_matrix[row2][row2] = 0
    elementary_matrix[row2][row1] = 1
    return elementary_matrix

def solveMatrix(matrix, size):
    for i in range(size):
        # preprocess the matrix
        currentPivot = abs(matrix[i][i])
        maxRow = i
        row = i + 1
        while row < size:
            if abs(matrix[row][i]) > currentPivot:
                currentPivot = abs(matrix[row][i])
                maxRow = row
            row += 1
        matrix = swapRows(matrix, i, maxRow, size)

        matrix = setPivotToOne(matrix, i, i, size)
        row = i + 1
        while row < size:
            matrix = nullify(matrix, row, i, size, matrix[i][i])
            row += 1
    for i in range(1, size):
        row = i - 1
        while row >= 0:
            matrix = nullify(matrix, row, i, size, matrix[i][i])
            row -= 1

    return matrix


def identityMatrix(size):
    matrix = []
    for i in range(size):
        matrix.append([])
        for j in range(size):
            matrix[i].append(0)

    for i in range(size):
        matrix[i][i] = 1
    return matrix


def swapRows(matrix, row1, row2, size):
    identity = identityMatrix(size)
    identity[row1], identity[row2] = identity[row2], identity[row1]
    return multiplyMatrices(identity, matrix, size)


def editMatrix(matrix, x, y, val):
    matrix[x][y] = val
    return matrix

# assuming nxn+1 matrix
def multiplyMatrices(matrix1, matrix2, size):
    mat = []
    for i in range(size+1):
        mat.append([])

    for i in range(size):
        for j in range(size + 1):
            sum = 0
            for k in range(size):
                sum += matrix1[i][k] * matrix2[k][j]
            mat[i].append(sum)
    return mat


def nullify(matrix, x, y, size, pivot):
    identity = identityMatrix(size)
    return multiplyMatrices(editMatrix(identity, x, y, -1*matrix[x][y]/pivot), matrix, size)


def setPivotToOne(matrix, x, y, size):
    identity = identityMatrix(size)
    return multiplyMatrices(editMatrix(identity, x, y, 1/matrix[x][y]), matrix, size)


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def cubic_spline_interpolation(table, target):
    gamma = []
    mu = []
    d = []
    h = []

    for i in range(0, len(table) - 1):
        h.append(table[i + 1].x - table[i].x)

    for i in range(1, len(table) - 1):
        gamma.append(h[i]/(h[i] + h[i - 1]))
        mu.append(1 - h[i]/(h[i] + h[i - 1]))
        d.append((6/(h[i] + h[i - 1])*((table[i + 1].y - table[i].y)/h[i] - (table[i].y - table[i - 1].y)/h[i - 1])))

    # build matrix
    mat = identityMatrix(len(d))
    for i in range(len(d)):
        mat[i][i] = 2
        if i != 0:
            mat[i][i - 1] = mu[i]
        if i != len(d) - 1:
            mat[i][i + 1] = gamma[i]
        mat[i].append(d[i])

    # extract result
    m = [0]
    res = solveMatrix(mat, len(d))
    for x in range(len(res) - 1):
        m.append(res[x][-1])
    m.append(0)

    for y in range(len(table) - 1):
        if target > table[y].x:
            if target < table[y + 1].x:
                return ((table[y + 1].x - target)**3 * m[y] + (target - table[y].x) ** 3 * m[y + 1])/(6*h[y]) \
                        + ((table[y + 1].x - target) * table[y].y + (target - table[y].x) * table[y + 1].y)/(h[y]) \
                        - h[y]*((table[y + 1].x - target) * m[y] + (target - table[y].x) * m[y + 1])/6
    print("Target out of bounds")

def lagrange(X, Y, xr):
    n = len(X)
    Lmult = 1 #calculating multipliper in L calculation
    Lsum = 0 # calculating sum in P calculation
    for i in range(n):
        Lmult=1
        for j in range(n):
            if i != j:
                Lmult *= ((xr - X[j]) / (X[i] - X[j]))

        Lsum += Lmult*Y[i]

    return Lsum

def neville(datax, datay, x):

    n = len(datax)
    p = n*[0]
    for k in range(n):
        for i in range(n-k):
            if k == 0:
                p[i] = datay[i]
            else:
                p[i] = ((x-datax[i+k])*p[i] + (datax[i]-x)*p[i+1]) / (datax[i]-datax[i+k])
    return p[0]

def Linear_Interpolation (table,point_Value):
    x1, x2, y1, y2 = 0, 0, 0, 0
    flag = False

    sum = 0
    x1, x2, y1, y2 = 0, 0, 0, 0
    length = len(table)
    for key in table:
        if key < point_Value:
            x1 = key
            y1 = table[key]
        if key > point_Value:
            x2 = key
            y2 = table[key]
            flag = True  # if a bigger point found
            break
        if key == point_Value:
            return table[key]

    try:
        z = (((y1 - y2) / (x1 - x2)) * point_Value) + (((y2 * x1) - (y1 * x2)) / (x1 - x2))
    except ZeroDivisionError:
        z = 0


    if not flag:
        return "the x is out of the boundaries of the given table"
    return z

def Polynomial_Interpolation (table, xf):
    y_vactor = []
    x_vector = []
    result_vector = []
    for key in table:
        x_vector.append(key)
        y_vactor.append([table[key]])

    result_mat = CalcMatrix(create_matrix(x_vector), y_vactor)
    for col in result_mat:
        result_vector.append(col[0])

    return create_polinom(result_vector)(xf)

def create_polinom(param_list):
    def polinom(x):
        pow = 0
        val = 0
        for param in param_list:
            val += param * math.pow(x, pow)
            pow += 1
        return val
    return polinom


table = {}
point_table = [Point(0, 0)]
x_table = []
y_table = []

inp1 = input("Enter new X: ")
while inp1 != "":
    inp2 = float(input("Enter new Y: "))
    inp1 = float(inp1)
    table[inp1] = inp2
    x_table.append(inp1)
    y_table.append(inp2)
    point_table.append(Point(inp1, inp2))
    print("New point ({}, {}) created.".format(inp1, inp2))
    inp1 = input("Enter new X: (Press ENTER to exit)")
target_point = float(input("Enter target point: "))

print("Linear interpolation:", Linear_Interpolation(table, target_point))
print("Polynomial interpolation :", Polynomial_Interpolation (table, target_point))
print("Lagrange interpolation:", lagrange(x_table, y_table, target_point))
print("Neville interpolation:", neville(x_table, y_table, target_point))
print("Cubic spline interpolation:", cubic_spline_interpolation(point_table, target_point))
