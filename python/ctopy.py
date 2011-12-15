import copasi

def toItems(array):
    n = len(array);
    A = copasi.tc_createItemsArray(n);
    for i in range(0, n):
        copasi.tc_setItem(A, i, array[i]);

    return A;

def fromItems(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = copasi.tc_getItem(array,i);

    #copasi.tc_deleteItemsArray(array);
    return A;

def toStrings(array):
    n = len(array);
    A = copasi.tc_createStringsArray(n);
    for i in range(0, n):
        copasi.tc_setString(A, i, array[i]);

    return A;

def fromStrings(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = copasi.tc_getString(array,i);

    #copasi.tc_deleteStringsArray(array);
    return A;

def fromMatrix(matrix, row_wise = False):
    n = matrix.rows;
    m = matrix.cols;
    cols = fromStrings(matrix.colnames);
    rows = fromStrings(matrix.rownames);
    if row_wise:
        A = range(0,n);
        for i in range(0, n):
            A[i] = range(0,m);
            for j in range(0,m):
                A[i][j] = copasi.tc_getMatrixValue(matrix,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = copasi.tc_getMatrixValue(matrix,j,i);

    #copasi.tc_deleteMatrix(matrix);
    return [rows, cols, A];

def toMatrix(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    A = copasi.tc_createMatrix(0,0);
    if row_wise:
        A = copasi.tc_createMatrix(n,m);
    else:
        A = copasi.tc_createMatrix(m,n);
    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                copasi.tc_setMatrixValue(A,i,j,lists[i][j]);
            else:
                copasi.tc_setMatrixValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        copasi.tc_setRowName(A,i,rows[i]);

    for i in range(0,m):
        copasi.tc_setColumnName(A,i,cols[i]);

    return A;

def fromTable(table, row_wise = False):
    n = table.rows;
    m = table.cols;
    cols = fromStrings(table.colnames);
    rows = fromStrings(table.rownames);
    if row_wise:
        A = range(0,n);
        for i in range(0, n):
            A[i] = range(0,m);
            for j in range(0,m):
                A[i][j] = copasi.tc_getTableValue(table,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = copasi.tc_getTableValue(table,j,i);

    return [rows, cols, A];

def toTable(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    
    A = copasi.tc_createTable(0,0);
    if row_wise:
        A = copasi.tc_createTable(n,m);
    else:
        A = copasi.tc_createTable(m,n);

    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                copasi.tc_setTableValue(A,i,j,lists[i][j]);
            else:
                copasi.tc_setTableValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        copasi.tc_setString(A.rownames,i,rows[i]);

    for i in range(0,m):
        copasi.tc_setString(A.colnames,i,cols[i]);

    return A;

def toHex(r,g,b):
    hexchars = "0123456789ABCDEF0";
    return "#" + hexchars[r / 16] + hexchars[r % 16] + hexchars[g / 16] + hexchars[g % 16] + hexchars[b / 16] + hexchars[b % 16];

def fromC(x):
    if type(x) == type(copasi.tc_createMatrix(0,0)): return fromMatrix(x)
    if type(x) == type(copasi.tc_createStringsArray(0)): return fromStrings(x)
    if type(x) == type(copasi.tc_createItemsArray(0)): return fromItems(x)
    if type(x) == type(copasi.tc_createTable(0,0)): return fromTable(x)
    return x

def toC(x, rows = [], cols = []):
    if type(x) == type([]) and len(x) > 0 and type(x[0]) == type([]):
        if (type(x[0][0]) == type(1.0)):
            return toMatrix(x,False, rows,cols)
        elif (type(x[0][0]) == type('hello')):
            return toTable(x,False, rows,cols)
    if type(x) == type([]) and len(x) > 0 and type(x[0]) == type('hello'): return toStrings(x)
    if type(x) == type([]) and len(x) > 0 and type(x[0]) == type(1): return toItems(x)
    return x

