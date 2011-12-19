from copasi import *

def toItems(array):
    n = len(array);
    A = c_createItemsArray(n);
    for i in range(0, n):
        c_setItem(A, i, array[i]);

    return A;

def fromItems(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = c_getItem(array,i);

    #c_deleteItemsArray(array);
    return A;

def toStrings(array):
    n = len(array);
    A = c_createStringsArray(n);
    for i in range(0, n):
        c_setString(A, i, array[i]);

    return A;

def fromStrings(array):
    n = array.length;
    A = range(0,n);
    for i in range(0, n):
        A[i] = c_getString(array,i);

    #c_deleteStringsArray(array);
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
                A[i][j] = c_getMatrixValue(matrix,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = c_getMatrixValue(matrix,j,i);

    #c_deleteMatrix(matrix);
    return [rows, cols, A];

def toMatrix(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    A = c_createMatrix(0,0);
    if row_wise:
        A = c_createMatrix(n,m);
    else:
        A = c_createMatrix(m,n);
    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                c_setMatrixValue(A,i,j,lists[i][j]);
            else:
                c_setMatrixValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        c_setRowName(A,i,rows[i]);

    for i in range(0,m):
        c_setColumnName(A,i,cols[i]);

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
                A[i][j] = c_getTableValue(table,i,j);
    else:
        A = range(0,m);
        for i in range(0, m):
            A[i] = range(0,n);
            for j in range(0,n):
                A[i][j] = c_getTableValue(table,j,i);

    return [rows, cols, A];

def toTable(lists, row_wise = False , rows = [], cols = []):
    n = len(lists);
    m = len(lists[0]);
    
    A = c_createTable(0,0);
    if row_wise:
        A = c_createTable(n,m);
    else:
        A = c_createTable(m,n);

    for i in range(0, n):
        for j in range(0,m):
            if row_wise:
                c_setTableValue(A,i,j,lists[i][j]);
            else:
                c_setTableValue(A,j,i,lists[i][j]);
    n = len(rows);
    m = len(cols);

    for i in range(0,n):
        c_setString(A.rownames,i,rows[i]);

    for i in range(0,m):
        c_setString(A.colnames,i,cols[i]);

    return A;

