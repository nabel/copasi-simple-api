#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cstructs.h"

CAPIEXPORT c_matrix c_createMatrix(int rows, int cols)
{
	int i;
	c_matrix M;
	M.rows = rows;
	M.cols = cols;
	M.colnames = c_createStringsArray(cols);
	M.rownames = c_createStringsArray(rows);
	
	if (rows > 0 && cols > 0)
	{
		M.values = (double*)malloc( rows * cols * sizeof(double) );
		for (i=0; i < (rows*cols); ++i)
			M.values[i] = 0.0;
	}
	else
		M.values = 0;
	return M;
}

CAPIEXPORT c_strings c_createStringsArray(int len)
{
	int i;
	c_strings A;
	if (len < 1)
	{
		A.length = 0;
		A.strings = 0;
	}
	else
	{
		A.length = len;
		A.strings = (char**)malloc(len * sizeof(char*));
		for (i=0; i < len; ++i)
			A.strings[i] = 0;
	}
	return A;
}

CAPIEXPORT double c_getMatrixValue(c_matrix M, int i, int j)
{ 
	if (M.values && i >= 0 && j >= 0 && i < M.rows && j < M.cols)
		return M.values[ i*M.cols + j ];
	return 0.0;
}

CAPIEXPORT void c_setMatrixValue(c_matrix M, int i, int j, double d)
{ 
	if (M.values && i >= 0 && j >= 0 && i < M.rows && j < M.cols)
		M.values[ i*M.cols + j ] = d;
}

CAPIEXPORT const char * c_getRowName(c_matrix M, int i)
{ 
	return c_getString(M.rownames,i);
}

CAPIEXPORT void c_setRowName(c_matrix M, int i, const char * s)
{
	c_setString(M.rownames,i,s);
}

CAPIEXPORT const char * c_getColumnName(c_matrix M, int i)
{ 
	return c_getString(M.colnames,i);
}

CAPIEXPORT void c_setColumnName(c_matrix M, int i, const char * s)
{
	c_setString(M.colnames,i,s);
}

CAPIEXPORT const char* c_getString(c_strings S, int i)
{
	if (S.strings && i >= 0 && i < S.length)
		return S.strings[ i ];
	return 0;
}

CAPIEXPORT void c_setString(c_strings S, int i, const char * s)
{
	int n=0;
	char * str;
	if (i >= 0 && i < S.length)
	{
		while (s && s[n]) ++n;
		
		if (n > 0)
		{
			str = (char*)malloc((n+1)*sizeof(char));
			sprintf(str,"%s\0",s);
			if (S.strings[i])
				free(S.strings[i]);
			S.strings[ i ] = str;
		}
	}
}

CAPIEXPORT void c_deleteMatrix(c_matrix M)
{
	if (M.values)
		free(M.values);
	M.rows = M.cols = 0;	
	M.values = 0;
	c_deleteStringsArray(M.rownames);
	c_deleteStringsArray(M.colnames);
}

CAPIEXPORT void c_deleteStringsArray(c_strings C)
{
	int i;
	if (C.strings)
	{
		for (i=0; i < C.length; ++i) 
			if (C.strings[i]) 
				free(C.strings[i]);
		free(C.strings);
	}
	C.length = 0;
	C.strings = 0;
}

CAPIEXPORT c_matrix c_appendColumns(c_matrix A, c_matrix B)
{
	int i,j,k=0;
	c_matrix C;
	int fromA = 0, toA = A.cols, fromB = 0, toB = B.cols;

	C.colnames.length = C.rownames.length = 0;
	C.colnames.strings = C.rownames.strings = 0;
	C.rows = C.cols = 0;
	C.values = 0;

	if (A.rows != B.rows) return C;
	if (fromA < 0 || toA < 0 || fromA > A.cols || toA > A.cols ||
		fromB < 0 || toB < 0 || fromB > B.cols || toB > B.cols ||
		fromA >= toA || fromB >= toB)
		return C;

	C.rows = A.rows;
	C.cols = ((toA - fromA) + (toB - fromB));
	C.values = (double*)malloc( C.rows * C.cols * sizeof(double) );
	
	for (i=0; i < (C.rows*C.cols); ++i)
		C.values[i] = 0.0;

	if (A.colnames.strings && B.colnames.strings)
	{
		C.colnames.length = C.cols;
		C.colnames.strings = (char**)malloc( C.cols * sizeof(char*) );
		for (i=0; i < A.cols; ++i)
		{
			k = 0;
			while (A.colnames.strings[i] && A.colnames.strings[i][k]) ++k;
			C.colnames.strings[i] = (char*)malloc((1+k) * sizeof(char));
			C.colnames.strings[i][k] = 0;
			for (j=0; j < k; ++j)
				C.colnames.strings[i][j] = A.colnames.strings[i][j];
		}
		for (i=0; i < B.cols; ++i)
		{
			k = 0;
			while (B.colnames.strings[i] && B.colnames.strings[i][k]) ++k;
			C.colnames.strings[i+A.cols] = (char*)malloc((1+k) * sizeof(char));
			C.colnames.strings[i+A.cols][k] = 0;
			for (j=0; j < k; ++j)
				C.colnames.strings[i+A.cols][j] = B.colnames.strings[i][j];
		}
	}

	if (A.rownames.strings && B.rownames.strings)
	{
		C.rownames.length = C.rows;
		C.rownames.strings = (char**)malloc( C.rows * sizeof(char*) );
		for (i=0; i < A.rows; ++i)
		{
			k = 0;
			while (A.rownames.strings[i] && A.rownames.strings[i][k]) ++k;
			C.rownames.strings[i] = (char*)malloc((1+k) * sizeof(char));
			C.rownames.strings[i][k] = 0;
			for (j=0; j < k; ++j)
				C.rownames.strings[i][j] = A.rownames.strings[i][j];
		}
	}
	
	k = (toA - fromA);
	for (i=fromA; i < toA; ++i)
	{
		for (j=0; j < C.rows; ++j)
			C.values[ j*C.cols + i - fromA ] = A.values[ j * A.cols + i ];
	}
	
	for (i=fromB; i < toB; ++i)
	{
		for (j=0; j < C.rows; ++j)
			C.values[ j*C.cols + k + i - fromB ] = B.values[ j * B.cols + i ];
	}
	return C;
}

CAPIEXPORT c_matrix c_appendRows(c_matrix A, c_matrix B)
{
	int i,j,k=0;
	c_matrix C;
	int fromA = 0, toA = A.rows, fromB = 0, toB = B.rows;

	C.colnames.strings = C.rownames.strings = 0;
	C.colnames.length = C.rownames.length = 0;
	C.rows = C.cols = 0;
	C.values = 0;

	if (A.cols != B.cols) return C;
	if (fromA < 0 || toA < 0 || fromA > A.cols || toA > A.cols ||
		fromB < 0 || toB < 0 || fromB > B.cols || toB > B.cols ||
		fromA >= toA || fromB >= toB)
		return C;

	C.cols = A.cols;
	C.rows = ((toA - fromA) + (toB - fromB));
	C.values = (double*)malloc( C.rows * C.cols * sizeof(double) );
	
	for (i=0; i < (C.rows*C.cols); ++i)
		C.values[i] = 0.0;

	if (A.rownames.strings && B.rownames.strings)
	{
		C.rownames.length = C.rows;
		C.rownames.strings = (char**)malloc( C.rows * sizeof(char*) );
		for (i=0; i < A.rows; ++i)
		{
			k = 0;
			while (A.rownames.strings[i] && A.rownames.strings[i][k]) ++k;
			C.rownames.strings[i] = (char*)malloc((k+1) * sizeof(char));
			C.rownames.strings[i][k] = 0;
			for (j=0; j < k; ++j)
				C.rownames.strings[i][j] = A.colnames.strings[i][j];
		}
		for (i=0; i < B.rows; ++i)
		{
			k = 0;
			while (B.rownames.strings[i] && B.rownames.strings[i][k]) ++k;
			C.rownames.strings[i+A.cols] = (char*)malloc((1+k) * sizeof(char));
			C.rownames.strings[i+A.cols][k] = 0;
			for (j=0; j < k; ++j)
				C.rownames.strings[i+A.cols][j] = B.rownames.strings[i][j];
		}
	}



	if (A.colnames.strings && B.colnames.strings)
	{
		C.colnames.length = C.cols;
		C.colnames.strings = (char**)malloc( (C.cols + 1) * sizeof(char*) );
		C.colnames.strings[C.cols] = 0;
		for (i=0; i < A.cols; ++i)
		{
			k = 0;
			while (A.colnames.strings[i] && A.colnames.strings[i][k]) ++k;
			C.colnames.strings[i] = (char*)malloc((1+k) * sizeof(char));
			C.colnames.strings[i][k] = 0;
			for (j=0; j < k; ++j)
				C.colnames.strings[i][j] = A.colnames.strings[i][j];
		}
	}
	
	k = (toA - fromA);
	for (i=fromA; i < toA; ++i)
	{
		for (j=0; j < C.cols; ++j)
			C.values[ (i-fromA)*C.cols + j ] = A.values[ i * A.cols + j ];
	}
	
	for (i=fromB; i < toB; ++i)
	{
		for (j=0; j < C.cols; ++j)
			C.values[ (i + k - fromB)*C.cols + j ] = B.values[ i * B.cols + j ];
	}

	return C;
}

CAPIEXPORT void c_printMatrixToFile(const char* s, c_matrix output)
{
	int i,j;
	FILE * outfile = fopen(s,"w+");
	if (output.colnames.strings)
	{
		fprintf(outfile, "#");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				fprintf(outfile, "%s\t", c_getColumnName(output, j));
			else
				fprintf(outfile, "%s\n", c_getColumnName(output, j));
	}
	for (i=0; i < output.rows; ++i)
	{
		if (c_getRowName(output,i))
			fprintf(outfile, "%s\t", c_getRowName(output, i));
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				fprintf(outfile, "%lf\t", c_getMatrixValue(output, i, j));
			else
				fprintf(outfile, "%lf\n", c_getMatrixValue(output, i, j));
	}
	
	fclose(outfile);
}

CAPIEXPORT void c_printOutMatrix(c_matrix output)
{
	int i,j;
	if (output.colnames.strings)
	{
		printf("\t");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				if (c_getColumnName(output, j))
					printf("%s\t", c_getColumnName(output, j));
				else
					printf("\t");
			else
				if (c_getColumnName(output, j))
					printf("%s\n", c_getColumnName(output, j));
				else
					printf("\n");
	}

	for (i=0; i < output.rows; ++i)
	{
		if (c_getRowName(output,i))
			printf("%s\t", c_getRowName(output, i));
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				printf("%lf\t", c_getMatrixValue(output, i, j));
			else
				printf("%lf\n", c_getMatrixValue(output, i, j));
	}
}

CAPIEXPORT 
int c_getStringIndex(c_strings A, const char * s)
{
	int i=0;
	if (A.length == 0 || !A.strings) return -1;

	for (i=0; i < A.length; ++i)
		if (c_getString(A,i) && strcmp( c_getString(A,i) , s ) == 0)
			return i;
	return -1;
}

CAPIEXPORT 
int c_getRowIndex(c_matrix m, const char * s)
{
	return c_getStringIndex( m.rownames, s );
}

CAPIEXPORT 
int c_getColumnIndex(c_matrix m, const char * s)
{
	return c_getStringIndex( m.colnames, s );
}

