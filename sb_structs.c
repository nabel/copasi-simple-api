#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "sb_structs.h"

SBAPIEXPORT sb_matrix sb_createMatrix(int rows, int cols)
{
	int i;
	sb_matrix M;
	M.rows = rows;
	M.cols = cols;
	M.colnames = sb_createStringsArray(cols);
	M.rownames = sb_createStringsArray(rows);
	
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

SBAPIEXPORT sb_table sb_createTable(int rows, int cols)
{
	int i;
	sb_table M;
	M.rows = rows;
	M.cols = cols;

	M.colnames = sb_createStringsArray(cols);
	M.rownames = sb_createStringsArray(rows);
	if (rows > 0 && cols > 0)
	{
		M.strings = (char**)malloc( rows * cols * sizeof(char*) );
		for (i=0; i < (rows*cols); ++i)
			M.strings[i] = 0;
	}
	else
		M.strings = 0;
	return M;
}

SBAPIEXPORT sb_strings sb_createStringsArray(int len)
{
	int i;
	sb_strings A;
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

SBAPIEXPORT sb_items sb_createItemsArray(int len)
{
	int i;
	sb_items A;
	if (len < 1)
	{
		A.length = 0;
		A.items = 0;
	}
	else
	{
		A.length = len;
		A.items = (long*)malloc(len * sizeof(long));
		for (i=0; i < len; ++i)
			A.items[i] = 0;
	}
	return A;
}

SBAPIEXPORT double sb_getMatrixValue(sb_matrix M, int i, int j)
{ 
	if (M.values && i >= 0 && j >= 0 && i < M.rows && j < M.cols)
		return M.values[ i*M.cols + j ];
	return 0.0;
}

SBAPIEXPORT void sb_setMatrixValue(sb_matrix M, int i, int j, double d)
{ 
	if (M.values && i >= 0 && j >= 0 && i < M.rows && j < M.cols)
		M.values[ i*M.cols + j ] = d;
}

SBAPIEXPORT const char * sb_getRowName(sb_matrix M, int i)
{ 
	return sb_getString(M.rownames,i);
}

SBAPIEXPORT void sb_setRowName(sb_matrix M, int i, const char * s)
{
	sb_setString(M.rownames,i,s);
}

SBAPIEXPORT const char * sb_getColumnName(sb_matrix M, int i)
{ 
	return sb_getString(M.colnames,i);
}

SBAPIEXPORT void sb_setColumnName(sb_matrix M, int i, const char * s)
{
	sb_setString(M.colnames,i,s);
}

SBAPIEXPORT const char* sb_getTableValue(sb_table S, int i, int j)
{
	if (S.strings && i >= 0 && j >= 0 && i < S.rows && j < S.cols)
		return S.strings[ i*S.cols + j ];
	return 0;
}

SBAPIEXPORT void sb_setTableValue(sb_table S, int i, int j, const char * s)
{
	int n=0;
	char * str;
	if (i >= 0 && j >= 0 && i < S.rows && j < S.cols)
	{
		while (s && s[n]) ++n;
		str = (char*)malloc((n+1)*sizeof(char));
		sprintf(str,"%s\0",s);
	
		S.strings[ i*S.cols + j ] = str;
	}
}

SBAPIEXPORT const char* sb_getString(sb_strings S, int i)
{
	if (S.strings && i >= 0 && i < S.length)
		return S.strings[ i ];
	return 0;
}

SBAPIEXPORT void sb_setString(sb_strings S, int i, const char * s)
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

SBAPIEXPORT long sb_getItem(sb_items A, int i)
{
	if (i >= 0 && i < A.length)
		return A.items[ i ];
	return 0;
}

SBAPIEXPORT void sb_setItem(sb_items A, int i, long o)
{
	if (i >= 0 && i < A.length)
		A.items[ i ] = o;
}

SBAPIEXPORT void sb_deleteMatrix(sb_matrix M)
{
	if (M.values)
		free(M.values);
	M.rows = M.cols = 0;	
	M.values = 0;
	sb_deleteStringsArray(M.rownames);
	sb_deleteStringsArray(M.colnames);
}

SBAPIEXPORT void sb_deleteTable(sb_table M)
{
	if (M.strings)
		free(M.strings);
	M.rows = M.cols = 0;
	M.strings = 0;
	sb_deleteStringsArray(M.rownames);
	sb_deleteStringsArray(M.colnames);
}

SBAPIEXPORT void sb_deleteItemsArray(sb_items A)
{
	if (A.items) 
		free(A.items);
	A.length = 0;
	A.items = 0;
}

SBAPIEXPORT void sb_deleteStringsArray(sb_strings C)
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

SBAPIEXPORT sb_matrix sb_appendColumns(sb_matrix A, sb_matrix B)
{
	int i,j,k=0;
	sb_matrix C;
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

SBAPIEXPORT sb_matrix sb_appendRows(sb_matrix A, sb_matrix B)
{
	int i,j,k=0;
	sb_matrix C;
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

SBAPIEXPORT void sb_printMatrixToFile(const char* s, sb_matrix output)
{
	int i,j;
	FILE * outfile = fopen(s,"w+");
	if (output.colnames.strings)
	{
		fprintf(outfile, "#\t");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				fprintf(outfile, "%s\t", sb_getColumnName(output, j));
			else
				fprintf(outfile, "%s\n", sb_getColumnName(output, j));
	}
	for (i=0; i < output.rows; ++i)
	{
		if (sb_getRowName(output,i))
			fprintf(outfile, "%s\t", sb_getRowName(output, i));
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				fprintf(outfile, "%lf\t", sb_getMatrixValue(output, i, j));
			else
				fprintf(outfile, "%lf\n", sb_getMatrixValue(output, i, j));
	}
	
	fclose(outfile);
}

SBAPIEXPORT void sb_printOutMatrix(sb_matrix output)
{
	int i,j;
	if (output.colnames.strings)
	{
		printf("\t");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				printf("%s\t", sb_getColumnName(output, j));
			else
				printf("%s\n", sb_getColumnName(output, j));
	}

	for (i=0; i < output.rows; ++i)
	{
		if (sb_getRowName(output,i))
			printf("%s\t", sb_getRowName(output, i));
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				printf("%lf\t", sb_getMatrixValue(output, i, j));
			else
				printf("%lf\n", sb_getMatrixValue(output, i, j));
	}
}


SBAPIEXPORT void sb_printTableToFile(const char* s, sb_table output)
{
	int i,j;
	FILE * outfile = fopen(s,"w+");
	if (output.colnames.strings)
	{
		fprintf(outfile, "#\t");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				fprintf(outfile, "%s\t", sb_getString(output.colnames, j));
			else
				fprintf(outfile, "%s\n", sb_getString(output.colnames, j));
	}

	if (output.strings)
	{
		for (i=0; i < output.rows; ++i)
		{
			if (sb_getString(output.rownames,i))
				fprintf(outfile, "%s\t", sb_getString(output.rownames,i));
			for (j=0; j < output.cols; ++j)
				if (j < (output.cols-1))
					fprintf(outfile, "%s\t", sb_getTableValue(output, i, j));
				else
					fprintf(outfile, "%s\n", sb_getTableValue(output, i, j));
		}
	}
	fclose(outfile);
}

SBAPIEXPORT void sb_printOutTable(sb_table output)
{
	int i,j;
	if (output.colnames.strings)
	{
		printf("\t");
		for (j=0; j < output.cols; ++j)
			if (j < (output.cols-1))
				printf("%s\t", sb_getString(output.colnames, j));
			else
				printf("%s\n", sb_getString(output.colnames, j));
	}
	
	if (output.strings)
	{
		for (i=0; i < output.rows; ++i)
		{
			if (sb_getString(output.rownames,i))
				printf("%s\t", sb_getString(output.rownames,i));
			for (j=0; j < output.cols; ++j)
				if (j < (output.cols-1))
					printf("%s\t", sb_getTableValue(output, i, j));
				else
					printf("%s\n", sb_getTableValue(output, i, j));
		}
	}
}

SBAPIEXPORT 
int sb_getStringIndex(sb_strings A, const char * s)
{
	int i=0;
	if (A.length == 0 || !A.strings) return -1;

	for (i=0; i < A.length; ++i)
		if (sb_getString(A,i) && strcmp( sb_getString(A,i) , s ) == 0)
			return i;
	return -1;
}

SBAPIEXPORT 
int sb_getRowIndex(sb_matrix m, const char * s)
{
	return sb_getStringIndex( m.rownames, s );
}

SBAPIEXPORT 
int sb_getColumnIndex(sb_matrix m, const char * s)
{
	return sb_getStringIndex( m.colnames, s );
}

