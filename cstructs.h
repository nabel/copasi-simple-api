 /**
  * @file    C_structs.h
  * @brief   Additional structures and methods used with the Simple C API

 */

#ifndef TINKERCELL_CSTRUCTS_H
#define TINKERCELL_CSTRUCTS_H

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#        define BEGIN_C_DECLS extern "C" {
#        define END_C_DECLS }
#   else
#        define BEGIN_C_DECLS
#        define END_C_DECLS
#endif
#endif

# ifndef CAPIEXPORT
#  if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#    if defined(STATIC_LINKED)
#          define CAPIEXPORT
#    else
#   if defined(C_EXPORTS) || defined(tinkercellapi_EXPORTS)
#              if defined(USE_STDCALL)
#                   define CAPIEXPORT __stdcall __declspec(dllexport)
#              else
#                   define CAPIEXPORT __declspec(dllexport)
#              endif
#          else
#              define CAPIEXPORT
#          endif
#     endif
#  else
#    if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#      define CAPIEXPORT __attribute__ ((visibility("default")))
#    else
#      define CAPIEXPORT
#    endif
#  endif
# endif //CAPIEXPORT

BEGIN_C_DECLS

/*!\brief An array of strings with length information. Use c_getString(M,i) to get the i-th string*/
typedef struct
{
	int length;
	char ** strings;
} c_strings;

/*!\brief A 2D table of doubles with row and column names. 
Use c_getMatrixValue(M,i,j) to get the i,j-th value in c_matrix M

Use c_setMatrixValue(M,i,j,value) to set a value */
typedef struct
{
	int rows;            /*!< Number of rows in the matrix */
	int cols;            /*!< Number of columns in the matrix */
	double * values;     /*!< Pointer to the values contained in the matrix */
	c_strings rownames; /*!< Pointer to list of row labels */
	c_strings colnames; /*!< Pointer to list of column labels */
} c_matrix;


/*!\brief Create a matrix with the given number of rows and columns
 \param int number of rows
 \param int number of columns
 \return c_matrix
 \ingroup Basic
*/
CAPIEXPORT c_matrix c_createMatrix(int rows, int cols);


/*!\brief Create an array of strings
 \param int length
 \return c_strings
 \ingroup Basic
*/
CAPIEXPORT c_strings c_createStringsArray(int len);


/*!\brief Get i,jth value from a c_matrix
 \param c_matrix matrix
 \param int row
 \param int column
 \return double value at the given row, column
 
 Use c_setMatrixValue(M,i,j,value) to set a value

 \ingroup Basic
*/
CAPIEXPORT double c_getMatrixValue(c_matrix M, int i, int j);


/*!\brief Set i,jth value of a c_matrix
 \param c_matrix matrix
 \param int row
 \param int column
 \param double value at the given row, column
 \ingroup Basic
*/
CAPIEXPORT void c_setMatrixValue(c_matrix M, int i, int j, double d);


/*!\brief Get ith row name from a c_matrix
 \param c_matrix matrix
 \param int row
 \return string row name
 \ingroup Basic
*/
CAPIEXPORT const char * c_getRowName(c_matrix M, int i);


/*!\brief Set ith row name for a c_matrix
 \param c_matrix matrix
 \param int row
 \param string row name
 \ingroup Basic
*/
CAPIEXPORT void c_setRowName(c_matrix M, int i, const char * s);


/*!\brief Get jth column name of a c_matrix
 \param c_matrix matrix
 \param int column
 \return string column name
 \ingroup Basic
*/
CAPIEXPORT const char * c_getColumnName(c_matrix M, int j);


/*!\brief Set jth column name of a c_matrix
 \param c_matrix matrix
 \param int column
 \param string column name
 \ingroup Basic
*/
CAPIEXPORT void c_setColumnName(c_matrix M, int j, const char * s);

/*!\brief Get ith string in array of strings
 \param c_strings array
 \param int index
 \return string value
 \ingroup Basic
*/
CAPIEXPORT const char* c_getString(c_strings S, int i);


/*!\brief Set ith string in array of strings
 \param c_strings array
 \param int index
 \param string value
 \ingroup Basic
*/
CAPIEXPORT void c_setString(c_strings S, int i, const char * c);

/*!\brief Get the index of a string in the array
 \param c_strings array
 \param char* a string in the array
 \return int index of that string
 \ingroup Basic
*/
CAPIEXPORT int c_getStringIndex(c_strings A, const char * s);


/*!\brief get the row number of a row name
 \param c_matrix matrix
 \param char* a string in the matrix
 \return int index of that string
 \ingroup Basic
*/
CAPIEXPORT int c_getRowIndex(c_matrix, const char * s);


/*!\brief Get the column number of a column name
 \param c_matrix matrix
 \param char* a string in the matrix
 \return int index of that string
 \ingroup Basic
*/
CAPIEXPORT int c_getColumnIndex(c_matrix, const char * s);


/*!\brief Delete a matrix
 \param &c_matrix pointer to matrix
 \ingroup Basic
*/
CAPIEXPORT void c_deleteMatrix(c_matrix M);


/*!\brief Delete an array of strings
 \param &c_strings pointer to array
 \ingroup Basic
*/
CAPIEXPORT void c_deleteStringsArray(c_strings C);


/*!\brief Combine two matrices by appending their columns. Row size must be equal for both matrices
 \param c_matrix first matrix
 \param c_matrix fsecond matrix
 \return c_matrix new combined matrix
 \ingroup Basic
*/
CAPIEXPORT c_matrix c_appendColumns(c_matrix A, c_matrix B);


/*!\brief Combine two matrices by appending their rows. Column sizes must be equal for both matrices
 \param c_matrix first matrix
 \param c_matrix fsecond matrix
 \return c_matrix new combined matrix
 \ingroup Basic
*/
CAPIEXPORT c_matrix c_appendRows(c_matrix A, c_matrix B);


/*!\brief Print a matrix to file
 \param char* file name
 \param c_matrix
 \ingroup Basic
*/
CAPIEXPORT void c_printMatrixToFile(const char* file, c_matrix M);


/*!\brief Print a matrix to stdout
 \param char* file name
 \param c_matrix
 \ingroup Basic
*/
CAPIEXPORT void c_printOutMatrix(c_matrix M);


END_C_DECLS
#endif

