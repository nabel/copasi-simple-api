 /**
  * @file    SBStructs.h
  * @brief   Additional structures and methods used with the SBW C API

 */

#ifndef SB_CSTRUCTS_H
#define SB_CSTRUCTS_H

#ifndef BEGIN_C_DECLS
#ifdef __cplusplus
#        define BEGIN_C_DECLS extern "C" {
#        define END_C_DECLS }
#   else
#        define BEGIN_C_DECLS
#        define END_C_DECLS
#endif
#endif

# ifndef SBAPIEXPORT
#  if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#    if defined(STATIC_LINKED)
#          define SBAPIEXPORT
#    else
#   if defined(SB_EXPORTS) || defined(sbapi_EXPORTS)
#              if defined(USE_STDCALL)
#                   define SBAPIEXPORT __stdcall __declspec(dllexport)
#              else
#                   define SBAPIEXPORT __declspec(dllexport)
#              endif
#          else
#              define SBAPIEXPORT
#          endif
#     endif
#  else
#    if defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#      define SBAPIEXPORT __attribute__ ((visibility("default")))
#    else
#      define SBAPIEXPORT
#    endif
#  endif
# endif // SBAPIEXPORT

BEGIN_C_DECLS

/*!\brief An array of strings with length information. Use sb_getString(M,i) to get the i-th string*/
typedef struct
{
	int length;
	char ** strings;
} sb_strings;


/*!\brief An array of int objects with length information. Use sb_getItem(M,i) to get the i-th item.*/
typedef struct
{
	int length;
	long* items;
} sb_items;


/*!\brief A 2D table of doubles with row and column names. 
Use sb_getMatrixValue(M,i,j) to get the i,j-th value in sb_matrix M

Use sb_setMatrixValue(M,i,j,value) to set a value */
typedef struct
{
	int rows;            /*!< Number of rows in the matrix */
	int cols;            /*!< Number of columns in the matrix */
	double * values;     /*!< Pointer to the values contained in the matrix */
	sb_strings rownames; /*!< Pointer to list of row labels */
	sb_strings colnames; /*!< Pointer to list of column labels */
} sb_matrix;


/*!\brief A 2D table of strings with row and column names. Use sb_getTableValue(M,i,j) to get the i,j-th value in sb_matrix M.*/
typedef struct
{
	int rows, cols;
	char ** strings;
	sb_strings rownames;  /*!< Pointer to list of column labels */
	sb_strings colnames;  /*!< Pointer to list of column labels */
} sb_Table;


/*!\brief Create a matrix with the given number of rows and columns
 \param int number of rows
 \param int number of columns
 \return sb_matrix
 \ingroup Basic
*/
SBAPIEXPORT sb_matrix sb_createMatrix(int rows, int cols);


/*!\brief Create a strings table with the given number of rows and columns
 \param int number of rows
 \param int number of columns
 \return sb_table
 \ingroup Basic
*/
SBAPIEXPORT sb_table sb_createTable(int rows, int cols);


/*!\brief Create an array of strings
 \param int length
 \return sb_strings
 \ingroup Basic
*/
SBAPIEXPORT sb_strings sb_createStringsArray(int len);


/*!\brief Create an array of items
 \param int number of items
 \return sb_items
 \ingroup Basic
*/
SBAPIEXPORT sb_items sb_createItemsArray(int len);

/*!\brief Get i,jth value from a sb_matrix
 \param sb_matriix matrix
 \param int row
 \param int column
 \return double value at the given row, column
 
 Use sb_setMatrixValue(M,i,j,value) to set a value

 \ingroup Basic
*/
SBAPIEXPORT double sb_getMatrixValue(sb_matrix M, int i, int j);


/*!\brief Set i,jth value of a sb_matrix
 \param sb_matrix matrix
 \param int row
 \param int column
 \param double value at the given row, column
 \ingroup Basic
*/
SBAPIEXPORT void sb_setMatrixValue(sb_matrix M, int i, int j, double d);


/*!\brief Get ith row name from a sb_matrix
 \param sb_matrix matrix
 \param int row
 \return string row name
 \ingroup Basic
*/
SBAPIEXPORT const char * sb_getRowName(sb_matrix M, int i);


/*!\brief Set ith row name for a sb_matrix
 \param sb_matrix matrix
 \param int row
 \param string row name
 \ingroup Basic
*/
SBAPIEXPORT void sb_setRowName(sb_matrix M, int i, const char * s);


/*!\brief Get jth column name of a sb_matrix
 \param sb_matrix matrix
 \param int column
 \return string column name
 \ingroup Basic
*/
SBAPIEXPORT const char * sb_getColumnName(sb_matrix M, int j);


/*!\brief Set jth column name of a sb_matrix
 \param sb_matrix matrix
 \param int column
 \param string column name
 \ingroup Basic
*/
SBAPIEXPORT void sb_setColumnName(sb_matrix M, int j, const char * s);


/*!\brief Get i,j-th string in a table
 \param sb_table table
 \param int row
 \param int column
 \return string value at row,column
 \ingroup Basic
*/
SBAPIEXPORT const char* sb_getTableValue(sb_table S, int i, int j);


/*!\brief Set i,jth string in a table
 \param sb_table table
 \param int row
 \param int column
 \param string value at row,column
 \ingroup Basic
*/
SBAPIEXPORT void sb_setTableValue(sb_table S, int i, int j, const char * s);


/*!\brief Get ith string in array of strings
 \param sb_strings array
 \param int index
 \return string value
 \ingroup Basic
*/
SBAPIEXPORT const char* sb_getString(sb_strings S, int i);


/*!\brief Set ith string in array of strings
 \param sb_strings array
 \param int index
 \param string value
 \ingroup Basic
*/
SBAPIEXPORT void sb_setString(sb_strings S, int i, const char * c);


/*!\brief Get ith long item in array of items
\param sb_items array
 \param int index
 \return long value
 \ingroup Basic
*/
SBAPIEXPORT long sb_getItem(sb_items A, int i);


/*!\brief Set ith long item in array of items
 \param sb_items array
 \param int index
 \param long value
 \ingroup Basic
*/
SBAPIEXPORT void sb_setItem(sb_items A, int i, long o);


/*!\brief Get the index of a string in the array
 \param sb_strings array
 \param char* a string in the array
 \return int index of that string
 \ingroup Basic
*/
SBAPIEXPORT int sb_getStringIndex(sb_strings A, const char * s);


/*!\brief get the row number of a row name
 \param sb_matrix matrix
 \param char* a string in the matrix
 \return int index of that string
 \ingroup Basic
*/
SBAPIEXPORT int sb_getRowIndex(sb_matrix, const char * s);


/*!\brief Get the column number of a column name
 \param sb_matrix matrix
 \param char* a string in the matrix
 \return int index of that string
 \ingroup Basic
*/
SBAPIEXPORT int sb_getColumnIndex(sb_matrix, const char * s);


/*!\brief Delete a matrix
 \param &sb_matrix pointer to matrix
 \ingroup Basic
*/
SBAPIEXPORT void sb_deleteMatrix (sb_matrix M);


/*!\brief Delete a strings table
 \param &sb_table pointer to table
 \ingroup Basic
*/
SBAPIEXPORT void sb_deleteTable(sb_table M);


/*!\brief Delete an array of items
 \param &sb_items pointer to array
 \ingroup Basic
*/
SBAPIEXPORT void sb_deleteItemsArray(sb_items A);


/*!\brief Delete an array of strings
 \param &sb_strings pointer to array
 \ingroup Basic
*/
SBAPIEXPORT void sb_deleteStringsArray(sb_strings C);


/*!\brief Combine two matrices by appending their columns. Row size must be equal for both matrices
 \param sb_matrix first matrix
 \param sb_matrix fsecond matrix
 \return sb_matrix new combined matrix
 \ingroup Basic
*/
SBAPIEXPORT sb_matrix sb_appendColumns(sb_matrix A, sb_matrix B);


/*!\brief Combine two matrices by appending their rows. Column sizes must be equal for both matrices
 \param sb_matrix first matrix
 \param sb_matrix fsecond matrix
 \return sb_matrix new combined matrix
 \ingroup Basic
*/
SBAPIEXPORT sb_matrix sb_appendRows(sb_matrix A, sb_matrix B);


/*!\brief Print a matrix to file
 \param char* file name
 \param sb_matrix
 \ingroup Basic
*/
SBAPIEXPORT void sb_printMatrixToFile(const char* file, sb_matrix M);


/*!\brief Print a matrix to stdout
 \param char* file name
 \param sb_matrix
 \ingroup Basic
*/
SBAPIEXPORT void sb_printOutMatrix(sb_matrix M);


/*!\brief Print a table to file
 \param char* file name
 \param sb_table
 \ingroup Basic
*/
SBAPIEXPORT void sb_printTableToFile(const char* file, sb_table M);


/*!\brief Print a table to stdout
 \param sb_table
 \ingroup Basic
*/
SBAPIEXPORT void sb_printOutTable(sb_table M);

END_C_DECLS
#endif

