#ifndef __SURFMAT_HPP__
#define __SURFMAT_HPP__

// TODO: template on ordinal type and scalar type

//#include "/home/briadam/dakota/trunk_cg/packages/teuchos/src/Teuchos_SerialDenseMatrix.hpp"
#include "/home/kdalbey/MatrixWrapperTest/teuchos_install/include/Teuchos_SerialDenseMatrix.hpp"

namespace nkm {

/** SurfMat_Teuchos class: an implementation of SurfMat that wraps a
 Teuchos SerialDenseMatrix<int, T>
*/
template<typename T> 
class SurfMat //SMTeuchos: public SurfMat<T>
{

public:
  
  // ---------
  // structors
  // ---------

  /// default constructor; does not allocate memory
  SurfMat();
  
  /// typical matrix constructor
  SurfMat(int Nrows_in, int Ncols_in = 1);
  
  /// copy constructor, make a deep copy
  SurfMat(const SurfMat<T>& other);
  
  /// destructor 
  ~SurfMat();


  // ----------------
  // sizing operators
  // ----------------

  /** Change this matrix to a matrix with nrows_new rows and ncols_new
      columns. It does not attempt to preserve values in memory and is
      a no-op if the size hasn't changed. */
  void newSize(int nrows_new, int ncols_new = 1);

  /** Enlarge or shrink the matrix while keeping contigous in memory
      elements contiguous in memory; no-op is no change in size.  If
      matrix grows, optionally zero entries */
  void reshape(int nrows_new, int ncols_new = 1);
  
  /** Enlarge or shrink the matrix while adding zeros after the last
      row and/or last column and/or chopping off the rows and/or
      columns after the newly requested last row and/or column.  No-op
      if size unchanged */
  void resize(int nrows_new, int ncols_new = 1);
    

  // ------------
  // initializers
  // ------------

  /// set all of the matrix's elements to zero
  void zero();

  /// Set all of the matrix's elements to zero, except for the
  /// diagonal (i=j) which is set to 1 (including for rectangular matrices)
  void identity();


  // -------
  // copiers
  // -------
  
  /// Make a deep copy, returning self
  SurfMat<T>& copy(const SurfMat<T>& other);

  /// Assignment operator.  Performs deep copy.
  SurfMat<T>& operator=(const SurfMat<T>& other);
  
  // ------------
  // data set/get
  // ------------

  // TODO: add convenience functions to avoid returning elements by
  // pointer/non-const ref

  /// return the number of rows
  int getNRows() const;
  
  /// return the number of columns
  int getNCols() const;
  
  /// return the number of elements
  int getNElems() const;

  /// get the comparison tolerance
  const T getTol() const;

  /// set the comparison tolerance
  void putTol(T tol_in);
    
  /// vector style retrieve of a single element from the matrix by value
  const T& operator()(int k) const;
  
  /// vector style retrieve of a single element from the matrix by reference
  T& operator()(int k);

  /// matrix style retrieve of a single element from the matrix by value
  const T& operator()(int i, int j) const;
  
  /// matrix style retrieve of a single element from the matrix by reference
  T& operator()(int i, int j);

  /// vector style retrieve of pointer to element, for passing to BLAS &
  /// LAPACK convenience
  T* ptr(int k);

  /// vector style retrieve of pointer to element, for passing to BLAS &
  /// LAPACK convenience
  const T* ptr(int k) const;

  /// matrix style retrieve of pointer to element, for passing to BLAS &
  /// LAPACK convenience
  T* ptr(int i, int j);
  
  /// matrix style retrieve of pointer to element, for passing to BLAS &
  /// LAPACK convenience
  const T* ptr(int i, int j) const;


  // ----------------
  // member functions
  // ----------------

  /// returns the smallest value contained in this matrix
  T minElem() const;
  
  /// returns the location of the smallest value contained in this
  /// matrix, if there is a tie it returns the location of the first
  /// one encountered
  int iMinElem() const;
  
  /// returns the smallest value contained in this matrix and its
  /// location in loc, if there is a tie it records the location of the
  /// first smallest value encountered
  void minElem(T& val, int& loc) const;
  
  /// returns the largest value contained in this matrix
  T maxElem() const;
  
  /// returns the location of the largest value contained in this
  /// matrix, if there is a tie it returns the location of the first
  /// one encountered
  int iMaxElem() const;

  /// returns the largest value contained in this matrix and its
  /// location in loc, if there is a tie it records the location of the
  /// first largest value encountered
  void maxElem(T& val, int& loc) const;
  
  /// returns the smallest and largest values contained in this matrix
  void minMaxElem(T& smallest, T& largest) const;


  // ---------------
  // row/col get/set
  // ---------------

  
  /// replace the values in row "irow" of this matrix with the values
  /// stored in "row"; this function will NOT expand the size of this
  /// matrix if you specify a row index larger than NRows-1.
  void putRows(SurfMat<T>& row, int irow);
 
  /// replace the values in multiple rows (which rows they are is
  /// specified in irows) of this matrix with the values stored in
  /// "rows"; this function will NOT expand the size of this matrix if
  /// you specify a row index larger than NRows-1.
  void putRows(SurfMat<T>& rows, SurfMat<int> irows);

  /// replace the values in column "jcol" of this matrix with the
  /// values stored in "col"; this function will NOT expand the size of
  /// this matrix if you specify a column index larger than NCols-1.
  void putCols(SurfMat<T>& col, int jcol);

  ///replace the values in multiple columns (which columns they are is
  ///specified in jcols) of this matrix with the values stored in
  ///"cols"; this function will NOT expand the size of this matrix if
  ///you specify a column index larger than NCols-1.
  void putCols(SurfMat<T>& cols, SurfMat<int> jcols);

  /// get one row (index stored in irow) of this matrix and return as a
  /// "row vector"
  SurfMat<T>& getRows(SurfMat<T>& result, int irow) const;

  /// get multiple rows (indices stored in matrix irows) of this matrix
  /// and return as a new matrix
  SurfMat<T>& getRows(SurfMat<T>& result, SurfMat<int>& irows) const;

  /// get one column (index stored in jcol) of this matrix and return
  /// as a vector
  SurfMat<T>& getCols(SurfMat<T>& result, int jcol) const;

  /// get multiple columns (indices stored in matrix jcols) of this
  /// matrix and return as a new matrix
  SurfMat<T>& getCols(SurfMat<T>& result, SurfMat<int>& jcols) const;
  
  /// returns a copy of the matrix excluding 1 row whose index is
  /// stored in irow, this isn't an inline because you will need to
  /// copy all but 1 row
  SurfMat<T>& excludeRows(SurfMat<T>& result, int irow) const;
  
  /// returns a copy of the matrix excluding 1 row whose index is
  /// stored in irow, this isn't an inline because for the typical use
  /// you will need to copy most of the matrix
  SurfMat<T>& excludeRows(SurfMat<T>& result, SurfMat<int>& irows) const;
  
  /// returns a copy of the matrix excluding 1 column whose index is
  /// stored in jcol, this isn't an inline because you will need to
  /// copy all but 1 column
  SurfMat<T>& excludeCols(SurfMat<T>& result, int jcol) const;
  
  /// return a copy of the matrix that excludes all columns whose
  /// indices are stored in the matrix irows, this isn't an inline
  /// because for the typical use case you will need to copy most of
  /// the matrix
  SurfMat<T>& excludeCols(SurfMat<T>& result, SurfMat<int>& jcols) const;
  
  
  // ------------
  // compare/sort
  // ------------

  // TODO: decide which of these should be public (probably few)

  /// used in the generic matrix quicksort; returns -1 a-b<-tol, 0 if
  /// -tol<=a-b<=tol, +1 if tol<a-b
  int compareElemAElemB(int ia, int ib) const;
  
  /// used in the generic matrix quicksort; swap the ia-th and ib-th
  /// elements of the "vector"
  void swapElems(int ia, int ib);
  
  /// ascending sort the elements of the matrix as if it were a vector,
  /// could easily add an integer argument (+1 or -1) that indicates
  /// whether you want to sort into ascending or descending order,
  /// would need to add a field to the matrix though
  void sortElems();  

  void qsortElems(int istart, int istop);
  void qsortRows(int istart, int istop);
  void qsortCols(int istart, int istop);

  /// performs an ascending unique sort of the elements, eliminates
  /// duplicates, and reshapes it into a (shruken, if appropriate)
  /// vector
  void uniqueElems();
  
  /// used in the generic matrix quicksort; compares 2 rows (a and b)
  /// starting with column 0; it returns -1 if a-b<-tol, or +1 if
  /// tol<a-b, if -tol<=a-b<=tol it moves onto column to and so on, iff
  /// every column of the 2 rows are equal then it returns 0, could add
  /// a column order field to the matrix to add the capability to
  /// change the order in which the rows are sorted and whether or not
  /// you wanted to do an ascending or descending sort
  int compareRowARowB(int irowa, int irowb) const;

  /// used in the generic matrix quicksort; swap rows a and b of the matrix
  void swapRows(int irowa, int irowb);
  
  /// ascending sort the rows of the matrix (comparing from left to
  /// right), you could easily add a MtxInt argument with 2 rows the
  /// first indicating order of columns being compared, the second
  /// being +1 or -1 to indicate to sort that row into ascending or
  /// descending order, would need to add a field to the matrix though
  void sortRows();
  
  /// performs an ascending sort of the rows, eliminates duplicates,
  /// and shrinks the matrix if it is appropriate to do so
  void uniqueRows();
  
  /// used in the generic matrix quicksort; compares 2 columns starting
  /// with row 0; it returns -1 if a-b<-tol, or +1 tol<a-b, if
  /// -tol<=a-b<=tol it moves onto the next row, iff every row of the 2
  /// columns have -tol<=a-b<=tol then it returns 0, could add a row
  /// order field to the matrix to add the capability to change the
  /// order in which the columns are sorted and whether or not you
  /// wanted to do an ascending or descending sort
  int compareColAColB(int jcola, int jcolb) const;

  /// used in the generic matrix quicksort; swap columns a and b of the matrix
  void swapCols(int jcola, int jcolb);
  
  /// ascending sort the coluns of the matrix (comparing from row zero
  /// onward), you could easily add a MtxInt argument with 2 columns
  /// the first indicating order of columns being compared (some rows
  /// can be omitted if you don't want to include them as sorting
  /// criteria), the second column being +1 or -1 to indicate to sort
  /// that row into ascending or descending order, would need to add a
  /// field to the matrix though
  void sortCols();
  
  /// performs an ascending sort of the columns, eliminates duplicates,
  /// and shrinks the matrix if it is appropriate to do so
  void uniqueCols();


private:

  /// resize to zero
  void clear();

  Teuchos::SerialDenseMatrix<int, T> tsdm;

  T tol; // an inequaltiy tolerance for equality checking should be 0
	 // for integers, WARNING TO WHOEVER FINISHES IMPLEMENTING
	 // THIS, you will have to decide what to do with tol when
	 // you clear(), newsize(), reshape(), or resize() the
	 // matrix, you will likely need to modify the contructors
  
  // this shouldn't be necessary but do it for good measure
  friend void mtxqsort(int istart, int istop,
		       int (* compare_a_b)(int ia, int ib),
		       void (* swap_a_b)(int ia, int ib));

};


// Definitions
#include <algorithm>


// ---------
// structors
// ---------

template< typename T >
SurfMat<T>::SurfMat(): 
  tol(0) //of type T, so one type conversion every construction
{ /* empty constructor */ }


template< typename T >
SurfMat<T>::SurfMat(int nrows_in, int ncols_in)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_in>=0)&&(ncols_in>=0));
#endif

#ifdef __SURFMAT_ZERO_MEM__
  tsdm.shape(nrows_in, ncols_in);
#else
  tsdm.shapeUninitialized(nrows_in, ncols_in);
#endif
}


template< typename T >
SurfMat<T>::SurfMat(const SurfMat<T>& other)
{
  copy(other);
}


template< typename T >
SurfMat<T>::~SurfMat() 
{
  /* empty destructor */  
}


// ----------------
// sizing operators
// ----------------

template< typename T >
void SurfMat<T>::clear() 
{ 
  newSize(0,0); 
}


template< typename T >
void SurfMat<T>::newSize(int nrows_new, int ncols_new) 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif

  if( (getNRows() != nrows_new) || (getNCols() != ncols_new) ) {
#ifdef __SURFMAT_ZERO_MEM__
    tsdm.shape(nrows_new, ncols_new);
#else
    tsdm.shapeUninitialized(nrows_new, ncols_new);
#endif
  }
  // TODO: if same total data, avoid realloc?  Probably not since
  // I'd need a new TSDM.
}


// TODO: suspect code here (need to preserve continguous memory)
template< typename T >
void SurfMat<T>::reshape(int nrows_new, int ncols_new) 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif

  // TODO: not sure how best to do this with TSDM
  if( (getNRows()!=nrows_new) || (getNCols()!=ncols_new) ) {
    int nelem_new = nrows_new*ncols_new;
    if (getNElems() >= nelem_new) {
      tsdm = Teuchos::SerialDenseMatrix<int, T>(Teuchos::Copy, 
						tsdm.values(), 
						nrows_new, //stride 
						nrows_new, ncols_new);
    }
    else {
      // matrix is growing; want trailing zeros
      // tmp will be zeroed out
      Teuchos::SerialDenseMatrix<int, T> tmp(nrows_new, ncols_new, true);
      T* ldata = tmp.values();
      for (int k=0; k<nelem_new; ++k)
	ldata[k] = operator()(k);
      tsdm.shapeUninitialized(nrows_new, ncols_new);
      //      tsdm = tmp;
      tsdm.assign(tmp);
    }
  }

}


template< typename T >
void SurfMat<T>::resize(int nrows_new, int ncols_new) 
{
  if( (getNRows() != nrows_new) || (getNCols() != ncols_new) ) 
    tsdm.reshape(nrows_new, ncols_new);
}


// ------------
// initializers
// ------------

template< typename T >
void SurfMat<T>::zero() 
{ 
  tsdm = 0; 
}


template< typename T >
void SurfMat<T>::identity() 
{
  zero(); 
  int n = (getNRows() < getNCols()) ? getNRows() : getNCols(); 
  for (int i=0; i<n; ++i)
    tsdm(i, i) = 1;
}


// -------
// copiers
// -------
  
template< typename T >
SurfMat<T>& SurfMat<T>::copy(const SurfMat<T>& other)
{
  
  // need to size (can't do tsdm = other.tsdm;)
  tsdm.shapeUninitialized(other.tsdm.numRows(),
			  other.tsdm.numCols());
  tsdm = other.tsdm;

  tol = other.getTol();
  
  return *this;
}


template< typename T >
SurfMat<T>& SurfMat<T>::operator=(const SurfMat<T>& other) 
{
  return copy(other);
}


// ------------
// data set/get
// ------------

template< typename T >
int SurfMat<T>::getNRows() const 
{ 
  return tsdm.numRows(); 
}
  

template< typename T >
int SurfMat<T>::getNCols() const 
{ 
  return tsdm.numCols(); 
}
  

template< typename T >
int SurfMat<T>::getNElems() const 
{ return getNRows()*getNCols(); 
}


template< typename T >
const T SurfMat<T>::getTol() const 
{ 
  return tol; 
}


template< typename T >
void SurfMat<T>::putTol(T tol_in) 
{
  tol = tol_in; 
}


template< typename T >
const T& SurfMat<T>::operator()(int k) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=k)&&(k<(getNRows()*getNCols())));
#endif
  return tsdm.values()[k];
}
  

template< typename T >
T& SurfMat<T>::operator()(int k)
{
#ifdef __SURFMAT_ERR_CHECK__
  if(!((0<=k)&&(k<getNRows()*getNCols())))
    assert((0<=k)&&(k<getNRows()*getNCols()));
#endif
  return tsdm.values()[k];
}
  

template< typename T >
const T& SurfMat<T>::operator()(int i, int j) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=i)&&(i<getNRows())&&(0<=j)&&(j<getNCols()));
#endif
  return tsdm[j][i];
}
  

template< typename T >
T& SurfMat<T>::operator()(int i, int j) 
{
#ifdef __SURFMAT_ERR_CHECK__
  if( !( (0<=i) && (i<getNRows()) &&
	 (0<=j) && (j<getNCols())    ) ){
    printf("i=%d NRows=%d j=%d NCols=%d\n",i,getNRows(),j,getNCols());
    assert((0<=i)&&(i<getNRows())&&(0<=j)&&(j<getNCols()));
  }
#endif
  return  tsdm[j][i];
}
  

template< typename T >
T* SurfMat<T>::ptr(int k) 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=k)&&(k<getNRows()*getNCols()));
#endif
  return tsdm.values()+k;
}


template< typename T >
const T* SurfMat<T>::ptr(int k) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=k)&&(k<getNRows()*getNCols()));
#endif
  return tsdm.values()+k;
};


template< typename T >
T* SurfMat<T>::ptr(int i, int j) 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=i)&&(i<getNRows())&&(0<=j)&&(j<getNCols()));
#endif
  return tsdm[j]+i;
}


template< typename T >
const T* SurfMat<T>::ptr(int i, int j) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=i)&&(i<getNRows())&&(0<=j)&&(j<getNCols()));
#endif
  return tsdm[j]+i;
}


// ----------------
// member functions
// ----------------

template< typename T >
inline T SurfMat<T>::minElem() const 
{
  T val;
  int loc;
  minElem(val,loc);
  return val;
}


template< typename T >
inline int SurfMat<T>::iMinElem() const 
{
  T val;
  int loc;
  minElem(val,loc);
  return loc;
}


template< typename T >
inline void SurfMat<T>::minElem(T& val, int& loc) const 
{
  int nelem = getNElems();
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem>0);
#endif
  val = operator()(0);
  loc = 0;
  for(int k=1; k<nelem; ++k)
    if(operator()(k) < val){
      val = operator()(k);
      loc = k;
    }
}


template< typename T >
inline T SurfMat<T>::maxElem() const 
{
  T val;
  int loc;
  maxElem(val,loc);
  return val;
}


template< typename T >
inline int SurfMat<T>::iMaxElem() const
{
  T val;
  int loc;
  maxElem(val,loc);
  return loc;
}


template< typename T >
inline void SurfMat<T>::maxElem(T& val, int& loc) const 
{
  int nelem = getNElems();
#ifdef __SURFMAT_ERR_CHECK__
  assert(nelem>0);
#endif
  val = operator()(0);
  loc = 0;
  for(int k=1; k<nelem; ++k)
    if(operator()(k) > val){
      val = operator()(k);
      loc = k;
    }
}


template< typename T >
void SurfMat<T>::minMaxElem(T& smallest, T& largest) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert(getNElems()>0);
#endif
  T* ldata = tsdm.values();
  smallest = largest = ldata[0];
  for(int k=1; k<getNElems(); ++k) {
    // note I did not use an else if for the largest so that both if
    // statements could be evaluated simultaneous by different
    // processing units
    if(ldata[k] < smallest)
      smallest = ldata[k];
    if(ldata[k] > largest)
      largest = ldata[k];
  }    
}
  

// ---------------
// row/col get/set
// ---------------

template< typename T >
void SurfMat<T>::putRows(SurfMat<T>& row, int irow)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=irow)&&(irow<getNRows())&&
	 (row.getNRows()==1)&&(row.getNCols()==getNCols()));
#endif
  for(int j=0; j<getNCols(); ++j)
    tsdm(irow, j) = row(j);
}
  

template< typename T >
void SurfMat<T>::putRows(SurfMat<T>& rows, SurfMat<int> irows)
{
  int nrows_put = irows.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=irows.minElem())&&(irows.maxElem()<getNRows())&&
	 (rows.getNRows()==nrows_put)&&(rows.getNCols()==getNCols()));
#endif
  for(int j=0; j<getNCols(); ++j)
    for(int k=0; k<nrows_put; ++k)
      tsdm(irows(k), j) = rows(k,j);
}
 
 
template< typename T >
void SurfMat<T>::putCols(SurfMat<T>& col, int jcol)
{
#ifdef __SURFMAT_ERR_CHECK__
  if(!((0<=jcol)&&(jcol<getNCols())&&
       (col.getNRows()==getNRows())&&(col.getNCols()==1))){
    //assert((0<=jcol)&&(jcol<NCols)&&
    //(col.getNRows()==NRows)&&(col.getNCols()==1));
    assert(0<=jcol);
    assert(jcol<getNCols());
    assert(col.getNRows()==getNRows());
    assert(col.getNCols()==1);
  }
#endif
  for(int i=0; i<getNRows(); ++i)
    tsdm[jcol][i] = col(i);
}


template< typename T >
void SurfMat<T>::putCols(SurfMat<T>& cols, SurfMat<int> jcols) 
{
  int ncols_put=jcols.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=jcols.minElem())&&(jcols.maxElem()<getNCols())&&
	 (cols.getNRows()==getNRows())&&(cols.getNCols()==ncols_put));
#endif
  for(int k=0; k<ncols_put; ++k)
    for(int i=0; i<getNRows(); ++i)
      tsdm(i, jcols(k)) = cols(i, k);
}
  

template< typename T >
SurfMat<T>& SurfMat<T>::getRows(SurfMat<T>& result, int irow) const
{
#ifdef __SURFMAT_ERR_CHECK__
  if(!((0<=irow)&&(irow<getNRows()))) {
    printf("irow=%d NRows=%d\n",irow,getNRows()); fflush(stdout);
    assert((0<=irow)&&(irow<getNRows()));
  }
#endif
  result.newSize(1,getNCols());
  result.putTol(tol);
  for(int j=0; j<getNCols(); ++j)
    result(j) = tsdm[j][irow];
  return result;
}
  

template< typename T >
SurfMat<T>& SurfMat<T>::getRows(SurfMat<T>& result, SurfMat<int>& irows) const
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=irows.minElem())&&(irows.maxElem()<getNRows()));
#endif
  int nrows_res = irows.getNElems();
  result.newSize(nrows_res, getNCols());
  result.putTol(tol);
  if(nrows_res > 0) 
    for(int j=0; j<getNCols(); ++j)
      for(int i=0; i<nrows_res; ++i)
	result(i,j) = tsdm[j][irows(i)];	
  return result;
}
  

template< typename T >
SurfMat<T>& SurfMat<T>::getCols(SurfMat<T>& result, int jcol) const 
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=jcol)&&(jcol<getNCols()));
#endif
  result.newSize(getNRows());
  result.putTol(tol);
  for (int i=0; i<getNRows(); ++i) 
    result(i) = tsdm[jcol][i];
  return result;
}
 
 
template< typename T >
SurfMat<T>& SurfMat<T>::getCols(SurfMat<T>& result, SurfMat<int>& jcols) const
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=jcols.minElem())&&(jcols.maxElem()<getNCols()));
#endif
  int ncols_res = jcols.getNElems();
  result.newSize(getNRows(), ncols_res); 
  result.putTol(tol);
  if(ncols_res>0)
    for(int j=0; j<ncols_res; ++j)
      for(int i=0; i<getNRows(); ++i)
	result(i,j) = tsdm[jcols(j)][i];	
    
  return result;
}


template< typename T >
SurfMat<T>& SurfMat<T>::excludeRows(SurfMat<T>& result, int irow) const
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=irow)&&(irow<getNRows()));
#endif
  if(getNRows()==1)
    result.clear();
  else {
    result.newSize(getNRows()-1, getNCols());
    result.putTol(tol);
    int isrc, ikeep, j;
    for(j=0; j<getNCols(); ++j) {      
      for(isrc=ikeep=0; isrc<irow; ++isrc, ++ikeep)
	result(ikeep, j) = tsdm[j][isrc];
      isrc=irow+1;
      for(;isrc<getNRows();++isrc,++ikeep)
	result(ikeep, j) = tsdm[j][isrc];
    }
  }
  return result;
}


template< typename T > SurfMat<T>& SurfMat<T>::
excludeRows(SurfMat<T>& result, SurfMat<int>& irows) const
{
  int j;
  int nexclude=irows.getNElems();
  if(nexclude<1) {
    // the list of rows to exclude is empty so copy over the whole matrix
    // TODO: use copy() method
    result.newSize(getNRows(),getNCols());
    result.putTol(tol);
    for(j=0; j<getNCols(); ++j)
      for(int i=0; i<getNRows(); ++i)
	    result(i,j) = tsdm[j][i];
  }
  else{
    irows.uniqueElems(); //sort the rows to exclude into ascending
			 //order and eliminate duplicate listings
    nexclude=irows.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irows(0))&&(irows(nexclude-1)<getNRows()));
#endif
    if(nexclude==getNRows()) {
      //the user wants us to eliminate _all_ rows
      result.clear();
    }
    else {
      //the user wants us to eliminate some but not all rows
      result.newSize(getNRows()-nexclude, getNCols());
      result.putTol(tol);
      int iexclude, ikeep, isrc;
      for(j=0;j<getNCols();++j) {
	iexclude=ikeep=isrc=0;
	while(isrc < getNRows()) {
	  if(iexclude<nexclude) {
	    for( ; isrc<irows(iexclude) ; ++isrc, ++ikeep)
	      result(ikeep, j) = tsdm[j][isrc];
	    //at this point isrc=irows(iexclude)
	    ++iexclude;
	    ++isrc;
	  }
	  else{
	    for(;isrc<getNRows();++isrc, ++ikeep)
	      result(ikeep,j) = tsdm[j][isrc];
	    //at this point isrc=NRows and the while loop will terminate
	  }
	}//while loop terminates
      }//do the same thing with the next column  
    }
  }
  return result;
}


template< typename T >
SurfMat<T>& SurfMat<T>::excludeCols(SurfMat<T>& result, int jcol) const
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((0<=jcol)&&(jcol<getNCols()));
#endif
  if(getNCols()==1)
    result.clear();
  else{
    result.newSize(getNRows(), getNCols()-1);
    result.putTol(tol);
    int jsrc, jkeep, i;
    for(jsrc=jkeep=0; jsrc<jcol; ++jsrc, ++jkeep)
      for(i=0; i<getNRows(); ++i)
	result(i, jkeep) = tsdm[jsrc][i];
    ++jsrc;
    for(; jsrc<getNRows(); ++jsrc, ++jkeep)
      for(i=0; i<getNRows(); ++i)
	result(i, jkeep) = tsdm[jsrc][i];
  }
  return result;
}


template< typename T > SurfMat<T>& SurfMat<T>::
excludeCols(SurfMat<T>& result, SurfMat<int>& jcols) const
{
  int nexclude=jcols.getNElems();
  int i;
  if(nexclude<1) {
    result.newSize(getNRows(), getNCols());
    result.putTol(tol);
    for(int j=0;j<getNCols();++j)
      for(i=0;i<getNRows();++i)
	result(i,j) = tsdm[j][i];
  }
  else{
    jcols.uniqueElems();
    nexclude=jcols.getNElems();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcols(0))&&(jcols(nexclude-1)<getNCols()));
#endif
    if(nexclude==getNCols()) {
      //the user wants us to eliminate _all_ columns
      result.clear(); 
    }
    else {
      //the user wants us to eliminate some but not all columns
      result.newSize(getNRows(),getNCols()-nexclude);
      result.putTol(tol);
      int jexclude, jkeep, jsrc;
      jexclude=jkeep=jsrc=0;
      while(jsrc<getNCols()) {
	if(jexclude<nexclude) {
	  for(;jsrc<jcols(jexclude);++jsrc,++jkeep)
	    for(i=0;i<getNRows();++i)
	      result(i,jkeep) = tsdm[jsrc][i];
	  //at this point jsrc=jrows(jexclude)
	  ++jexclude;
	  ++jsrc;
	}
	else{
	  for(;jsrc<getNCols();++jsrc,++jkeep)
	    for(i=0;i<getNRows();++i)
	      result(i,jkeep) = tsdm[jsrc][i];
	  //at this point jsrc=NCols and the while loop will terminate
	}
      }//while loop terminates
    }
  }
  return result;
}


// ------------
// compare/sort
// ------------

template< typename T >
int SurfMat<T>::compareElemAElemB(int ia, int ib) const
{
  const T diff = operator()(ia) - operator()(ib);
  return(-(diff<-tol)+(tol<diff));
}

  
template< typename T >
void SurfMat<T>::swapElems(int ia, int ib)
{
  std::swap(operator()(ia), operator()(ib));
}

  
template< typename T >
void SurfMat<T>::sortElems()
{
  qsortElems(0, getNElems()-1);
}


template< typename T >
void SurfMat<T>::qsortElems(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapElems(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareElemAElemB(i,istart)<=0))
	++i;
      while((j > istart) && (compareElemAElemB(istart,j)<0))
	--j;
      if( i < j)
	swapElems(i,j);
    }
    // swap two elements
    swapElems(istart,j);
    // recursively sort the lesser list
    qsortElems(istart,j-1);
    qsortElems(j+1,istop);
  }
  return;
}
  

template< typename T >
void SurfMat<T>::qsortRows(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapRows(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareRowARowB(i,istart)<=0))
	++i;
      while((j > istart) && (compareRowARowB(istart,j)<0))
	--j;
      if( i < j)
	swapRows(i,j);
    }
    // swap two elements
    swapRows(istart,j);
    // recursively sort the lesser list
    qsortRows(istart,j-1);
    qsortRows(j+1,istop);
  }
  return;
}


template< typename T >
void SurfMat<T>::qsortCols(int istart, int istop)
{
  int i,j,k;
  if( istart < istop) {
    k = (istart+istop)/2;
    swapCols(istart,k);
    i = istart+1;
    j = istop;
    while(i <= j){
      while((i <= istop) && (compareColAColB(i,istart)<=0))
	++i;
      while((j > istart) && (compareColAColB(istart,j)<0))
	--j;
      if( i < j)
	swapCols(i,j);
    }
    // swap two elements
    swapCols(istart,j);
    // recursively sort the lesser list
    qsortCols(istart,j-1);
    qsortCols(j+1,istop);
  }
  return;
}


template< typename T >
void SurfMat<T>::uniqueElems()
{
  int nelems=getNElems();
  if(nelems>0) {
    sortElems();    
    int i,k;
    i=1; 
    while(i<nelems) {
      if(operator()(i) == operator()(i-1)) {
	for(k=i+1; k<nelems; ++k)
	  operator()(k-1)=operator()(k);
	--nelems;}
      else ++i;
    }
    reshape(nelems);
  }
}


template< typename T >
int SurfMat<T>::compareRowARowB(int irowa, int irowb) const
{
  T diff;
  int j=0;
  do {
    diff = tsdm[j][irowa] - tsdm[j][irowb];
    ++j;
  } while((-tol<=diff)&&(diff<=tol)&&(j<getNCols()));
  return (-(diff<-tol)+(tol<diff));  //could multiply outermost ()
  //by +1/-1 for
  //ascending/descending
}
  

template< typename T >
void SurfMat<T>::swapRows(int irowa, int irowb)
{
  T swap;
  for(int j=0; j<getNCols(); ++j) {
    std::swap(tsdm[j][irowa], tsdm[j][irowb]);
  }
}
 

template< typename T >
void SurfMat<T>::sortRows()
{
  qsortRows(0, getNRows()-1);
}

template< typename T >
void SurfMat<T>::uniqueRows()
{
  int nrows = getNRows();
  if(nrows>0) {
    sortRows();
    int i,j,k;
    i=1;
    while(i<nrows) {
      if(compareRowARowB(i,i-1)) {
	for(j=0;j<getNCols();++j) 
	  for(k=i+1;k<nrows;++k)
	    tsdm[j][k-1] = tsdm[j][k];
	--nrows;
      }
      else 
	++i;
    }
    resize(nrows, getNCols()); //not an error, want non-contiguous
    //memory chopping
  }
}
  

template< typename T >
int SurfMat<T>::compareColAColB(int jcola, int jcolb) const
{
  T diff;
  int i=0;
  do {
    diff=tsdm[jcola][i] - tsdm[jcolb][i];
    ++i;
  } while ((-tol<=diff)&&(diff<=tol)&&(i<getNRows()));
  return (-(diff<-tol)+(tol<diff));
}


template< typename T >
void SurfMat<T>::swapCols(int jcola, int jcolb)
{
  T swap;
  for(int i=0; i<getNRows(); ++i)
    std::swap(tsdm[jcola][i], tsdm[jcolb][i]);
}

template< typename T >
void SurfMat<T>::sortCols()
{
  qsortCols(0, getNCols()-1);
}

template< typename T >
void SurfMat<T>::uniqueCols()
{
  int ncols = getNCols();
  if(ncols>0) {
    sortCols();
    int i,j,k;
    j=1;
    while(j<ncols) {
      if(compareColAColB(j,j-1)) {
	for(k=j+1;k<ncols; ++k) 
	  for(i=0;i<getNRows(); ++i)
	    tsdm[k-1][i] = tsdm[k][i];
	--ncols;
      }
      else ++j;
    }
    reshape(getNRows(),ncols); // not an error, want contiguous memory
    // chopping, reshape will call reshape2,
    // resize would call resize2 which would
    // call reshape2
  }
}

}

#endif
