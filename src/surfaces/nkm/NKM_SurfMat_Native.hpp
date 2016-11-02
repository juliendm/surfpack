#ifndef __SURFMAT_HPP__
#define __SURFMAT_HPP__
#include <cstdlib>
#include <vector>
//#include <iostream>
#include <string>
#include <sstream>
#include "surfpack_system_headers.h"

namespace nkm {


  //    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols));


/**
 * The SurfMat Class template
 */
template<typename T>
class SurfMat
{
public:
  /**
   * Default constructor, does not allocate memory
   */
  SurfMat() {
    tol  =static_cast<T>(0);  /**< of type T, so one type conversion every construction */
    NRowsAct=NColsAct=NRows=NCols=0;
    return;
  };
  
  /*
   * Typical matrix constructor, defaults to 1 column (a vector)
   * if the number of columns is not specified
   *
   * @param nrows_in
   *            the number of rows
   * @param ncols_in
   *            the number of columns
   */
  SurfMat(int nrows_in, int ncols_in);
  
  /**
   * Copy constructor, makes a deep copy
   *
   * @param other
   *            object to copy from
   */
  SurfMat(const SurfMat<T>& other);
  
  /**
   * Copy from vector constructor
   *
   * @param vecT
   *            Vector to copy from
   * @param nrows_in
   *            number of rows to copy from vector
   * @param ncols_in
   *            number of columns
   */
  SurfMat(const std::vector<T> &vecT, int nrows_in=-99999, int ncols_in=1);

  /**
   * Destructor
   */
  ~SurfMat() {
    clear();
    return;
  };
  
  /**
   * Deallocates memory
   */
  void clear();
  
  /**
   * Change this matrix to a matrix with nrows_new apparent rows and ncols_new apparent columns.
   *
   * It does not attempt to preserve values in memory which,
   * if a change is needed, makes it faster than reshape2 or resize2:
   * if the matrices actual number of rows and actual number of columns it will
   * just change the apparent number of rows and apparent number of columns without
   * actually changing the memory foot print unless the user forces memory reallocation
   * by specifying if_force=true.  Normally, this function won't do anything if
   * you request a new size the same as the current size. The user can force the
   * memory footprint of this matrix to shrink to exactly the requested amount by
   * specifying if_force=true
   *
   * @param nrows_new
   *        the number of rows to change to
   * @param ncols_new
   *        the number of columns to change to
   * @param if_force
   *        forces memory footprint to shrink if true
   */
  inline void newSize(int nrows_new, int ncols_new, bool if_force=false){
    if((NRows!=nrows_new)||(NCols!=ncols_new)||
       ((if_force==true)&&((nrows_new!=NRowsAct)||(ncols_new!=NColsAct)))) {      
      if((nrows_new<=NRowsAct)&&(ncols_new<=NColsAct)&&
	 (if_force==false)) {
	NRows=nrows_new;
	NCols=ncols_new;
      }
      else
	newSize2(nrows_new,ncols_new,if_force);
    }
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
    return;};

  /**
   * Enlarge or shrink the matrix while keeping contiguous in memory elements contiguous in memory
   *
   * This function doesn't do anything if you request a new size the same as the current size
   *
   * @param nrows_new
   *        the number of rows to change to
   * @param ncols_new
   *        the number of columns to change to
   * @param if_force
   *        forces memory footprint to shrink if true
   */
  inline void reshape(int nrows_new, int ncols_new, bool if_force=false){
    if((NRows!=nrows_new)||(NCols!=ncols_new)||
       ((if_force==true)&&((nrows_new!=NRowsAct)||(ncols_new!=NColsAct)))) {
      if((NRows==nrows_new)&&(ncols_new<=NColsAct)&&
	 (if_force==false))
	NCols=ncols_new;
      else
	reshape2(nrows_new,ncols_new,if_force);
    }
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  return;};

  /**
   * Enlarge or shrink the matrix while adding zeros
   *
   * Enlarge or shrink the matrix while adding zeros after the last row and/or last column
   * and/or chopping off the rows and/or columns after the newly requested last row
   * and/or column: this function doesn't do anything if you request a new size the same
   * as the current size.
   *
   * @param nrows_new
   *            the number of rows to change to
   * @param ncols_new
   *            the number of columns to change to
   * @param if_force
   *            forces memory footprint to shrink if true
   */
  inline void resize(int nrows_new, int ncols_new, bool if_force=false){
    if((NRows!=nrows_new)||(NCols!=ncols_new)||
       ((if_force==true)&&((nrows_new!=NRowsAct)||(ncols_new!=NColsAct)))) {
      if((nrows_new<=NRowsAct)&&(ncols_new<=NColsAct)&&
	 (if_force==false)) {
	NRows=nrows_new;
	NCols=ncols_new;
      }
      else
	resize2(nrows_new,ncols_new,if_force);
    }
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
    return;
  };
  
  /**
   * Set all of the matrix's elements to zero
   */
  //void zero(){int nelems=getNElems(); for(int k=0; k<nelems; k++) data[k]=0; return;};
  inline void zero(){int nelems=getNElemsAct(); for(int k=0; k<nelems; k++) data[k]=0; return;};

  /**
   * Set all of the matrix's elements to zero, except for the diagonal (i=j) which is set to 1
   *
   * @note Works for rectangular matrices
   */
  inline void identity(){
    zero(); 
    int n=(NRows<NCols)?NRows:NCols; 
    for(int i=0; i<n; i++) 
      data[jtoi[i]+i]=static_cast<T>(1);
    return;
  }
  
  /**
   * Make a deep copy
   *
   * @param other
   *        object to copy
   * @param if_force
   *        forces memory footprint to shrink if true
   *
   * @return new object that's a copy of the matrix
   */
  SurfMat<T>& copy(const SurfMat<T>& other, bool if_force=false);


  /**
   *  Assignment operator.  Performs deep copy.
   *
   *  @return new object that's a copy
   */
  inline SurfMat<T>& operator=(const SurfMat<T>& other){return (copy(other));};

  /**
   * Get the actual number of rows
   *
   * @return the actual number of rows
   */
  inline int getNRowsAct() const {return NRowsAct;};
  
  /**
   * Get the apparent number of rows
   *
   * @return the apparent number of rows
   */
  inline int getNRows() const {return NRows;};

  /**
   * Get the actual number of columns
   *
   * @return the actual number of columns
   */
  inline int getNColsAct() const {return NColsAct;};
  
  /**
   * Get the apparent number of columns
   *
   * @return the apparent number of columns
   */
  inline int getNCols() const {return NCols;};
  
  /**
   * Get the actual number of elements
   *
   * @return the actual number of elements
   */
  inline int getNElemsAct() const {return (data.size());};

  /**
   * Get the apparent number of elements
   *
   * @return the apparent number of elements
   */
  inline int getNElems() const {return (NRows*NCols);};

  /**
   * Get the equality tolerance
   *
   * @return the equality tolerance
   */
  inline const T getTol() const {return tol;};

  /**
   * Set the equality tolerance
   *
   * @param tol_in
   *            the equality tolerance
   */
  inline void putTol(T tol_in){tol=tol_in; return;};
  
  /**
   * Matrix style retrieval operator of a single element from the matrix by value
   *
   * @param i
   *        row of the element
   * @param j
   *        column of the element
   *
   * @return the matrix value
   */
  inline const T& operator()(int i, int j) const {
#ifdef __SURFMAT_ERR_CHECK__
    if(!((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols))){
      printf("NRowsAct=%d NRows=%d\n",NRowsAct,NRows);
      printf("NColsAct=%d NCols=%d\n",NColsAct,NCols);
      printf("data.size()=%d  NRowsAct*NColsAct=%d\n",data.size(),NRowsAct*NColsAct);
      assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));}
    if(!((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols))) {
      printf("ERROR: need 0<=i=%d<NRows=%d & 0<=j=%d<NCols=%d",
	     i,NRows,j,NCols);
      fflush(stdout);
      printf("\n");
      assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
    }
#endif
    return const_cast<T&> (data[jtoi[j]+i]);};

  
  /**
   * Matrix style retrieval operator of a single element from the matrix by reference
   *
   * @param i
   *        row of the element
   * @param j
   *        column of the element
   *
   * @return reference to the value
   */
#ifndef  __SURFMAT_ERR_CHECK__
  inline T& operator()(int i, int j) {
#else
  T& operator()(int i, int j) {
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    if(!((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols))){
      printf("i=%d NRows=%d j=%d NCols=%d",i,NRows,j,NCols);
      printf("\n"); //so can breakpoint on this line in debugger
      assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
    }
#endif
    return data[jtoi[j]+i];
  };


  /*  
  //vector style retrieve of pointer to element, for passing to BLAS & LAPACK convenience
  inline T* ptr(int k){
    //T* ptr(int k){
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    if(!((0<=k)&&(k<NRows*NCols))) {
      printf("ERROR: SurfMat T* ptr(k): need 0<k=%d, k<NRows*NCols=%d, NRows=%d, NCols=%d\n",k,NRows*NCols,NRows,NCols);
      assert((0<=k)&&(k<NRows*NCols));
    }
#endif
    return &(data[jtoi[k/NRows]+(k%NRows)]);
  };

  inline const T* ptr(int k) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=k)&&(k<NRows*NCols));
#endif
    return &(data[jtoi[k/NRows]+(k%NRows)]);
  };
  */

  /**
   * Matrix style retrieve of pointer to element, for passing to BLAS & LAPACK convenience
   *
   * @param i
   *        row of the element
   * @param j
   *        column of the element
   *
   * @return pointer to matrix element
   */
  inline T* ptr(int i, int j){
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
#endif
    return &(data[jtoi[j]+i]);
  };

  /**
   * Matrix style retrieve of pointer to element, for passing to BLAS & LAPACK convenience
   *
   * @param i
   *        row of the element
   * @param j
   *        column of the element
   *
   * @return pointer to matrix element
   */

  inline const T* ptr(int i, int j) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=i)&&(i<NRows)&&(0<=j)&&(j<NCols));
#endif
    return &(data[jtoi[j]+i]);
  };

  /**
   * Get the smallest value contained in this matrix, wraps T& minElem(T& val ,int& loc)
   *
   * @return the smallest value in the matrix
   */
  inline T minElem() const;
  
  /**
   * Gets the location of the smallest value contained in this matrix.
   *
   * If there is a tie it returns the location of the first one encountered,
   * wraps T& minElem(T& val ,int& loc)
   *
   * @return location of the smallest value
   */
  inline int iMinElem() const;
  
  /**
   * Get the smallest value contained in this matrix and its location in loc.
   *
   * If there is a tie it records the location of the first smallest value encountered
   *
   * @param val
   *        smallest value
   * @param loc
   *        smallest value location
   */
  inline void minElem(T& val, int& loc) const;
  
  /**
   * Gets the largest value contained in this matrix, wraps T& maxElem(T& val ,int& loc)
   *
   * @return largets value
   */
  inline T maxElem() const;
  
  /**
   * Gets the location of the largest value contained in this matrix
   *
   * If there is a tie it returns the location of the first one encountered, wraps T& maxElem(T& val ,int& loc)
   *
   * @return location of the largest value
   */
  inline int iMaxElem() const;

  /**
   * Get the largest value contained in this matrix and its location in loc.
   *
   * If there is a tie it records the location of the first largest value encountered
   *
   * @param val
   *        largest value
   * @param loc
   *        largest value location
   */
  inline void maxElem(T& val, int& loc) const;
  
  /**
   * Get the smallest and largest values contained in this matrix
   *
   * @param smallest
   *            smallest value
   * @param largest
   *            largest value
   */
  inline void minMaxElem(T& smallest, T& largest) const {
    int nelem=NRows*NCols;
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert(nelem>0);
#endif
    smallest=largest=data[0];
    for(int k=1; k<nelem; k++){
      //note I did not use an else if for the largest so that both if statements could be evaluated simultaneous by different processing units
      if(data[k]<smallest) smallest=data[k];
      if(data[k]>largest ) largest =data[k];
    }    
    return;
  };
  
  /**
   * replace the values in row "irow" of this matrix with the values stored in vecT_row
   *
   * This method will NOT expand the size of this matrix if you specify a row index larger
   * than NRows-1.
   *
   * @param vecT_row
   *        the vector to use
   * @param irow
   *        the row to replace
   */
  inline void putRows(const std::vector<T>& vecT_row, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=irow)&&(irow<NRows)&&(vecT_row.size()==NCols));
#endif
    for(int j=0; j<NCols; ++j) 
      data[jtoi[j]+irow]=vecT_row[j];
    return;    
  };

  /**
   * Replace the values in row "irow" of this matrix with the entries listed in string s
   *
   * @param s
   *        the string to use
   * @param irow
   *        the row to replace
   *
   * @return 0 if the string and row have the same number of entries and 1 otherwise
   */
  inline int putRows(const std::string& s, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=irow)&&(irow<NRows));
#endif
    std::istringstream is(s);
    int j;
    //T temp;
    for(j=0; j<NCols; ++j) {
      if(is.eof())
	break;
      //is >> temp; data[jtoi[j]+irow]=temp;
      //printf("irow=%d NRows=%d j=%d NCols=%d temp=%g\n",irow,NRows,j,NCols,temp);
      //printf("data[jtoi[j=%d]+(i=%d)]=%g\n",j,irow,data[jtoi[j]+irow]);
      is >> data[jtoi[j]+irow];
    }
    return(!((j==NCols)&&(is.eof())));
  };    

  /**
   *  Replace the values in row "irow" of this matrix with the values stored in "row"
   *
   *  This method will NOT expand the size of this matrix if you specify a row index
   *  larger than NRows-1.
   *
   *  @param row
   *        the row to use
   *  @param irow
   *        the row to replace
   */
  inline void putRows(const SurfMat<T>& row, int irow) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=irow)&&(irow<NRows)&&
	   (row.getNRows()==1)&&(row.getNCols()==NCols));
#endif
    for(int j=0; j<NCols; j++) 
      data[jtoi[j]+irow]=row(0,j);
    return;    
  };
  
  /**
   *  Replace the values in multiple rows of this matrix with the values stored in "rows"
   *
   *  This method will NOT expand the size of this matrix if you specify a row index
   *  larger than NRows-1.
   *
   *  @param rows
   *        the rows to use
   *  @param irows
   *        the rows to replace
   */
  inline void putRows(const SurfMat<T>& rows, SurfMat<int> irows) {
    int nrows_put=irows.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
    assert((irows.getNCols()==1)&&
	   (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=irows.minElem())&&(irows.maxElem()<NRows)&&
	   (rows.getNRows()==nrows_put)&&(rows.getNCols()==NCols));
#endif
    for(int j=0; j<NCols; j++)
      for(int k=0; k<nrows_put; k++)
	data[jtoi[j]+irows(k,0)]=rows(k,j);
    return;
  };
  
  /**
   * Replace the values in column "jcol" of this matrix with the values stored in vector "vecT_col"
   *
   * This method will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
   *
   * @param vectT_col
   *        vector to use
   * @param jcol
   *        column to replace
   */
  inline void putCols(const std::vector<T>& vecT_col, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=jcol)&&(jcol<NCols)&&(vecT_col.size()==NRows));
#endif
    for(int i=0; i<NRows; i++)
      data[jtoi[jcol]+i]=vecT_col[i];
    return;
  };

  /**
   * Replace the values in column "jcol" of this matrix with the entries listed in string s
   *
   * This method will NOT expand the size of this matrix if you specify a column index larger than NCols-1.
   *
   * @param s
   *        string to use
   * @param jcol
   *        column to replace
   *
   * @return 0 if the string and column have the same number of entries and 1 otherwise
   */
  inline int putCols(const std::string& s, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=jcol)&&(jcol<NCols));
#endif
    std::istringstream is(s);
    int i;
    for(i=0; i<NRows; ++i) {
      if(is.eof())
	break;
      is >> data[jtoi[jcol]+i];
    }    
    return(!((i==NRows)&&(is.eof()))); 
  };

  /**
   * Replace the values in column "jcol" of this matrix with the values stored in "col"
   *
   * This function will NOT expand the size of this matrix if you specify a column index
   * larger than NCols-1.

   * @param col
   *        column to use
   * @param jcol
   *        column to replace
   */
  inline void putCols(const SurfMat<T>& col, int jcol) {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    if(!((0<=jcol)&&(jcol<NCols)&&
	 (col.getNRows()==NRows)&&(col.getNCols()==1))){
      //assert((0<=jcol)&&(jcol<NCols)&&
      //(col.getNRows()==NRows)&&(col.getNCols()==1));
      assert(0<=jcol);
      assert(jcol<NCols);
      assert(col.getNRows()==NRows);
      assert(col.getNCols()==1);
    }
#endif
    for(int i=0; i<NRows; i++)
      data[jtoi[jcol]+i]=col(i,0);
    return;
  };

  /**
   * Replace the values in multiple columns of this matrix with the values stored in "cols"
   *
   * This function will NOT expand the size of this matrix if you specify a column index
   * larger than NCols-1.
   *
   * @param cols
   *        column(s) to use
   * @param jcol
   *        column to replace
   */
  inline void putCols(const SurfMat<T>& cols, SurfMat<int> jcols) {
    int ncols_put=jcols.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
    if(!((jcols.getNCols()==1)&&
	 (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct))) {
      printf("jcols.NCols=%d",jcols.getNCols());
      fflush(stdout);
      printf("\n");
    assert((jcols.getNCols()==1)&&
	   (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    }
    assert((0<=jcols.minElem())&&(jcols.maxElem()<NCols)&&
	   (cols.getNRows()==NRows)&&(cols.getNCols()==ncols_put));
#endif
    /*
    printf("putCols(): min(jcols)=%d max(jcols)=%d nelem_act=%d\n",jcols.minElem(),jcols.maxElem(),jcols.getNElemsAct());
    fflush(stdout);
    printf("  this mat\n         NRowsAct=%d NRows=%d\n         NColsAct=%d NCols=%d\n",NRowsAct,NRows,NColsAct,NCols);
    fflush(stdout);
    printf("         nelem_act=%d\n",data.size());
    fflush(stdout);
    printf("  cols mat\n         NRowsAct=%d NRows=%d\n         NColsAct=%d NCols=%d\n         nelem_act=%d\n",cols.NRowsAct,cols.NRows,cols.NColsAct,cols.NCols,cols.data.size());
    fflush(stdout);
    */

    for(int k=0; k<ncols_put; k++)
      for(int i=0; i<NRows; i++)
	data[jtoi[jcols(k,0)]+i]=cols(i,k); 
    return;
  };
  
  /**
   * Get one row (index stored in irow) of this matrix and return as a "row vector"
   *
   * @param result
   *        the resulting row vector
   * @param irow
   *        the row index to get
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return
   *        the resulting row vector
   */
  inline SurfMat<T>& getRows(SurfMat<T>& result, int irow, bool if_force=false) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    if(!((0<=irow)&&(irow<NRows))) {
      printf("irow=%d NRows=%d\n",irow,NRows); fflush(stdout);
      assert((0<=irow)&&(irow<NRows));
    }
#endif
    result.newSize(1,NCols,if_force);
    result.tol=tol;
    for(int j=0; j<NCols; j++) 
      result(0,j)=data[jtoi[j]+irow];
    return result;
  };
  
  /**
   * Get multiple rows (indices stored in matrix irows) of this matrix and return as a new matrix
   *
   * @param result
   *        the resulting matrix
   * @param irows
   *        the row(s) to get
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return
   *        the resulting matrix
   */
  inline SurfMat<T>& getRows(SurfMat<T>& result, SurfMat<int>& irows, bool if_force=false) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((irows.getNCols()==1)&&
	   (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=irows.minElem())&&(irows.maxElem()<NRows));
#endif
    int nrows_res=irows.getNRows();
    result.newSize(nrows_res,NCols,if_force);
    result.tol=tol;
    if(nrows_res>0) 
      for(int j=0; j<NCols; j++)
	for(int i=0; i<nrows_res; i++)
	  result(i,j)=data[jtoi[j]+irows(i,0)];	
    
    return result;
  };
  
  /**
   * Get one column (index stored in jcol) of this matrix and return as a vector
   *
   * @param result
   *        the resulting vector
   * @param jcol
   *        the column to get
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return
   *        the resulting vector
   */
  inline SurfMat<T>& getCols(SurfMat<T>& result, int jcol, bool if_force=false) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=jcol)&&(jcol<NCols));
#endif
    result.newSize(NRows,1,if_force);
    result.tol=tol;
    for (int i=0; i<NRows; i++) result(i,0)=data[jtoi[jcol]+i];
    return result;
  };
  
  /**
   * Get multiple columns (indices stored in matrix jcols) of this matrix and return as a new matrix
   *
   * @param result
   *        the resulting matrix
   * @param jcols
   *        the column(s) to get
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return
   *        the resulting column
   */
  inline SurfMat<T>& getCols(SurfMat<T>& result, SurfMat<int>& jcols, bool if_force=false) const {
#ifdef __SURFMAT_ERR_CHECK__
    assert((jcols.getNCols()==1)&&
	   (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
    assert((0<=jcols.minElem())&&(jcols.maxElem()<NCols));
#endif
    int ncols_res=jcols.getNRows();
    result.newSize(NRows,ncols_res,if_force); 
    result.tol=tol;
    if(ncols_res>0)
      for(int j=0; j<ncols_res; j++)
	for(int i=0; i<NRows; i++)
	  result(i,j)=data[jtoi[jcols(j,0)]+i];	
    
    return result;
  };
  
  /**
   * Gets a copy of the matrix excluding one row whose index is stored in irow
   *
   * @param result
   *        the resulting matrix
   * @param irow
   *        the row index to exclude
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return the resulting matrix
   */
  SurfMat<T>& excludeRows(SurfMat<T>& result, int irow, bool if_force=false);

  /**
   * Gets a copy of the matrix that excludes all the rows whose index is stored in irows
   *
   * @param result
   *        the resulting matrix
   * @param irows
   *        the row index(es) to exclude
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return the resulting matrix
   */
  ///returns a copy of the matrix that excludes all rows whose indices are stored in the matrix irows, this isn't an inline because for the typical use you will need to copy most of the matrix
  SurfMat<T>& excludeRows(SurfMat<T>& result, SurfMat<int>& irows, bool if_force=false);
  
  /**
   * Gets a copy of the matrix excluding one column whose index is stored in jcol
   *
   * @param result
   *        the resulting matrix
   * @param jcol
   *        the column index to exclude
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return the resulting matrix
   */
  SurfMat<T>& excludeCols(SurfMat<T>& result, int jcol, bool if_force=false);
  
  /**
   * Gets a copy of the matrix that excludes all the columns whose index is stored in jcols
   *
   * @param result
   *        the resulting matrix
   * @param jcols
   *        the column index(es) to exclude
   * @param if_force
   *        forces memory footprint to shrink if true
   * @return the resulting matrix
   */
  SurfMat<T>& excludeCols(SurfMat<T>& result, SurfMat<int>& jcols, bool if_force=false);
  
  /**
   * Compare 2 elements
   *
   * Used in the generic matrix quicksort
   *
   * @param ia
   *        index of first element to compare
   * @param ib
   *        index of 2nd element to compare
   * @return -1 if a-b<-tol, 0 if -tol<=a-b<=tol, +1 if tol<a-b
   */
  inline int compareElemAElemB(int ia, int ib) {
    T diff=data[ia]-data[ib];
    return(-(diff<-tol)+(tol<diff));
  }
  
  /**
   * Swap the ia-th and ib-th elements of the "vector"
   *
   * Used in the generic matrix quicksort
   *
   * @param ia
   *        index of first element to swap
   * @param ib
   *        index of 2nd element to swap
   */
  inline void swapElems(int ia, int ib) {
    T swap=data[ia];
    data[ia]=data[ib];
    data[ib]=swap;
    return;
  };
  
  /**
   * Ascending sort the elements of the matrix as if it were a vector
   *
   * @todo Could easily add an integer argument (+1 or -1) that indicates whether you
   *       want to sort into ascending or descending order, would need to add a field
   *       to the matrix though
   */
  inline void sortElems() {
    qsortElems(0,getNElems()-1);
    return;
  };
  

  /**
   * Quicksort elements of the matrix
   *
   * @param istart
   *        index of element to start at
   * @param istop
   *        index of element to stop at
   */
  void qsortElems(int istart, int istop);

  /**
   * Quicksort rows of the matrix
   *
   * @param istart
   *        index of row to start at
   * @param istop
   *        index of row to stop at
   */
  void qsortRows(int istart, int istop);

  /**
   * Quicksort columns of the matrix
   *
   * @param istart
   *        index of column to start at
   * @param istop
   *        index of column to stop at
   */
  void qsortCols(int istart, int istop);

/**
   * Performs an ascending unique sort of the elements.
   *
   * Eliminates duplicates, and reshapes it into a (shrunken, if
   * appropriate) vector
   */
  inline void uniqueElems() {
    int nelems=getNElems();
    reshape(nelems,1); // make sure the data we are about to quick
		       // sort is contiguous
    if(nelems>1) {
      sortElems();
      int i=0;
      for(int k=1; k<nelems; ++k) {
	if(data[i]!=data[k]) {
	  ++i;
	  if(i<k) data[i]=data[k];
	}              
      }                            
      ++i;
      if(i<nelems) {
	// i is the index of the element after the last one we are
	// keeping, i.e. the number of elements that we are keeping,
	// this is just record keeping of the number of unique
	// elements, it won't actually change the size of memory.
	reshape(i,1); 
      }
    }
    return;
  };
  
  /**
   * Compares 2 rows (a and b) starting with column 0
   *
   * Used in the generic matrix quicksort.
   *
   * @todo Could add a column order field to the matrix to add the capability to change
   *       the order in which the rows are sorted and whether or not you wanted to do an ascending
   *       or descending sort
   * @param irowa
   *        index of first row to compare
   * @param irowb
   *        index of 2nd row to compare
   * @return -1 if a-b<-tol, or +1 if tol<a-b; if -tol<=a-b<=tol it moves onto
   *         column 1, then 2 and so on, if every column of the 2 rows are equal
   *         then it returns 0.
   */
  inline int compareRowARowB(int irowa, int irowb) {
    T diff;
    int j=0;
    do{
      diff=data[jtoi[j]+irowa]-data[jtoi[j]+irowb];
      j++;
    }while((-tol<=diff)&&(diff<=tol)&&(j<NCols));
    return (-(diff<-tol)+(tol<diff));  //could multiply outermost () by +1/-1 for ascending/descending
  };
  
  /**
   * Swap rows a and b of the matrix
   *
   * used in the generic matrix quicksort
   *
   * @param irowa
   *        index of first row to swap
   * @param irowb
   *        index of 2nd row to swap
   */
  inline void swapRows(int irowa, int irowb) {
    T Tswap;
    for(int j=0; j<NCols; j++) {
      Tswap=data[jtoi[j]+irowa];
      data[jtoi[j]+irowa]=data[jtoi[j]+irowb];
      data[jtoi[j]+irowb]=Tswap;
    }
    return;
  };
  
  /**
   * Ascending sort the rows of the matrix (comparing from left to right)
   *
   * @todo You could easily add a MtxInt argument with 2 rows the first indicating order
   *       of columns being compared, the second being +1 or -1 to indicate to sort that
   *       row into ascending or descending order, would need to add a field to the matrix
   *       though
   */
  inline void sortRows() {
    qsortRows(0,NRows-1);
    return;
  };
  
  /**
   * Performs an ascending sort of the rows
   *
   * Eliminates duplicates, and shrinks the matrix if it is appropriate to do so
   */
  inline void uniqueRows() {
    if(NRows>1) {
      sortRows();
      int i,k,j;
      for(i=1; i<NRows; ++i)
	{
	  if(compareRowARowB(i,i-1) == 0)
	    break;
	}
      for(k=i+1; k<NRows; ++k) 
	if(compareRowARowB(k,i)) {
	  //only copy in the next row that's different from our "current marker" rather than move all rows one row earlier
	  for(j=0; j<NCols; ++j)
	    data[jtoi[j]+i]=data[jtoi[j]+k];
	  ++i; //make i the index of the row after the one we are keeping
	}   
      resize(i,NCols); //i is the index of the rows after the last one we are keeping, i.e. the number of rows that we are keeping, this is just record keeping of the number of unique rows, it won't actually change the size of memory.
    }
    return;
  };
  
  /**
   * Compares 2 columns starting with row 0
   *
   * Used in the generic matrix quicksort
   *
   * @todo could add a row order field to the matrix to add the capability to change the
   *       order in which the columns are sorted and whether or not you wanted to do an
   *       ascending or descending sort
   * @param jcola
   *        index of first column to compare
   * @param jcolb
   *        index of 2nd column to compare
   * @return -1 if a-b<-tol, or +1 if tol<a-b, if -tol<=a-b<=tol it moves onto the next row
   *          if every row of the 2 columns have -tol<=a-b<=tol then it returns 0
   */
  inline int compareColAColB(int jcola, int jcolb) {
    T diff;
    int i=0;
    do{
      diff=data[jtoi[jcola]+i]-data[jtoi[jcolb]+i];
      i++;
    }while((-tol<=diff)&&(diff<=tol)&&(i<NRows));
    return (-(diff<-tol)+(tol<diff));
  };

  /**
   * Swap columns a and b of the matrix
   *
   * Used in the generic matrix quicksort
   *
   * @param jcola
   *        index of first column to swap
   * @param jcolb
   *        index of 2nd column to swap
   */
  inline void swapCols(int jcola, int jcolb) {
    T Tswap;
    for(int i=0; i<NRows; i++) {
      Tswap=data[jtoi[jcola]+i];
      data[jtoi[jcola]+i]=data[jtoi[jcolb]+i];
      data[jtoi[jcolb]+i]=Tswap;
    }
    return;
  };
  
  /**
   * Ascending sort the columns of the matrix (comparing from row zero onward)
   *
   * @todo You could easily add a MtxInt argument with 2 columns the first indicating
   *       order of columns being compared (some rows can be omitted if you don't want
   *       to include them as sorting criteria), the second column being +1 or -1 to
   *       indicate to sort that row into ascending or descending order, would need to
   *       add a field to the matrix though
   */
  inline void sortCols() {
    qsortCols(0,NCols-1);
    return;
  };
  
  /**
   * Performs an ascending sort of the columns
   *
   * Eliminates duplicates, and shrinks the matrix if it is appropriate to do so
   */
  inline void uniqueCols() {
    if(NCols>1) {
      sortCols();
      int i,j,k,j2i,k2i;
      for(j=1; j<NCols; ++j)
	{
	  if(compareColAColB(j,j-1) == 0)
	    {
	      break;
	    }
	}
      j2i = jtoi[j];
      for(k=j+1; k<NCols; ++k) 
	if(compareColAColB(k,j-1)) {
	  //only copy in the next coll that's different from our "current marker" rather than move all rows one row earlier
	  k2i = jtoi[k];
	  for(int i=0; i<NRows; ++i)
	    {
	      
	      data[j2i+i]=data[k2i+i];
	    }
	  ++j; //make i the index of the row after the one we are keeping
	  j2i = jtoi[j];
	}   

      reshape(NRows,j); //not an error, want contiguous memory chopping, reshape will call reshape2, resize would call resize2 which would call reshape2  
    }
    return;
  };

private:

  /**
   * Enlarge or shrink the matrix without preserving it's contents
   *
   * This is here for speed (resize2 and reshape2 are slower because they involve copying data,
   * inherently done by realloc when it can't grab extra space at the end of the array)
   *
   * newSize2 should not be called by anything except newSize
   *
   * @param nrows_new
   *            the number of rows to change to
   * @param ncols_new
   *            the number of columns to change to
   * @param if_force
   *            forces memory footprint to shrink if true
   */
  void newSize2(int nrows_new, int ncols_new, bool if_force=false);

  /**
   * Enlarge or shrink the matrix while keeping contigous in memory elements contiguous in memory
   *
   * Reshape2() should not be called by anything but reshape() and resize2()
   *
   * @param nrows_new
   *        the number of rows to change to
   * @param ncols_new
   *        the number of columns to change to
   * @param if_force
   *        forces memory footprint to shrink if true
   */
  void reshape2(int nrows_new, int ncols_new, bool if_force=false);

  /**
   * Enlarge or shrink the matrix while keeping contigous in memory elements contiguous in memory
   *
   * Resize2 should not be called by anything but resize
   *
   */
  void resize2(int nrows_new, int ncols_new, bool if_force=false);
  int NRowsAct; /**< number of allocated rows, passed to LAPACK as LDA */
  int NColsAct; /**< number of allocated columns */
  int NRows; /**< number of used rows, passed to LAPACK as the number of rows */
  int NCols; /**< number of used columns */
  std::vector<T>  data;
  std::vector<int> jtoi;

  /**< An inequaltiy tolerance for equality checking. Should be 0 for integers
   *
   * @todo WARNING TO WHOEVER FINISHES IMPLEMENTING THIS, you will have to decide what to do with tol
   *       when you clear(), newsize(), reshape(), or resize() the matrix, you will likely need to modify
   *       the contructors
   */
  T   tol;

  
  //this shouldn't be necessary but do it for good measure
  friend void mtxqsort(int istart, int istop,
                       int (* compare_a_b)(int ia, int ib),
                       void (* swap_a_b)(int ia, int ib));

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class SurfPoint data
  template<class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif
};



//typical matrix constructor
template< typename T >
SurfMat<T>::SurfMat(int nrows_in, int ncols_in)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_in>=0)&&(ncols_in>=0));
#endif
  tol=static_cast<T>(0);
  NRowsAct=NColsAct=NRows=NCols=0;
  if((nrows_in>0)&&(ncols_in>0)) {
    NRowsAct=NRows=nrows_in;
    NColsAct=NCols=ncols_in;
    data.resize(NRowsAct*NColsAct);
    jtoi.resize(NColsAct);
    for(int j=0, J=0; j<NCols; ++j, J+=NRowsAct)
      jtoi[j]=J;
  }

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif

  return;
}

//copy constructor, make a deep copy
template< typename T >
SurfMat<T>::SurfMat(const SurfMat<T>& other){
  tol=other.tol;
  NRowsAct=other.NRowsAct;
  NRows=other.NRows;
  NColsAct=other.NColsAct;
  NCols=other.NCols;
  if(!((NRowsAct>0)&&(NRows>0)&&(NColsAct>0)&&(NCols>0))) {
    return;
    //printf("stop me\n");
    //assert((NRowsAct>0)&&(NRows>0)&&(NColsAct>0)(NCols>0)); //for debug only
  }

  data=other.data;
  jtoi.resize(NColsAct);
  for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
    jtoi[j]=J;

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif

  return;
}

//copy from vector constructor
template< typename T >
SurfMat<T>::SurfMat(const std::vector<T>& vecT, int nrows_in, int ncols_in)
{
  tol=static_cast<T>(0);
  NRowsAct=NRows=NColsAct=NCols=0;
  
  int size_vecT=vecT.size();

  if(nrows_in==-99999) 
    nrows_in=size_vecT;

#ifdef __SURFMAT_ERR_CHECK__
  assert((nrows_in>=0)&&(ncols_in>=0));
#endif

  if((nrows_in>0)&&(ncols_in>0)) {
    NRowsAct=NRows=nrows_in;
    NColsAct=NCols=ncols_in;
    for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
      jtoi[j]=J;    
  }  
  
  int n_to_copy=NRowsAct*NColsAct;
  if(n_to_copy>size_vecT)
    n_to_copy=size_vecT;
  
  for(int k=0; k<n_to_copy; ++k)
    data[k]=vecT[k];

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif

  return;
}


//deallocate the memory
template< typename T >
void SurfMat<T>::clear()
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  if(!(((NRows==0)&&(NCols==0)&&(NRowsAct==0)&&(NColsAct==0))||
       (NRowsAct && NColsAct))) {
    printf("NRows=%d NCols=%d NRowsAct=%d NColsAct=%d",NRows,NCols,NRowsAct,NColsAct);
    assert(((NRows==0)&&(NCols==0)&&(NRowsAct==0)&&(NColsAct==0))||
	   (NRowsAct && NColsAct));
  }
#endif
  if(NRowsAct) {
    jtoi.clear();
    data.clear();
    NRowsAct=NColsAct=NRows=NCols=0;
  }
  return;
}

//enlarge or shrink the matrix without preserving it's contents, this is here for speed (resize2 and reshape2 are slower because they involve copying data, inherently done by realloc when it can't grab extra space at the end of the array) 
template< typename T >
void SurfMat<T>::newSize2(int nrows_new, int ncols_new, bool if_force)
{

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif

  int nelem_act=data.size();
  int nelem_new=nrows_new*ncols_new;

  if((NRows==nrows_new)&&(NCols==ncols_new)&&
     ((if_force==false)||
      ((if_force==true)&&(nelem_act==nelem_new)&&
       (NRowsAct==nrows_new)&&(NColsAct==ncols_new)))) {
    //neither the actual or apparent size need to change
    return;
  }

  if(((if_force==false)&&(nelem_act>=nelem_new)&&
      (NRowsAct>=nrows_new)&&(NColsAct>=ncols_new)) ||
     ((if_force==true )&&(nelem_act==nelem_new)&&
      (NRowsAct==nrows_new)&&(NColsAct==ncols_new))) {
    //the actual size doesn't need to change, so just change the apparent size 
    NRows=nrows_new;
    NCols=ncols_new;
    return;
  }

  if(nelem_new==0){
    if(if_force==true) 
      clear(); //only deallocate if they force you to, clear() is the way to deallocate when they know they want to deallocate
    else
      NRows=NCols=0;
    return;
  }

  if(((if_force==false)&&(nelem_new<=nelem_act))||
     ((if_force==true )&&(nelem_new==nelem_act))) {
    //don't need to change the number of elements so don't resize data
  }
  else{
    data.resize(nelem_new);
    nelem_act=nelem_new;
  }

#ifdef __SURFMAT_ZERO_MEM__
  for(int j=0; j<nelem_act; ++j) 
    data[j]=0;   
#endif

  NRowsAct=NRows=nrows_new;  

  NCols=nelem_act/NRowsAct;  //what NColsAct needs to be
  //NCols is the number of FULL columns with nrows_new rows that can fit 
  //into nelem_act 
  //if if_force==true  then nelem_act=nelem_new
  //if if_force==false then nelem_act>=nelem_new
  
  if(NColsAct!=NCols) {
    //we need to change the number of Columns
    NColsAct=NCols;
    jtoi.resize(NColsAct);
    NCols=ncols_new;
  }

#ifdef __SURFMAT_ERR_CHECK__
  assert((nelem_act==data.size())&&
	 ((if_force==false)&&
	  (nelem_act>=NRowsAct*NColsAct)&&(nelem_act>=nelem_new)&&
	  (NRowsAct>=nrows_new)&&(NColsAct>=ncols_new))||
	 ((if_force==true )&&
	  (nelem_act==NRowsAct*NColsAct)&&(nelem_act==nelem_new)&&
	  (NRowsAct==nrows_new)&&(NColsAct==ncols_new)));
#endif

  for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
    jtoi[j]=J;

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif

  return;
}

//enlarge or shrink the matrix while keeping contigous in memory elements contiguous in memory: reshape2() should not be called by anything but reshape() and resize2()
template< typename T >
void SurfMat<T>::reshape2(int nrows_new, int ncols_new, bool if_force)
{
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif

  int nelem_act=data.size();
  int nelem_new=nrows_new*ncols_new;

  if((NRows==nrows_new)&&(NCols==ncols_new)&&
     ((if_force==false)||
      ((if_force==true)&&(nelem_act==nelem_new)&&
       (NRowsAct==nrows_new)&&(NColsAct==ncols_new)))) {
    //neither the actual or apparent size need to change
    return;
  }

  if((NRows==nrows_new)&&
     (((if_force==false)&&(NColsAct>=ncols_new))||
      ((if_force==true)&&(NRowsAct==nrows_new)&&(NColsAct==ncols_new)))) {
    //only the apparent number of columns needs to change
    NCols=ncols_new;
    return;
  }

  if(nelem_new==0){
    if(if_force==true) 
      clear(); //only deallocate if they force you to, clear() is the way to deallocate when they know they want to deallocate
    else
      NRows=NCols=0;
    return;
  }



  int k;
  if(NRowsAct!=NRows) {
    //we know that NRows < NRowsAct, so we need to copy the data towards the beginning so that it's contiguous in memory before we do the actual reshaping, 
    k=NRows;
    for(int j=1; j<NCols; ++j)
      for(int i=0; i<NRows; ++i, ++k) 
	data[k]=data[jtoi[j]+i];
    //NRows=NRowsAct; //unnesecary
  }
  else{
    k=NRows*NCols;
  }
  
  if(((if_force==false)&&(nelem_act<nelem_new))||
     ((if_force==true )&&(nelem_act!=nelem_new))) {
    data.resize(nelem_new);
    nelem_act=nelem_new;
  }

  NRowsAct=NRows=nrows_new;
  NCols=nelem_act/NRowsAct; //the number of columns that we need
  //if if_force==true then NCols=ncols_new
  
  if(NColsAct!=NCols) {
    NColsAct=NCols;
    jtoi.resize(NColsAct);
  }
  NCols=ncols_new;

  for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
    jtoi[j]=J;
  return;
  
#ifdef __SURFMAT_ZERO_MEM__
  for(; k<nelem_act; ++k) 
    data[k]=0;
#endif      
  
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  return;
}

//enlarge or shrink the matrix while adding "zeros" after the last row and/or last column and/or chopping off the rows and/or columns after the newly requested last row and/or column, behavior "should be" the same as Teuchos' Serial Dense Matrix's reshape function, this function should only be called by resize()
template< typename T >
void SurfMat<T>::resize2(int nrows_new, int ncols_new, bool if_force) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert((nrows_new>=0)&&(ncols_new>=0));
#endif

  int nelem_act=data.size();
  int nelem_new=nrows_new*ncols_new;

  if((NRows==nrows_new)&&(NCols==ncols_new)&&
     ((if_force==false)||
      ((if_force==true)&&(nelem_act==nelem_new)&&
       (NRowsAct==nrows_new)&&(NColsAct==ncols_new)))) {
    //neither the actual or apparent size need to change
    return;
  }

  if(((if_force==false)&&(nelem_act>=nelem_new)&&
      (NRowsAct>=nrows_new)&&(NColsAct>=ncols_new)) ||
     ((if_force==true )&&(nelem_act==nelem_new)&&
      (NRowsAct==nrows_new)&&(NColsAct==ncols_new))) {
    //the actual size doesn't need to change, so just change the apparent size 
    NRows=nrows_new;
    NCols=ncols_new;
    return;
  }

  if(nelem_new==0){
    if(if_force==true) 
      clear(); //only deallocate if they force you to, clear() is the way to deallocate when they know they want to deallocate
    else
      NRows=NCols=0;
    return;
  }

  
  if(NRowsAct==nrows_new) {
    //other than to add more columns there is no reason for us to copy the data and if we need to data.resize() will take care of it for us
    NRows=nrows_new;
    if(((if_force==false)&&(nelem_act<nelem_new))||
       ((if_force==true )&&(nelem_act!=nelem_new))) {
      data.resize(nelem_new);
      nelem_act=nelem_new;
    }
    int ncols_act_new=nelem_act/NRowsAct;
    if(NColsAct!=ncols_act_new) {
      jtoi.resize(ncols_act_new);
    }
    NColsAct=ncols_act_new;
    NCols=ncols_new;

    for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
      jtoi[j]=J;
  }
  else if(NRowsAct<nrows_new) {
    //if the matrix we need has more elements than the matrix we have, then we
    //will need to copy the data to a new array
    //if the matrix we need has fewer element than the matrix we have, we can copy data to later parts of the array we already have

    if(nelem_act<nelem_new) {
      //we need to copy the data to a new (larger) array

      std::vector<T> newdata(nelem_new);
      std::vector<int> newjtoi(ncols_new);

      for(int j=0, J=0; j<ncols_new; ++j, J+=nrows_new)
	newjtoi[j]=J;

      int min_cols=(NCols<=ncols_new)?NCols:ncols_new;
      int min_rows=(NRows<=nrows_new)?NRows:nrows_new;
      for(int j=0; j<min_cols; ++j)
	for(int i=0; i<min_rows; ++i)
	  newdata[newjtoi[j]+i]=data[jtoi[j]+i];

      data.swap(newdata);
      jtoi.swap(newjtoi);
      NRowsAct=NRows=nrows_new;
      NColsAct=NCols=ncols_new;
    }
    else{
      //we can keep the data in the current array 

      int min_cols=(NCols<=ncols_new)?NCols:ncols_new;
      int min_rows=(NRows<=nrows_new)?NRows:nrows_new;
      if((if_force==true)&&(nelem_act!=nelem_new)) {
	//truncate data to its new size
	data.resize(nelem_new);
	nelem_act=nelem_new;
	for(int j=0, J=0; j<min_cols; ++j, J+=NRowsAct)
	  jtoi[j]=J;
      }

      int ncols_act_new=nelem_act/nrows_new;
      std::vector<int> newjtoi(ncols_act_new);      
      for(int j=0, J=0; j<ncols_act_new; ++j, J+=nrows_new)
	newjtoi[j]=J;

      for(int j=min_cols-1; j>=0; --j)
	for(int i=min_rows-1; i>=0; --i)
	  data[newjtoi[j]+i]=data[jtoi[j]+i];

      jtoi.swap(newjtoi);
      NRowsAct=NRows=nrows_new;
      NColsAct=ncols_act_new;
      NCols=ncols_new;
    }
  }
  else if(NRowsAct>nrows_new){

    int nelem_act_new=nelem_act;
    int nrows_act_new=NRowsAct;
    if(nelem_new>nelem_act) {
      data.resize(nelem_new);
      nelem_act_new=nelem_new;
      nrows_act_new=nrows_new;
    }
    else if((if_force==false)&&(NRowsAct*ncols_new>nelem_act)) {
      nrows_act_new=nelem_act/ncols_new;
      //I don't think that you should be able to get here, unless changes were made since KRD wrote it on 2011.05.09;
#ifdef __SURFMAT_ERR_CHECK__
      assert(nrows_act_new>=nrows_new);
#endif
    }
    else if(if_force==true) {
      nelem_act_new=nelem_new;
      nrows_act_new=nrows_new;
    }

    int ncols_act_new=nelem_act_new/nrows_act_new;

    std::vector<int> newjtoi(ncols_act_new);
    for(int j=0, J=0; j<ncols_act_new; ++j, J+=nrows_act_new) 
      newjtoi[j]=J;

    //copy towards the front of the array
    for(int j=1; j<ncols_new; ++j)
      for(int i=0; i<nrows_new; ++i)
	data[newjtoi[j]+i]=data[jtoi[j]+i];
    
    if((if_force==true)&&(nelem_new!=nelem_act))
      data.resize(nelem_new);
    
    jtoi.swap(newjtoi);
    NRowsAct=nrows_act_new;
    NRows=nrows_new;
    NColsAct=ncols_act_new;
    NCols=ncols_new;
  }
#ifdef __SURFMAT_ERR_CHECK__
  else{
    printf("it should not have been possible for you to get here unless NRowsAct OR nrows_new was a NaN\n");
    assert(false);
  }
#endif

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  return;
}

///Perform a deep copy
template< typename T >
SurfMat<T>& SurfMat<T>::copy(const SurfMat<T>& other, bool if_force){
  
  int nrows_new=other.NRows, ncols_new=other.NCols;
  int nelem_act=data.size(), //nelem=NRows*NCols, 
    nelem_new=nrows_new*ncols_new;

  //printf("copy:\n  nrows_new=%d ncols_new=%d nelem_new=%d\n  NRowsAct=%d NRows=%d\n NColsAct=%d NCols=%d\n  nelem_act=%d\n",nrows_new,ncols_new,nelem_new,NRowsAct,NRows,NColsAct,NCols,nelem_act);

  if(((if_force==false)&&(nelem_act>=nelem_new)&&
      (NRowsAct>=nrows_new)&&(NColsAct>=ncols_new))||
     ((if_force==true )&&(nelem_act==nelem_new)&&
      (NRowsAct==nrows_new)&&(NColsAct==ncols_new))) {
    //don't resize anything
    //printf("*if path 1\n");
  }
  else{
    //printf("*if path 2\n");  

    bool if_size_changed=false;
    if(((if_force==false)&&(nelem_new>nelem_act))||
       ((if_force==true )&&(nelem_new!=nelem_act))) {
      if_size_changed=true;
      data.resize(nelem_new);
      nelem_act=nelem_new; 
      //printf("*if path 2.1\n");
    }

    if(((if_force==false)&&(ncols_new>NColsAct))||
       ((if_force==true )||(ncols_new!=NColsAct))) {
      if_size_changed=true;
      NColsAct=ncols_new;
      jtoi.resize(NColsAct);
      //printf("*if path 2.2\n");
    }

    NRowsAct=nelem_act/NColsAct;
#ifdef __SURFMAT_ERR_CHECK__
    assert(NRowsAct>=nrows_new);
#endif

    if(if_size_changed==true) 
      for(int j=0, J=0; j<NColsAct; ++j, J+=NRowsAct)
	jtoi[j]=J;
  }  

 NRows=nrows_new;
 NCols=ncols_new;

 tol=other.tol;
  
 for(int j=0; j<NCols; ++j)
   for(int i=0; i<NRows; ++i)
     data[jtoi[j]+i]=other(i,j);

#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  //printf("exiting copy\n");
  return *this;
}


///returns a copy of the matrix excluding 1 row whose index is stored in irow, this isn't an inline because you will need to copy all but 1 row
template< typename T >
SurfMat<T>& SurfMat<T>::excludeRows(SurfMat<T>& result, int irow, bool if_force) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert((0<=irow)&&(irow<NRows));
#endif
  if(NRows==1) {
    if(if_force)
      result.clear();
    else{
      result.NRows=NRows;
      result.NCols=NCols;
      result.tol=tol;
    }
  }
  else{
    result.newSize(NRows-1,NCols,if_force);
    result.tol=tol;
    int isrc, ikeep, j;
    for(j=0; j<NCols; j++) {      
      for(isrc=ikeep=0; isrc<irow; ++isrc, ++ikeep)
	result(ikeep,j)=data[jtoi[j]+isrc];
      isrc=irow+1;
      for(;isrc<NRows; ++isrc, ++ikeep)
	result(ikeep,j)=data[jtoi[j]+isrc];
    }
  }

#ifdef __SURFMAT_ERR_CHECK__
  assert((result.data.size()>=result.NRowsAct*result.NColsAct)&&
	 (result.NRowsAct>=result.NRows)&&(result.NColsAct>=result.NCols));
#endif

  return result;
}

///returns a copy of the matrix excluding all rows whose indices are stored in matrix irows, this isn't an inline because for the typical use you will need to copy most of the matrix
template< typename T >
SurfMat<T>& SurfMat<T>::excludeRows(SurfMat<T>& result, SurfMat<int>& irows, bool if_force) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((irows.getNCols()==1)&&
	 (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  int j;
  int nexclude=irows.getNRows();
  if(nexclude<1) {
    //the list of rows to exclude is empty so copy over the whole matrix
    result.copy(*this,if_force);
  }
  else{
    irows.uniqueElems(); //sort the rows to exclude into ascending order and eliminate duplicate listings
    nexclude=irows.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=irows(0,0))&&(irows(nexclude-1,0)<NRows));
#endif
    if(nexclude==NRows) {
      //the user wants us to eliminate _all_ rows
      if(if_force==true)
	result.clear();
      else {
	result.tol=tol;
	result.NRows=0;
	result.NCols=0;
      }
    }
    else{
      //the user wants us to eliminate some but not all rows
      result.newSize(NRows-nexclude,NCols,if_force);
      result.tol=tol;
      int iexclude, ikeep, isrc;
      for(j=0;j<NCols;j++) {
	iexclude=ikeep=isrc=0;
	while(isrc<NRows) {
	  if(iexclude<nexclude) {
	    for(;isrc<irows(iexclude,0); ++isrc, ++ikeep)
	      result(ikeep,j)=data[jtoi[j]+isrc];
	    //at this point isrc=irows(iexclude,0)
	    ++iexclude;
	    ++isrc;}
	  else{
	    for(;isrc<NRows; ++isrc, ++ikeep)
	      result(ikeep,j)=data[jtoi[j]+isrc];
	    //at this point isrc=NRows and the while loop will terminate
	  }
	}//while loop terminates
      }//do the same thing with the next column  
    }
  }
#ifdef __SURFMAT_ERR_CHECK__
  assert((result.data.size()>=result.NRowsAct*result.NColsAct)&&
	 (result.NRowsAct>=result.NRows)&&(result.NColsAct>=result.NCols));
#endif
  return result;
}

///returns a copy of the matrix excluding 1 column whose index is stored in jcol, this isn't an inline because you will need to copy all but 1 column
template< typename T >
SurfMat<T>& SurfMat<T>::excludeCols(SurfMat<T>& result, int jcol, bool if_force) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert((0<=jcol)&&(jcol<NCols));
#endif
  if(NCols==1) {
    if(if_force==true)
      result.clear();
    else{
      result.NRows=0;
      result.NCols=0;
      result.tol=tol;
    }
  }
  else{
    result.newSize(NRows,NCols-1);
    result.tol=tol;
    int jsrc, jkeep, i;
    for(jsrc=jkeep=0; jsrc<jcol; jsrc++, jkeep++)
      for(i=0; i<NRows; i++)
	result(i,jkeep)=data[jtoi[jsrc]+i];
    ++jsrc;
    for(; jsrc<NRows; jsrc++, jkeep++)
      for(i=0; i<NRows; ++i)
	result(i,jkeep)=data[jtoi[jsrc]+i];
  }
#ifdef __SURFMAT_ERR_CHECK__
  assert((result.data.size()>=result.NRowsAct*result.NColsAct)&&
	 (result.NRowsAct>=result.NRows)&&(result.NColsAct>=result.NCols));
#endif
  return result;
}

///return a copy of the matrix that excludes all columns whose indices are stored in the matrix jcols, this isn't an inline because for the typical use case you will need to copy most of the matrix
template< typename T >
SurfMat<T>& SurfMat<T>::excludeCols(SurfMat<T>& result, SurfMat<int>& jcols, bool if_force) {
#ifdef __SURFMAT_ERR_CHECK__
  assert((jcols.getNCols()==1)&&
	 (data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
#endif
  int nexclude=jcols.getNRows();
  int i;
  if(nexclude<1) 
    result.copy(*this, if_force);
  else{
    jcols.uniqueElems();
    nexclude=jcols.getNRows();
#ifdef __SURFMAT_ERR_CHECK__
    assert((0<=jcols(0,0))&&(jcols(nexclude-1,0)<NCols));
#endif
    if(nexclude==NCols) {
      //the user wants us to eliminate _all_ columns
      if(if_force)
	result.clear();
      else{
	result.NRows=0;
	result.NCols=0;
	result.tol=tol;
      }
    }
    else {
      //the user wants us to eliminate some but not all columns
      result.newSize(NRows,NCols-nexclude);
      result.tol=tol;
      int jexclude, jkeep, jsrc;
      jexclude=jkeep=jsrc=0;
      while(jsrc<NCols) {
	if(jexclude<nexclude) {
	  for(;jsrc<jcols(jexclude,0); ++jsrc, ++jkeep)
	    for(i=0;i<NRows;++i)
	      result(i,jkeep)=data[jtoi[jsrc]+i];
	  //at this point jsrc=jrows(jexclude)
	  ++jexclude;
	  ++jsrc;
	}
	else{
	  for(;jsrc<NCols; ++jsrc, ++jkeep)
	    for(i=0; i<NRows; ++i)
	      result(i,jkeep)=data[jtoi[jsrc]+i];
	  //at this point jsrc=NCols and the while loop will terminate
	}
      }//while loop terminates
    }
  }
#ifdef __SURFMAT_ERR_CHECK__
  assert((result.data.size()>=result.NRowsAct*result.NColsAct)&&
	 (result.NRowsAct>=result.NRows)&&(result.NColsAct>=result.NCols));
#endif
  return result;
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
	i++;
      while((j > istart) && (compareElemAElemB(istart,j)<0))
	j--;
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
	i++;
      while((j > istart) && (compareRowARowB(istart,j)<0))
	j--;
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
	i++;
      while((j > istart) && (compareColAColB(istart,j)<0))
	j--;
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

///returns the smallest value contained in this matrix, wraps T& minElem(T& val ,int& loc)
template< typename T >
inline T SurfMat<T>::minElem() const {
  T val;
  int loc;
  minElem(val,loc);
  return val;
}

///returns the location of the smallest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& minElem(T& val ,int& loc)
template< typename T >
inline int SurfMat<T>::iMinElem() const {
  T val;
  int loc;
  minElem(val,loc);
  return loc;
}

///returns the smallest value contained in this matrix and its location in loc, if there is a tie it records the location of the first smallest value encountered
template< typename T >
inline void SurfMat<T>::minElem(T& val, int& loc) const {
  int nelem=NRows*NCols; 
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert(nelem>0);
#endif
  val=data[0];
  loc=0;
  for(int k=1; k<nelem; k++)
    if(data[k]<val){
      val=data[k];
      loc=k;}
  return;
}


///returns the largest value contained in this matrix, wraps T& maxElem(T& val ,int& loc)
template< typename T >
inline T SurfMat<T>::maxElem() const {
  T val;
  int loc;
  maxElem(val,loc);
  return val;
}

///returns the location of the largest value contained in this matrix, if there is a tie it returns the location of the first one encountered, wraps T& maxElem(T& val ,int& loc)
template< typename T >
inline int SurfMat<T>::iMaxElem() const {
  T val;
  int loc;
  maxElem(val,loc);
  return loc;
}

///returns the largest value contained in this matrix and its location in loc, if there is a tie it records the location of the first largest value encountered
template< typename T >
inline void SurfMat<T>::maxElem(T& val, int& loc) const {
  int nelem=NRows*NCols; 
#ifdef __SURFMAT_ERR_CHECK__
  assert((data.size()>=NRowsAct*NColsAct)&&(NRowsAct>=NRows)&&(NColsAct>=NCols)&&(jtoi.size()>=NColsAct));
  assert(nelem>0);
#endif
  val=data[0];
  loc=0;
  for(int k=1; k<nelem; k++)
    if(data[k]>val){
      val=data[k];
      loc=k;}
  return;
}
} // end namespace nkm

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template<typename T>
template<class Archive>
void nkm::SurfMat<T>::serialize(Archive & archive, 
				const unsigned int version)
{
  archive & NRowsAct; 
  archive & NColsAct; 
  archive & NRows; 
  archive & NCols; 
  archive & data;
  archive & jtoi;
  archive & tol; 
}
#endif

#endif // __SURFMAT_HPP__
