/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef SURFPACK_MATRIX
#define SURFPACK_MATRIX

#include "surfpack_system_headers.h"

/// A SurfpackMatrix uses an STL vector-- where all elements reside in a 
/// contiguous block of memory-- to represent a two-dimensional matrix.
/// Overloading of operator[], together with the use of helper class
/// SurfpackVector, allows for the familiar, convenient, intuitive double-
/// subscripting.  A[i][j] is the element in the i-th row and j-th column of A.
/// In Surfpack, the data in such matrices are frequently passed through to
/// Fortran APIs (for use with BLAS and LAPACK routines).  Fortran expects the
/// data to be organized with elements adjacent to other elements in the column.
/// In C, on the other hand, doubly-subscripted arrays are organized so that
/// elements of the same row are adjacent to each other in memory.  The client
/// of SurfpackMatrix/Vector must specify at the time of object construction
/// whether the data are to be organized after the manner of C or Fortran, via
/// the constructor paramater "for_fortran."  After that, the client may trust
/// that A[i][j] means row-i and column-j, without worrying about the underlying
/// organization of the data.
template< typename T >
class SurfpackMatrix
{
public:
  /// Creates a n_rows x n_cols matrix
  SurfpackMatrix(unsigned n_rows, unsigned n_cols, bool for_fortran = true);
  SurfpackMatrix(unsigned n_rows, unsigned n_cols, const std::string& vals, bool for_fortran = true);
  
  /// Defaults to a 1x1 matrix for use with Fortran
  SurfpackMatrix(bool for_fortran = true);

  /// Copy constructor.  Makes deep copy.
  SurfpackMatrix(const SurfpackMatrix<T>& other);

  /// Resize the matrix. 
  ///\todo Deprecate this.  Use resize instead.  Or make it work like Matlab's.
  void reshape(unsigned n_rows, unsigned n_cols);

  /// Resize the matrix; retain previous values where possible
  void resize(unsigned n_rows, unsigned n_cols);

  /// Elements separated by spaces; rows delimited by newlines.
  const std::string asString() const;

  /// Returns the elements of a matrix as a list, in the order they appear in
  /// memory.
  const std::string asArrayString() const;

  /// retrieve a single element from the matrix by value
  const T& operator()(unsigned r, unsigned c) const;

  /// retrieve a single element from the matrix by reference
  T& operator()(unsigned r, unsigned c);

  /// Assignment operator.  Performs deep copy.
  SurfpackMatrix<T>& operator=(const SurfpackMatrix<T>& other);

  /// Assignment operator.  Performs deep copy.
  SurfpackMatrix<T>& operator=(const std::vector< std::vector<T> >);
  
  /// Sum of the elements of the matrix.  Data type T should support operator+.
  T sum() const;

  /// Number of rows in the matrix.
  unsigned getNRows(char transpose = 'N') const;

  /// Number of columns in the matrix.
  unsigned getNCols(char transpose = 'N') const;
private:
  /// True iff matrix can be safely passed to a Fortran API.
  bool forFortran;

  /// Number of rows in the matrix.
  unsigned nRows;

  /// Number of columns in the matrix.
  unsigned nCols;

  /// Contigous block of memory stores elements of the matrix.
  std::vector< T > rawData;

  /// Row of matrix most recently references by operator[]
  //mutable SurfpackVector< T > oneRow;

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
  // allow serializers access to private data
  friend class boost::serialization::access;
  /// serializer for derived class matrix data
  template< class Archive> 
  void serialize(Archive & archive, const unsigned int version);
#endif


};

// Implementation of functions
template< typename T >
SurfpackMatrix< T >::SurfpackMatrix(unsigned n_rows, unsigned n_cols,
  bool for_fortran)
  : forFortran(for_fortran), nRows(n_rows), nCols(n_cols)
{
  rawData.resize(nRows*nCols);
}

template< typename T >
SurfpackMatrix< T >::SurfpackMatrix(bool for_fortran)
  : forFortran(for_fortran), nRows(1), nCols(1) 
{
  rawData.resize(1);
}

template< typename T >
SurfpackMatrix< T >::SurfpackMatrix(unsigned n_rows, unsigned n_cols,
  const std::string& vals, bool for_fortran)
  : forFortran(for_fortran), nRows(n_rows), nCols(n_cols)
{
  rawData.resize(nRows*nCols);
  std::istringstream is(vals);
  for (unsigned i = 0; i < rawData.size(); i++) {
    is >> rawData[i];
  }
}

template< typename T >
SurfpackMatrix< T>::SurfpackMatrix(const SurfpackMatrix<T>& other)
  : forFortran(other.forFortran), nRows(other.nRows), nCols(other.nCols),
  rawData(other.rawData)
{

}


template< typename T >
void SurfpackMatrix< T >::reshape(unsigned n_rows, unsigned n_cols)
{
  nRows = n_rows;
  nCols = n_cols;
  rawData.resize(n_rows*n_cols);
}

template< typename T >
void SurfpackMatrix< T >::resize(unsigned n_rows, unsigned n_cols)
{
  unsigned old_rows = nRows;
  unsigned old_cols = nCols;
  nRows = n_rows;
  nCols = n_cols;
  std::vector< T > old_data = rawData;
  // Now go back and make sure that all of the elements in the new
  // matrix that were also present in the old matrix retain their values
  //unsigned row_max = min(old_rows,nRows);
  unsigned row_max = (old_rows < nRows) ? old_rows : nRows;
  //unsigned col_max = min(old_cols,nCols);
  unsigned col_max = (old_cols < nCols) ? old_cols : nCols;
  rawData.resize(nRows*nCols);
  for (unsigned i = 0; i < nRows; i++) {
    for (unsigned j = 0; j < nCols; j++) {
      unsigned dest_index = forFortran ? j*nRows+i : i*nCols+j;
      if (i >= old_rows || j >= old_cols) {
        rawData[dest_index] = 0;
      } else {
        unsigned src_index = forFortran ? j*old_rows+i : i*old_cols+j;
        rawData[dest_index] = old_data[src_index];
      } 
    }
  }
}

template< typename T >
SurfpackMatrix< T >& SurfpackMatrix< T >::operator=(const SurfpackMatrix<T>& other)
{
 ///\todo check for assignment to self
 forFortran = other.forFortran;
 nRows = other.nRows;
 nCols = other.nCols;
 rawData = other.rawData;
 return *this;
}

template< typename T >
SurfpackMatrix<T>& SurfpackMatrix< T >::operator=(const std::vector< std::vector<T> > other)
{
  if (other.empty()) {
    reshape(0,0);
  } else {
    reshape(other.size(),other[0].size());
  } 
  ///\todo Use assign instead of using this clunky double loop
  for (unsigned i = 0; i < nRows; i++) {
    for (unsigned j = 0; j < nCols; j++) {
      (*this)(i,j) = other[i][j];
    }
  }
  return *this;

}

template< typename T >
const std::string SurfpackMatrix< T >::asString() const
{
  std::ostringstream os;
  os.precision(3);
  unsigned index;
  // In C, members of the row are contiguous.  In Fortran, if the matrix has
  // n rows, the elements of any row are spaced n elements apart.
  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      os << std::setw(7) << (*this)(r,c) << " ";
    }
    os << "\n";
  }
  return os.str();
}

template< typename T >
const std::string SurfpackMatrix< T >::asArrayString() const
{
  std::ostringstream os;
  for (unsigned i = 0; i < rawData.size(); i++) {
    os << rawData[i] << " "; 
  }
  return os.str();
}

template< typename T >
unsigned SurfpackMatrix< T >::getNRows(char transpose) const
{
  assert(transpose == 'N' || transpose == 'T');
  return (transpose == 'N') ? nRows : nCols;
}

template< typename T >
unsigned SurfpackMatrix< T >::getNCols(char transpose) const
{
  assert(transpose == 'N' || transpose == 'T');
  return (transpose == 'N') ? nCols : nRows;
}

template< typename T >
const T& SurfpackMatrix<T>::operator()(unsigned r, unsigned c) const
{
  return const_cast< SurfpackMatrix<T>& >(*this)(r,c);
}

template< typename T >
T& SurfpackMatrix<T>::operator()(unsigned r, unsigned c)
{
  if (forFortran) {
    return rawData[c*nRows+r];
  }
  return rawData[r*nCols+c];
}

template< typename T >
T SurfpackMatrix<T>::sum() const
{
  return std::accumulate(rawData.begin(),rawData.end(),0.0);
}

#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
template< typename T >
template< class Archive >
void SurfpackMatrix<T>::serialize(Archive & archive, 
				  const unsigned int version)
{
  archive & forFortran;
  archive & nRows;
  archive & nCols;
  archive & rawData;
}

#endif

typedef SurfpackMatrix<double> MtxDbl;

#endif
