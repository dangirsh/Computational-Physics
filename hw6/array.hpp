/* ---------------- array.hpp -----------------------

   template to make a 1D and/or 2D data type (array, vector or matrix)
   in C++,  with optional array bounds checking

   usually the type (T) is meant to be int, char, float, double
   other data types may or may not work

   define the symbol ARRAYT_BOUNDS_CHECK before including this
   file to enable bounds checking

  ----------------------------------------------

   functions:
    n1() = return 1st dimension
	n2() = return 2nd dimension
	ndim() = return number of dimensions
	n()    = return total size
	resize( n1 ) = change size to n1 (old data is lost)
	resize( n1, n2 ) = change size to n1 x n2 (old data is lost)

  if a1 and a2 are arrays (1, or 2 dimensions)
  of type T type T (usually float or double),
  then the following operations are allowed:

  array<T>(n1)    : construct a 1D array of size n1
  array<T>(n1,n2) : construct a 2D array of size (n1 x n2)

  a1(i)     : reference to element i of 1D array (vector) a1()
                    i ranges from 0 to n1-1
  a1(i,j)   : reference to element i,j of 1D array (matrix) a1
                    i ranges from 0 to n1-1 and j from 0 to n2-1

  a1 = a2   : a1 gets a copy of a2 (must be the same type)

  a1 += a2  : add a2 to a1 (element by element)

  NOTE:  Binary operations such as a1=a2+a3 should NOT be used because
      they are extremely inefficient (they create temporary arrays)

  -------------------------------------------------

  this code has been tested under MS-Visual C/C++ 6.0
  and Red Hat Linux 7.2 (gcc/g++ 2.96)

  this source code is formatted for a TAB size of 4 characters

   references:

	[1] Dov Bulka and David Mayhew, Efficient C++, Performance
        Performance Programming Techniques, Addison-Wesley 2000

	[2] D. M. Capper, Introducing C++ for Scientists, Engineers and
         Mathematicians, Springer-Verlag, 1994

	[3] James T. Smith, C++ Toolkit for Engineers and Scientists,
	     2nd edition, Springer 1999

	[4] B. Stroustrup, The C++ Programming Language 2nd edit.
		Addison Wesley 1991

	[5] D. Yang, C++ and Object Oriented Numeric Comp. for Sci. and Engin.,
	          (Springer), 2000, chapter 6

   started from matrix.hpp 20-jun-2001 E. Kirkland
   add test of nndim with bounds_check 1-oct-2002 ejk
   fix typo in destructor, and == 28-feb-2007 ejk
*/

#ifndef ARRAYT_HPP	// only include this file if its not already

#define ARRAYT_HPP	// remember that this has been included

// define the following symbol to enable bounds checking
//   can be defined here or in main calling program
//#define ARRAYT_BOUNDS_CHECK

#include <stdio.h>
#include <string.h>  // for memcpy()


//--- class definition -------------------------------------------

template < class T >
class array {
public:
	// constructor functions
	array( const int n1 );					// for 1D vector style
	array( const int n1, const int n2 );	// for 2D matrix style
	array( const array<T> &a );

	//  destructor function
	inline ~array() { delete [] p; nn=nndim=nn1=nn2=0; }

	// member operations
	inline array<T>& operator=( const array<T> &m );
	inline T& operator()( const int i1, const int i2);		// matrix
	inline T& operator()( const int i );	               	// vector
	inline array<T>& operator+=( const array<T> &m );

	// extra functions
	inline int n1() const { return nn1; }
	inline int n2() const { return nn2; }
	inline int ndim() const { return nndim; }
	inline int n() const { return nn; }
	void resize( const int n1, const int n2 );	 // matrix
	void resize( const int n ); 				 // vector

private:	// keep these read-only so they can't be accidentally changed

	T *p;					// pointer to storage area
	int nn;					// total number of elements
	int nndim;				// number of dimensions
	int nn1, nn2;			// size of each dimension
};

//--- constructor functions -------------------------------------------

template < class T >
array<T>::array( const int n1 )			// 1D vector
{
	if( n1 <= 0  ) {
		printf("array initialized with size = %d, NOT ALLOWED\n", n1 );
		exit( EXIT_FAILURE );
	}
	p = new T [ n1 ];
	if( NULL == p ) {
		printf( "Cannot allocate memory for array size = %d\n", n1 );
		exit( EXIT_FAILURE );
	}
	nn = nn1 = n1;
	nn2 = 0;
	nndim = 1;
}

template < class T >
array<T>::array( const int n1, const int n2 )		// 2D matrix
{
	if( ( n1 <= 0 ) || ( n2 <= 0 ) ) {
		printf("array initialized with size = %d x %d, NOT ALLOWED\n",
			n1, n2 );
		exit( EXIT_FAILURE );
	}
	p = new T [ n1*n2 ];
	if( NULL == p ) {
		printf( "Cannot allocate memory for array size = %d x %d\n",
			n1, n2 );
		exit( EXIT_FAILURE );
	}
	nn1 = n1;
	nn2 = n2;
	nn = n1*n2;
	nndim = 2;
}

template < class T >				// required for misc. operations
array<T>::array( const array<T> &a )
{
	p = new T [ a.nn ];
	if( NULL == p ) {
		printf( "Cannot allocate memory for array size = %d\n", a.nn );
		exit( EXIT_FAILURE );
	}
	memcpy( p, a.p, a.nn*sizeof(T) );
	nn1 = a.nn1;
	nn2 = a.nn2;
	nndim = a.nndim;
	nn = a.nn;
}

// -------  member function resize() -------------------------

template < class T >
void array<T>::resize( const int n  )		// 1D resize
{
	if( n <= 0  ) {
		printf("array resize() with size = %d, NOT ALLOWED\n",
			n );
		exit( EXIT_FAILURE );
	}
	delete [] p;
	p = new T [ n ];
	if( NULL == p ) {
		printf( "Cannot allocate memory for array resize = %d\n",
			n );
		exit( EXIT_FAILURE );
	}
	nn1 = n;
	nn2 = 0;
	nndim = 1;
	nn = n;
}

template < class T >
void array<T>::resize( const int n1, const int n2 )		// 2D resize
{
	if( ( n1 <= 0 ) || ( n2 <= 0 ) ) {
		printf("array resize() with size = %d x %d, NOT ALLOWED\n",
			n1, n2 );
		exit( EXIT_FAILURE );
	}
	delete [] p;
	p = new T [ n1 * n2 ];
	if( NULL == p ) {
		printf( "Cannot allocate memory for array resize = %d x %d\n",
			n1, n2 );
		exit( EXIT_FAILURE );
	}
	nn1 = n1;
	nn2 = n2;
	nndim = 2;
	nn = n1 * n2;
}

// ------- operator functions ----------------------------------


// -------  member function operator =
template < class T >
array<T>& array<T>::operator=( const array<T> &m )
{
	if( (nn != m.nn) || (nndim != m.nndim)  ){
		printf( "array = operator invoked with unequal sizes\n"
			"  m1 size = %d, dim= %d\n"
			"  m2 size = %d, dim= %d\n", nn, nndim, m.nn, m.nndim );
		exit( EXIT_FAILURE );
	} else {
		memcpy( p, m.p, nn*sizeof(T) );	// fastest way to do this
		return *this;
	}
}

// ------- index operators ----------------------------------
//
//  remember: [] only allows one argument so can't be used for > 1D
//
// ------- member function operator () = 2D index
template < class T >
inline T& array<T>::operator()( const int i1, const int i2 )
{
#ifdef ARRAYT_BOUNDS_CHECK
	if( (i1<0) || (i1>=nn1) ||
		(i2<0) || (i2>=nn2) || (nndim != 2 ) ){
		printf( "out of bounds index in array\n"
			"  size = %d x %d, ndim= %d\n  access = (%d , %d)\n",
			nn1, nn2, nndim, i1, i2 );
		exit( EXIT_FAILURE );
	}
#endif

	// both should work but one may be faster
	//   for different operations
	return *(p + i2 + i1*nn2);
//	return *(p + i1 + i2*nn1);
}

// ------- member function operator () = 1D index
template < class T >
inline T& array<T>::operator()( const int i )
{
#ifdef ARRAYT_BOUNDS_CHECK
	if( (i<0) || (i>=nn) || (nndim != 1 )){
		printf( "out of bounds index in array\n"
			"  size = %d, ndim= %d\n  access = %d\n", nn, nndim, i );
		exit( EXIT_FAILURE );
	}
#endif

	return *(p + i);
}


// ------- member function operator +=
//  should work for any number of dimensions
template < class T >
inline array<T>& array<T>::operator+=( const array<T>& m  )
{
	if( (m.nn != nn)  ){
		printf( "array += operator invoked with unequal sizes:\n"
			"   %d and %d\n", nn,  m.n() );
		exit( EXIT_FAILURE );
	} else {
		register int i;
		for( i=0; i<nn; i++) p[i] += m.p[i];
		return *this;
	}
}

#endif  // ARRAYT_HPP

