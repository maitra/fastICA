/**
 * \file order.h
 * Functions to permute indices to put arrays in ascending order, like 
 * <a href="http://www.r-project.org/">R</a>'s order() function.
 * For more information, see order.c.
 *
 * @author David Faden, dfaden@gmail.com
 * @author Karin S. Dorman kdorman@iastate.edu
 * @date August 21, 2005
 * @date Sun Sep 11 23:58:31 CDT 2011
 */

#ifndef __ORDER_H__
#define __ORDER_H__

#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "constants.h"
#include "array.h"

size_t* orderInt(const int* base, size_t numElements);
size_t* orderSize_t(const size_t* base, size_t numElements);
size_t* orderDouble(const double* base, size_t numElements);
size_t* orderString(char* const* base, size_t numElements);

/**
 * Comparison function for two objects.
 * Examples of this type of function include #intCompare(), #doubleCompare(),
 * and #stringCompare().
 * @param v1 void pointer to first object
 * @param v2 void pointer to second object
 * @return
 * 	-1 if v1 is "less than" v2,
 * 	0 if v1 "equals" v2, though see RETURN_RAND_CMP
 * 	1 if v1 is "greater than" v2.
 */
typedef int (*ComparisonFunc)(const void* v1, const void* v2);

/* some valid comparison functions */
int intCompare(const void *, const void *);
int doubleCompare(const void *, const void *);
int stringCompare(const void *, const void *);

/**
 * Comparison macro for numeric data types.
 * Returns -1 if a < b, 1 if a > b, and 0 otherwise from a function.
 * Useful to return from #ComparisonFunc and #CompareVectorElts functions.
 */
#define RETURN_CMP(a,b) if ((a) < (b)) { return -1; } \
	else if ((a) == (b)) { return 0; }            \
	else { return 1; }

/**
 * Random comparison macro for numeric data types.
 * Returns -1 if a < b, 1 if a > b, and -1 or 1 if a == b.
 * Useful to return from #ComparisonFunc and #CompareVectorElts functions.
 */
#define RETURN_RAND_CMP(a, b) if ((a) < (b)) { return -1; } \
	else if ((a) > (b)) { return 1;}                    \
	else if (rand()/RAND_MAX < 0.5) {return -1; }       \
	else { return 1; }


size_t* order(const void* base, size_t numElements, size_t size,
	      ComparisonFunc compare);


/**
 * Macro to define an #order() wrapper function.
 *
 * This macro takes a TYPE argument and defines a wrapper function named
 * order_TYPE.  The function takes two arguments, a const TYPE * pointer for an
 * array of TYPE elements, and a size_t, indicating the length of the array.
 * The function body defines a #ComparisonFunc-compliant function called compare_TYPE and
 * then calls the #order() function for arrays of type TYPE.  For example, the
 * call MAKE_ORDER_FUNC(int) would define a function order_int(const int *,
 * size_t), which is identical to #orderInt().

 *
 * @param TYPE data type for which you want an order function
 */
#define MAKE_ORDER_FUNC(TYPE) size_t* order_ ## TYPE (const TYPE* a, size_t len) \
{                                                                                \
	int compare_ ## TYPE (const void* v1, const void* v2) {                  \
		TYPE t1 = *((const TYPE*)v1);                                    \
		TYPE t2 = *((const TYPE*)v2);                                    \
		if (t1 < t2) { return -1; }                                      \
		else if (t1 == t2) { return 0; }                                 \
		else { return 1; }                                               \
	}                                                                        \
	return order(a, len, sizeof( TYPE ), compare_ ## TYPE );                 \
}


/* sorting that uses less memory, but requires C99 or greater standard */
#if (defined C99 || defined C11)

/**
 * Comparison function for elements within an array of simple data types.
 * Examples of this type of function include #compare_int_elts(),
 * #compare_double_elts(), and #compare_string_elts().
 *
 * @param v vector
 * @param left left index
 * @param right right index
 * @return -1 if v[index[left]] < v[index[right]],
 *   1 if v[index[left]] > v[index[right]], 0 otherwise
 */
typedef int (*CompareVectorElts)(const void *v, size_t *index, size_t left, size_t right, va_list);

/* some order functions */
void with_index_quicksort(void *, size_t *, CompareVectorElts, size_t, size_t, ...);
size_t *index_quicksort(void *, CompareVectorElts, size_t, size_t, ...);
size_t *order_int_simple(int *, size_t);
size_t *order_size_t_simple(size_t *, size_t);
size_t *order_double_simple(double *, size_t);

/* some valid comparison functions */
int compare_int_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_size_t_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_double_elts(const void *, size_t *, size_t, size_t, va_list);
int compare_string_elts(const void *, size_t *, size_t, size_t, va_list);

#endif

#endif
