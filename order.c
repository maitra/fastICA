/**
 * \file order.c
 *
 * Routines that return a permutation of indices to put an array in sorted
 * order, like <a href="http://www.r-project.org">R</a>'s order function.
 *
 * For arrays of simple data types, see #simple_int_order() for integers and
 * #index_quicksort() for any array of simple data types, provided a
 * #CompareVectorElts function is provided.  You can write your own
 * #CompareVectorElts function, but see #compare_int_elts,
 * #compare_double_elts, and #compare_string_elts for comparing
 * elements of an integer, double, or const string array.
 *
 * See #order() for a generic function that can order indices for any type of
 * object for which the user can provide a #ComparisonFunc.  You can write
 * your own #ComparisonFunc function, but see #intCompare, #doubleCompare, and
 * #stringCompare for comparing int, double, or strings.  The functions
 * #orderInt(), #orderDouble(), and #orderString() are wrapper functions for
 * #order() using these #ComparisonFunc functions, respectively.
 *
 * @author David Faden, dfaden@gmail.com
 * @author Karin S. Dorman, kdorman@iastate.edu
 * @date August 21, 2005
 * @date Sun Sep 11 14:43:36 CDT 2011
 */

#include "order.h"

/**
 * Pair struct for agnostic sorting of objects.
 */
typedef struct {
	size_t index;			/*!< index of object in array */
	const void* datum;		/*!< object in array */
	ComparisonFunc compare;	/*!< pointer to comparison function */
} Pair;

/**
 * Compare two #Pair objects using user-provided comparison function.
 *
 * Generic comparison function for two objects.  The arguments are of type
 * #Pair, which provides a comparison function for doing the actual comparison.
 * If the comparison functions of the two objects match, then use it to compare
 * the two objects.
 *
 * @param v1 pointer to first object
 * @param v2 pointer to second object
 * @return -1 if v1 < v2, 0 if v1 == v2, and 1 otherwise
 */
int comparePairs(const void* v1, const void* v2)
{
	const Pair* p1 = v1;
	const Pair* p2 = v2;
	assert(p1->compare == p2->compare);
	return p1->compare(p1->datum, p2->datum);
} /* comparePairs */

/**
 * Function returning index order of an arbitrary array of objects.
 *
 * This function uses the C library function qsort to find index order
 * that would sort arbitrary objects in increasing order according to a
 * user-defined comparison.  This function is very flexible with regards
 * to the types of objects it can sort, but requires considerable memory.
 * In particular it needs to allocate two pointers and a size_t index
 * for each element in the array.  See #simple_int_order() for 
 * index ordering of simple integer arrays or #index_quicksort() for
 * index ordering of arrays of an arbitrary, simple data type.
 *
 * The caller is responsible for freeing the returned object.
 *
 * @param base array of objects
 * @param numElements length of array
 * @param size size in bytes of elements in array
 * @param compare a function that is capable of comparing objects in the array
 * @return indices of array in order to produce objects in increasing order
 */
size_t* order(const void* base, size_t numElements, size_t size,
				ComparisonFunc compare)
{
	Pair *pairs = malloc(sizeof(Pair) * numElements);
	size_t *indices = malloc(sizeof(size_t) * numElements);
	size_t i;
	const char *p = base;

	for (i = 0; i < numElements; ++i, p += size) {
		pairs[i].index = i;
		pairs[i].datum = p;
		pairs[i].compare = compare;
	}

	qsort(pairs, numElements, sizeof(Pair), comparePairs);

	for (i = 0; i < numElements; ++i)
		indices[i] = pairs[i].index;

	free(pairs);

	return indices;
} /* order */


/**
 * A #ComparisonFunc-compliant function for two integers.
 * @param v1 pointer to first integer
 * @param v2 pointer to second integer
 * @return comparison value, see #RETURN_CMP.
 */
int intCompare(const void* v1, const void* v2)
{
	int i1 = *((const int*)v1);
	int i2 = *((const int*)v2);
	RETURN_CMP(i1, i2);
}

/**
 * A #ComparisonFunc-compliant function for two size_t positive integers.
 * @param v1 pointer to first size_t
 * @param v2 pointer to second size_t
 * @return comparison value, see #RETURN_CMP.
 */
int size_tCompare(const void* v1, const void* v2)
{
	size_t i1 = *((const size_t*)v1);
	size_t i2 = *((const size_t*)v2);
	RETURN_CMP(i1, i2);
}

/**
 * A #ComparisonFunc-compliant function for two doubles.
 * @param v1 pointer to first double
 * @param v2 pointer to second double
 * @return comparison value, see #RETURN_CMP.
 */
int doubleCompare(const void* v1, const void* v2)
{
	double d1 = *((const double*)v1);
	double d2 = *((const double*)v2);
	RETURN_CMP(d1, d2);
}

/**
 * A #ComparisonFunc-compliant function for two strings.
 * @param v1 pointer to first string
 * @param v2 pointer to second string
 * @return comparison value, see #RETURN_CMP.
 */
int stringCompare(const void* v1, const void* v2)
{
	return strcmp(*(const char**)v1, *(const char**)v2);
}

/**
 * Wrapper for #order() to produce index order of integer array.
 * @param base array of integers
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* orderInt(const int* base, size_t numElements)
{
	return order(base, numElements, sizeof(int), intCompare);
}

/**
 * Wrapper for #order() to produce index order of size_t array.
 * @param base array of size_t
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* orderSize_t(const size_t* base, size_t numElements)
{
	return order(base, numElements, sizeof(size_t), size_tCompare);
}

/**
 * Wrapper for #order() to produce index order of double array.
 * @param base array of doubles
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* orderDouble(const double* base, size_t numElements)
{
	return order(base, numElements, sizeof(double), doubleCompare); 
}

/**
 * Wrapper for #order() to produce index order of string array.
 * @param base array of strings
 * @param numElements number of elements in array
 * @return sorted index array for array in increasing order
 */
size_t* orderString(char* const* base, size_t numElements)
{
	return order(base, numElements, sizeof(const char*), stringCompare);
}

/* sorting that uses less memory, but requires C99 or greater standard */
#if (defined C99 || defined C11)
size_t index_partition(void *, size_t *, CompareVectorElts, size_t, size_t, size_t, va_list);
void vindex_quicksort(void *, size_t *, CompareVectorElts, size_t, size_t, va_list);

/**
 * A #index_quicksort() wrapper for index sorting of integer arrays.
 *
 * This function allocates an index array, and then uses the quick sort 
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec integer vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_int_simple(int *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_int_elts, 0, n-1);
} /* order_int_simple */

/**
 * A #index_quicksort() wrapper for index sorting of size_t positive integer
 * arrays.
 *
 * This function allocates an index array, and then uses the quick sort 
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort directly with the appropriate comparisonFunc.
 *
 * @param vec size_t vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_size_t_simple(size_t *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_size_t_elts, 0, n-1);
} /* order_size_t_simple */

/**
 * A #index_quicksort() wrapper for index sorting of double arrays.
 *
 * This function allocates an index array, and then uses the quick sort 
 * algorithm to sort the index in increasing order.
 *
 * It is possible to write other functions like this for other simple data types
 * or simply call index_quicksort() directly with the appropriate comparisonFunc.
 *
 * @param vec double vector to index sort
 * @param n length of vector
 * @return indices sorted in order of increasing array value
 */
size_t *
order_double_simple(double *vec, size_t n)
{
	return index_quicksort((void *)vec, compare_double_elts, 0, n-1);
} /* order_double_simple */

/**
 * Quick sort by index for arrays of simple data types.
 * The permuted index it returns would put the array in first argument vec into
 * ascending order.  Order is determined by the user-provided last argument
 * compare.
 *
 * @param vec array of simple data types to order
 * @param left left most index (inclusive)
 * @param right right most index (inclusive)
 * @param compare #ComparisonFunc for comparing elements in array vec
 * @return ordered index array
 */
size_t *
index_quicksort(void *vec, CompareVectorElts compare,
	size_t left, size_t right, ...)
{
	size_t i, last;
	va_list vl;
	size_t *index;

	/* we will actuall only sort indices */
	last = right - left + 1;
	MAKE_VECTOR(index, last);

	if (!index)
		return NULL;

	for (i = 0; i < last; i++)
		index[i] = i;

	va_start(vl, right);

	vindex_quicksort(vec, index, compare, left, right, vl);

	va_end(vl);

	return index;
} /* index_quicksort */

void
with_index_quicksort(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, ...)
{
	va_list vl;
	va_start(vl, right);
	vindex_quicksort(vec, index, compare, left, right, vl);
	va_end(vl);
} /* with_index_quicksort */

/**
 * va_list version of index_quicksort
 */
void
vindex_quicksort(void *vec, size_t *index, CompareVectorElts compare,
	size_t left, size_t right, va_list vl)
{
	size_t pivot, middle;
	int compare1, compare2, compare3;
	va_list vl2;

	if (left < right) {
		/* choose median of first, middle, last elements (Wikipedia) */
		middle = (right - left)/2 + left;
		va_copy(vl2, vl);
		compare1 = compare(vec, index, left, middle, vl2);	/* !!! */
		va_end(vl2);
		va_copy(vl2, vl);
		compare2 = compare(vec, index, middle, right, vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		compare3 = compare(vec, index, left, right, vl2);
		va_end(vl2);
		if ((compare1 <= 0 && compare2 <= 0)		/* l < m < r */
			|| (compare1 >= 0 && compare2 >=0))	/* r < m < l */
			pivot = middle;
		else if((compare3 <= 0 && compare2 >= 0)	/* l < r < m */
			|| (compare2 <= 0 && compare3 >= 0))	/* m < r < l */
			pivot = right;
		else
			pivot = left;
		va_copy(vl2, vl);
		pivot = index_partition(vec, index, compare, left, right, pivot,
			vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		if (pivot) vindex_quicksort(vec, index, compare, left, pivot - 1, vl2);
		va_end(vl2);
		va_copy(vl2, vl);
		vindex_quicksort(vec, index, compare, pivot + 1, right, vl2);
		va_end(vl2);
	}
} /* vindex_quicksort */

/**
 * Quick sort partition algorithm to work with #index_quicksort().  See 
 * <a href="http://en.wikipedia.org/wiki/Quicksort">Wikipedia</a> for more
 * information.
 *
 * @param vec integer vector whose ordering is sought
 * @param index order vector
 * @param left left index (inclusive)
 * @param right right index (inclusive)
 * @param pivot input pivot index 
 * @param compare #ComparisonFunc
 * @return new pivot
 */
size_t
index_partition(void *vec, size_t *index, CompareVectorElts compare, 
	size_t left, size_t right, size_t pivot, va_list vl)
{
//const double **dvec = (const double **)vec;
	size_t i, new_pivot, swapper;
	va_list vl2;

	/* move pivot to end of array */
	swapper = index[right];
	index[right] = index[pivot];
	index[pivot] = swapper;

	/* sort this partition */
	new_pivot = left;
	for (i = left; i < right; i++) {
		va_copy(vl2, vl);
		if (compare(vec, index, i, right, vl2) < 0) {
			swapper = index[new_pivot];
			index[new_pivot] = index[i];
			index[i] = swapper;
			new_pivot++;
		}
		va_end(vl2);
	}

	/* move pivot to appropriate position */
	swapper = index[new_pivot];
	index[new_pivot] = index[right];
	index[right] = swapper;

	return new_pivot;
} /* index_partition */

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an integer array.
 * @param v vector of integers
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_int_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	const int *iv = (int *)v;
	RETURN_CMP(iv[index[left]], iv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in an size_t array.
 * @param v vector of size_t
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_size_t_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	const size_t *iv = (size_t *)v;
	RETURN_CMP(iv[index[left]], iv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in a double array.
 * @param v vector of doubles 
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_double_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	const double *dv = (double *)v;
	RETURN_CMP(dv[index[left]], dv[index[right]]);
}

/**
 * A #CompareVectorElts-compliant function for comparing two elements in a string array.
 * @param v vector of const strings
 * @param index index array
 * @param left index of first element
 * @param right index of second element
 * @return comparison value, see #RETURN_CMP.
 */
int compare_string_elts(const void *v, size_t *index, size_t left, size_t right,
	va_list vl)
{
	const char **sv = (const char **)v;
	return strcmp(sv[index[left]], sv[index[right]]);
}
#endif


