namespace Utils {

/**
* \brief Returns the unravelled index of the provided flat index.
*        Therefore is the inversion of flattening an ndims dimensional index.
* @param len_dims an int array of length ndims containing the lengths of the dimensions. (Input)
* @param ndims int denoting the number of dimensions. (Input)
* @flattened_index an int denoting the flat index. (Input)
* @unravelled_index_out an int array with length ndims where the unflat indices are written to. (Output)
*/
inline void unravel_index(const int* const len_dims, const int ndims, const int flattened_index, int* unravelled_index_out){
	//idea taken from http://codinghighway.com/2014/02/22/c-multi-dimensional-arrays-part-2-flattened-to-unflattened-index/
    std::vector<int> mul(ndims);
	mul[ndims-1]=1;
	for (int j = ndims-2; j >= 0; j--)
		mul[j] = mul[j+1]*len_dims[j+1];
	for (int j = 0; j < ndims; j++)
		unravelled_index_out[j]=(flattened_index/mul[j])%len_dims[j];
}
} // Namespace Utils
