import numpy as np

def rankWithTies(A, method, *args):
    '''
    rankWithTies - Ranking with tie-break control
     SYNTAX:
       R,B,I = rankWithTies(A,method,<sortOpts>)

    PURPOSE:
    This function performs ranking with additional control of tie-breaks.
    Ranking is fundamentally a sorting operator where one reports the relative
    position of each element within the sorted list.  When ties are present,
    these relative positions become ambiguous, and this function controls how
    those positions get reported.

    The outputs R and I are indices that describe permutation operators with the
    properties A = B[R] and B = A[I] for vector-valued inputs.  R and I will not
    represent inverse permutation operators when ties are present if one does
    not ensure unique[R] = unique[I].  For example, the 'min' and 'max' methods
    below will lead to R and I representing permutations that are not inverses
    of one another when ties are present.

     INPUT:
       A           - Array to be sorted

       method      - Method for breaking ties.  Common defaults are provided and
                     one can also define an arbitrary function that acts on the
                     stable ranks
           'min'       - Smallest index
           'max'       - Largest index
           'random'    - Random permutation index
           'stable'    - Order of occurance in A

                     A custom function defined on the stable ranks given by
                       lambda rMax, nMatch :
                     This function must have to prototype
                       ROut(ind) = methodFun(rMax,nMatch)
                     where rMax is the maximum stable rank and nMatch is the
                     number of matches within the set B[rMax-(nMatch-1):rMax].
                     For the aforementioned defaults this function is:
           'min'       - lambda rMax, nMatch : rMax-(nMatch-1)
           'max'       - lambda rMax, nMatch : rMax
           'random'    - lambda rMax, nMatch :
                        np.random.permutation(nMatch) + (rMax - (nMatch - 1))
           'stable'    - lambda rMax, nMatch:np.arange((rMax-(nMatch-1)),rMax+1)


       sortOpts    - Options to numpy built-in sort call

     OUTPUT:
       R       - Ranks.  Reverse sort index s.t. A = B[R] along dim
       B       - Sort result
       I       - Sort index with ties broken via method s.t. B = A[I] along dim

     ASSUMPTIONS:
       All input variables are of the correct type, valid (if applicable),
       and given in the correct order.

     NOTES:
       Vectorization doesn't explicitly parallelize replaceDuplicateRanks
    '''
    #stores calls to built in methods
    DEFAULT_METHODS = {'min': lambda rMax, nMatch : rMax - (nMatch-1),
                       'max': lambda rMax, nMatch : rMax,
                       'random': lambda rMax, nMatch : \
                        np.random.permutation(nMatch) + (rMax - (nMatch - 1))
                      }

    #if no dimension specified, assume ranking along last dimension of A
    if not len(args):
        args += (-1,)

    dim = args[0]

    I = np.argsort(A, *args)
    #avoids sorting twice
    #indexes A along 1D slices in dim using indices in respective I slices
    B = np.take_along_axis(A, I, dim)

    #initialization of inverse index array
    R = np.zeros_like(I)

    N = A.shape[dim]

    #reshapes the arange(N) array to be along the axis to insert in R
    indexer = tuple([None if i!=dim else Ellipsis for i in range(A.ndim)])

    np.put_along_axis(R, I, np.arange(N, dtype = float)[indexer], axis=dim)

    if method == 'stable':
        return R, B, I

    elif method in DEFAULT_METHODS.keys():
        fun = DEFAULT_METHODS[method]

    elif callable(method):
        fun = method

    else:
        raise ValueError(f"Invalid method type: {method}")

    #move axis to loop through to end for ease of iteration/ dim generalization
    B_shift = np.moveaxis(B, dim, -1)
    I_shift = np.moveaxis(I, dim, -1)
    R_shift = np.moveaxis(R, dim, -1)

    for index in np.ndindex(B_shift.shape[:-1]):
        #perform rank tie-breaking along 1D slices
        R_shift[index] = replaceDuplicateRanks(B_shift[index], \
                                        I_shift[index], R_shift[index], fun)

    return R, B, I

def replaceDuplicateRanks(B, I, R, fun):

    #arrays with the indices the unique values start/stop at and their distances
    uniq, uniq_ind, uniq_count = \
                        np.unique(B, return_index = True, return_counts = True)

    #makes sure the to check for the last value being a duplicate
    uniq_ind = np.append(uniq_ind, len(B))

    for i in range(len(uniq_ind)-1):
        if uniq_count[i] != 1:
            #replaces rank ties according to fun
            R[I[uniq_ind[i]:uniq_ind[i+1]]]= fun(uniq_ind[i+1]-1, uniq_count[i])

    return R
