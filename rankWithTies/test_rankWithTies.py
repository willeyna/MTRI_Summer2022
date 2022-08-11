#! python3
import unittest
import numpy as np
from rankWithTies import rankWithTies

class TestCases(unittest.TestCase):
    '''
     test_rankWithTies - Test suite

     PURPOSE:
       This code provides the test suite for rankWithTies through unittest.
       New bug reports for the code should not be closed with ensuring that a
       test in the suite both fails before the fix, and passes after it.
    '''

    def test_syntaxSupport(self):
        """
        Test that all of the described syntax is supported
        """
        AList = [np.array([np.pi, 1, 1, 2, 3, 3])]
        methodList = ['min', 'max', 'random', 'stable', \
            lambda rMax, nMatch : np.arange((rMax-(nMatch-1)), rMax+1)]
        sortOptsList = [[0], [0, 'mergesort'], []]

        for A in AList:
            for method in methodList:
                for opts in sortOptsList:
                    args = [A, method, *opts]

                    try:
                        R, B, I = rankWithTies(*args)
                    except:
                        self.fail(f"Unsupported argument set: {args}")

    def test_typeSupport(self):
        '''
        Tests input/output type support
        '''
        A = 'fooBar...'
        #different datatype input arrays for ranking
        AList = [np.array([float(ord(char)) for char in A]),
                 np.array([np.float32(ord(char)) for char in A]),
                 np.array([np.uint8(ord(char)) for char in A]),
                 np.array(list(A))]
        methodList = ['min', 'max', 'random', 'stable', \
            lambda rMax, nMatch : np.arange((rMax-(nMatch-1)), rMax+1)]
        sortOptsList = [[0], [0, 'mergesort'], []]

        for A in AList:
            for method in methodList:
                for opts in sortOptsList:
                    args = [A, method, *opts]

                    R, B, I = rankWithTies(*args)

                    self.assertEqual(R.dtype, np.int64)
                    self.assertEqual(B.dtype, args[0].dtype)
                    self.assertEqual(I.dtype, np.int64)

    def test_dimensionSupport(self):
        '''
        Tests input dimensionality support
        '''

        A = np.array(list('fooBar...'))
        nA = len(A)

        for dimInd in range(0,4):
            #reshape A to be along different dimensions
            aSize = np.ones(max(2, dimInd+1), dtype = 'int64')
            aSize[dimInd] = nA
            A = A.reshape(aSize)

            #Confirm the dimensions match
            R, B, I = rankWithTies( A, 'min')

            np.testing.assert_array_equal(R.shape, aSize)
            np.testing.assert_array_equal(B.shape, aSize)
            np.testing.assert_array_equal(I.shape, aSize)

    def test_permutationProperties(self):
        '''
        Tests whether the returned perumation arrays correctly sort and unsort A
        '''

        A = np.array(list('fooBar...'))
        methodList = ['min', 'max', 'random', 'stable', \
            lambda rMax, nMatch : np.arange((rMax-(nMatch-1)), rMax+1)]

        for method in methodList:
            R, B, I = rankWithTies(A, method)
            np.testing.assert_array_equal(B[R], A)
            np.testing.assert_array_equal(A[I], B)

    def test_methodFunctionalEquivalents(self):
        '''
        Tests built-in methods' functional equivalents
        Must carefully keep track of the numpy random state for method=random
        '''

        methodList = ['min', 'max', 'random', 'stable']
        equivalent = [lambda rMax, nMatch : rMax - (nMatch-1),
                      lambda rMax, nMatch : rMax,
                      lambda rMax, nMatch : np.random.permutation(nMatch) \
                       + (rMax - (nMatch - 1)),
                      lambda rMax, nMatch : \
                        np.arange((rMax-(nMatch-1)), rMax+1)]

        A = np.array(list('fooBar...'))

        for i in range(len(methodList)):

            np.random.seed(25)
            R1, B1, I1 = rankWithTies(A, methodList[i])

            np.random.seed(25)
            R2, B2, I2 = rankWithTies(A, equivalent[i])

            np.testing.assert_array_equal(R1, R2)
            np.testing.assert_array_equal(B1, B2)
            np.testing.assert_array_equal(I1, I2)


    #Skipped testing optional inputs for now as the cases are very different
    #   in numpy since asc/desc is not an arg and order matters in numpy's
    #   sort function inputs


    def test_vectorization(self):
        '''
        Test input array vectorization support
        '''

        methodList = ['min', 'max', 'random', 'stable', \
            lambda rMax, nMatch : np.arange((rMax-(nMatch-1)), rMax+1)]

        #make sure the vectorization holds for every method
        for method in methodList:
            aSize = [5,6]
            A = np.random.random(aSize)
            RTrue = np.zeros_like(A, dtype = int)
            BTrue = np.zeros_like(A)
            ITrue = np.zeros_like(A, dtype = int)

            #go down in column order and manually do 1D rankings
            for i in np.arange(0, aSize[1])[::-1]:
                RTrue[:,i], BTrue[:,i], ITrue[:,i]= rankWithTies(A[:,i], method)

            R, B, I = rankWithTies(A, method, 0)
            np.testing.assert_array_equal(R, RTrue) #explicit along dim 0
            np.testing.assert_array_equal(B, BTrue)
            np.testing.assert_array_equal(I, ITrue)

            for i in np.arange(0, aSize[0])[::-1]:
                RTrue[i,:], BTrue[i,:], ITrue[i,:]= rankWithTies(A[i,:], method)

            R, B, I = rankWithTies(A, method)
            np.testing.assert_array_equal(R, RTrue) #implicitly along dim 1
            np.testing.assert_array_equal(B, BTrue)
            np.testing.assert_array_equal(I, ITrue)

            R, B, I = rankWithTies(A, method)
            np.testing.assert_array_equal(R, RTrue, 1) #explicitly along dim 1
            np.testing.assert_array_equal(B, BTrue, 1)
            np.testing.assert_array_equal(I, ITrue, 1)


if __name__ == '__main__':
    unittest.main()
