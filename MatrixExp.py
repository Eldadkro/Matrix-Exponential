import numpy as np
import scipy as sp

b13 = []
teta13 = []

class Pade:
    """
    Pade represents a Pade approximation for a matrix exponension function as writen in [1].

    refernces:
    [1] THE SCALING AND SQUARING METHOD FOR THE MATRIX EXPONENTIAL REVISITED, Nicholas J. Higham,pages 1183-1184
    """
    def __init__(self, coafs) -> None:
        self.p_coafs = coafs
        self.degree = len(coafs)

    def evaluate(self,X):
        U,V = self._evalUV(X)
        p_m = U+V
        q_m = (-1)**(self.degree%2)*U + (-1)**(self.degree%2 + 1) * V
        #solver
        exp_X = np.linalg.solve(q_m,p_m)
        return exp_X

    def _evalUV(self, X:np.array, coafs: np.array):
        """
        returns the matrices: U,V for the polynomialy Pm, Qm
        """
        res = np.array((X.shape[0],X.shape[1],self.degree//2 + 1))
        res[:,:,0] = np.eye(X.shape[0],X.shape[1])
        X_square = X@X
        for i in range(1,self.degree//2 + 1):
            res[:,:,i] = res[:,:,-1]@X_square
        # for U
        U = np.zeros((X.shape[0],X.shape[1]))
        for i in range(self.degree//2 + 1):
            U += coafs[i] * res[:,:i]
        #for V
        V = np.zeros((X.shape[0],X.shape[1]))
        for i in range(self.degree//2 + self.degree%2):
            V += coafs[i] * res[:,:,i]
        return U,V
            







if __name__ == "main":
    pass