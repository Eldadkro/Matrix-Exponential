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
        """
        returns the exp(X)
        """
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
            

class LUsolver:
    """
    class using a builder Design pattern that decposition a matrix to LU decomp and solves the linear equation using the decomp
    """
    def __init__(self,A:np.array, b:np.array) -> None:
        self.A = A
        self.b = b
        self.L = None
        self.U = None
    
    def _decomp_rec(A):
        if A.shape[0] == 1:
            return 1,A
        schur_comp = A[1:,1:] - A[1:,0]@A[0,1:]/A[0:0]
        L_shur,U_shur = LUsolver._decomp_rec(schur_comp)
        L = np.eye(A.shape[0],A.shape[0])
        L[0,1:] = A[1:,0] / A[0,0]
        L[1:,1:] = L_shur
        U = A
        U[1:,0] = np.zeros((A.shjape[0] - 1,1))
        U[1:,1:] = U_shur
        return L,U

    def decomp(self):
        """
        decomp a matrix to L,U 
        """
        self.L,self.U = LUsolver._decomp_rec(self.A)
        return self

   

    def solve(self):
        y = np.empty((self.A.shape[0],1))
        y[0] = self.b[0]
        for i in range(1,self.A.shape[0]):
            y[i] = self.b[i] - self.L[0:i]@y[0:i] #need to check 
        
        x = np.empty((self.A.shape[0],1))
        x[-1] = y[-1]/self.U[-1,-1]
        for i in range(self.A - 2,-1,-1):
            x[i] = (y[i] - self.U[i,i+1:]@x[i+1:])/self.U[i,i]
        return x
        
    




if __name__ == "main":
    pass