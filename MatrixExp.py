#%%
import numpy as np
import scipy as sp
import math
# import warnings
# warnings.filterwarnings("ignore")
b13 = np.array([64764752532480000, 32382376266240000, 7771770303897600,
                1187353796428800, 129060195264000, 10559470521600,
                670442572800, 33522128640, 1323241920,
                40840800, 960960, 16380, 182, 1])

teta13 = np.array([1.495585217958292e-2,
                2.539398330063230e-1,
                9.504178996162932e-1,
                2.097847961257068e0,
                5.371920351148152e0])

class KazmaarchSolver:

    def __init__(self,A,b,sweeps=100) -> None:
        self.A = A
        self.b = b
        self.sweeps= sweeps
    
    def solve(self,X):
        n = self.A.shape[0]
        for sweep in range(self.sweeps):
            for i in range(n):
                X = X + ((self.b[i,0] - (self.A[i,:].reshape((1,-1))).dot(X)[0,0])/(self.A[i,:]@np.transpose(self.A[i,:])))*(np.transpose(self.A[i,:]).reshape((-1,1)))
        return X    

class LUsolver:
    """
    class using a builder Design pattern that decposition a matrix to LU decomp and solves the linear equation using the decomp

    reference:
    [1] Introduction to algorithms, Thired edition, Thomas H.Cormen et al, pages 813-827.
    """
    def __init__(self,A:np.array, b:np.array) -> None:
        self.A = A
        self.b = b
        self.L = None
        self.U = None
    
    def _decomp_rec(A):
        if A.shape[0] == 1:
            return 1,A
        schur_comp = A[1:,1:] - A[1:,0]@A[0,1:]/A[0,0]
        L_shur,U_shur = LUsolver._decomp_rec(schur_comp)
        L = np.eye(A.shape[0],A.shape[0])
        L[0,1:] = A[1:,0] / A[0,0]
        L[1:,1:] = L_shur
        U = A
        try:
            U[1:,0] = np.zeros((A.shape[0] - 1,1)).reshape(U[1:,0].shape)
        except:
            print(U.shape)
            print(U[1:,0].shape)
            SystemExit()
        U[1:,1:] = U_shur
        return L,U

    def _decomp(A):
        n = A.shape[0]
        U = np.zeros((n,n))
        L = np.eye(n,n)

        for k in range(n-1):
            U[k,k] = A[k,k]
            L[k+1:n,k] = A[k+1:n,k]/A[k,k]
            U[k,k+1:n] = A[k,k+1:n]
            A[k+1:,k+1:]  =  A[k+1:,k+1:] - A[k+1:n,k]@A[k,k+1:n] / A[k,k]
        return L,U


    def decomp(self):
        """
        decomp a matrix to L,U 
        """
        self.L,self.U = LUsolver._decomp(self.A)
        return self

   

    def solve(self):
        y = np.empty((self.A.shape[0],1))
        y[0] = self.b[0]
        for i in range(1,self.A.shape[0]):
            y[i] = self.b[i] - self.L[i,0:i].dot(y[0:i]) #need to check 
        
        x = np.empty((self.A.shape[0],1))
        x[-1] = y[-1]/self.U[-1,-1]
        for i in range(self.A.shape[0] - 2,-1,-1):
            x[i] = (y[i] - self.U[i,i+1:]@x[i+1:])/self.U[i,i]
        return x

class Solver:
    
    def __init__(self,A:np.array,res:np.array) -> None:
        self.A = A
        self.res = res
    
    def solve(self):
        """
        solve using LU decompositsion for each colum in res and X, where AX = res
        A,X,res have the same size of nXn
        """
        n = self.A.shape[0]
        X = np.empty((n,n))
        #for each colum of res
        for i in range(self.res.shape[1]):
            X[:,i] =  LUsolver(self.A,self.res[:,i]).decomp().solve().reshape(X[:,i].shape)
        return X

class Ksolver(Solver):
    def __init__(self, A: np.array, res: np.array) -> None:
        super().__init__(A, res)
    
    def solve(self):
        n = self.A.shape[0]
        X = np.empty((n,n))
        for i in range(self.res.shape[1]):
            K = KazmaarchSolver(self.A,self.res[:,i].reshape((-1,1)))
            sol = K.solve(X[:,i].reshape((-1,1)))
            X[:,i] =  sol.reshape(X[:,i].shape)
        return X

class Pade:
    """
    Pade represents a Pade approximation for a matrix exponension function as writen in [1].

    
    """
    def __init__(self, coafs:np.array) -> None:
        self.p_coafs = coafs
        self.degree = len(coafs) - 1
    
    def _evalUV(self, X:np.array, coafs: np.array):
        """
        returns the matrices: U,V for the polynomialy Pm, Qm
        """
        res = np.empty((X.shape[0],X.shape[1],self.degree + 1))
        res[:,:,0] = np.eye(X.shape[0],X.shape[1])
        X_square = X@X
        
        for i in range(1,self.degree+1):
            res[:,:,i] = res[:,:,i-1]@X

        #for even == U
        U = np.zeros((X.shape[0],X.shape[1]))
        for i in range(0,self.degree + 1,2):
            U += coafs[i] * res[:,:,i]
        V = np.zeros((X.shape[0],X.shape[1]))
        for i in range(1,self.degree + 1,2):
            V += coafs[i] * res[:,:,i]

        return U,V

    def evaluate(self,X):
        """
        returns the exp(X)
        """
        U,V = self._evalUV(X,self.p_coafs)
        p_m = U+V
        q_m = ((-1)**(self.degree%2))*U + ((-1)**(self.degree%2 + 1)) * V
        #solver

        exp_X = Ksolver(q_m,p_m).solve()
        return exp_X

class scaling_and_squaring:
    """
    class that handles the scaling and squaring algorithm and uses the previous classes to solve the linear equation systems
    to calculate the pade approxmation of the matricies.
    
    refernces:
    [2] THE SCALING AND SQUARING METHOD FOR THE MATRIX EXPONENTIAL REVISITED, Nicholas J. Higham,pages 1183-1184
    """
    def __init__(self, coafs:np.array,teta:np.array) -> None:
        self.coafs = coafs
        self.teta = teta

    def evaluate(self, X:np.array) -> np.array:
        teta = self.teta
        coafs = self.coafs
        norm1 = lambda X: np.max(np.sum(np.abs(X),axis=1))
        n = X.shape[0]
        miu = np.trace(X)/n
        X = X - miu * np.eye(n,n)
        X_norm = norm1(X)
        m = [3,5,7,9,13]
        #in case we can use exatly the correct pade approximation without the need to scale down
        for i in range(len(teta)):
            if X_norm <= teta[i]:
                return math.exp(miu)*Pade(coafs[:m[i]+1]).evaluate(X)
        s = math.ceil(math.log2(X_norm/teta[-1]))
        #scaling down
        A = X / (2**s)
        X = Pade(coafs).evaluate(A)
        #squaring 
        for i in range(s):
            X = X@X
        return math.exp(miu) * X


#%%

mat = np.loadtxt("inv_matrix(1000x1000).txt", delimiter=",", dtype=float)
exp_mat = scaling_and_squaring(b13,teta13).evaluate(mat)


#%%

mat = np.loadtxt("inv_matrix(1000x1000).txt", delimiter=",", dtype=float)
norm1 = lambda X: np.max(np.sum(np.abs(X),axis=1))
n = mat.shape[0]
s = math.ceil(math.log2(np.linalg.norm(mat,1)/teta13[-1]))
mat = mat/(2**s)
U,V = Pade(b13)._evalUV(mat,b13) 
np.savetxt("results\U.csv", U,delimiter=',')
np.savetxt("results\V.csv", V,delimiter=',')