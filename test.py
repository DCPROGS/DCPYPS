import numpy as np
import qmat
import time
import sys

def init_matrix_colqu(c_agonist, debug=False):
    """Returns the example matrix given by Colquhoun & Hawkes,
    chapter 18, p. 465, eq. 127, for a given agonist concentration.
    Example as in chapter 20, p. 593, eq. 4"""

    # Initialize rates
    beta_1 = 15.0
    beta_2 = 15000.0
    alpha_1 = 3000.0
    alpha_2 = 500.0
    k_m1 = 2000.0
    k_m2 = 2000.0
    k_p1 = 5.0e07 * c_agonist
    k_star_p2 = 5.0e08 * c_agonist
    k_p2 = 5.0e08 * c_agonist
    k_star_m2 = 1.0/3.0

    # Set up Q matrix
    Q = np.array([[             0, k_star_p2,       0,  alpha_1,    0],
                  [ 2.0*k_star_m2,         0, alpha_2,        0,    0],
                  [             0,    beta_2,       0, 2.0*k_m2,    0],
                  [        beta_1,         0,    k_p2,        0, k_m1],
                  [             0,         0,       0, 2.0*k_p1,    0]])
    Q *= 1e-3 # convert to 1/ms
    qmat.init_matrix(Q) # Update diagonal elements

    if debug: print "Initialised Colquhoun-Hawkes matrix:\n", Q

    return Q

# Below are some pure-python implementations of the corresponding
# C++ functions for benchmarking

def init_matrix(Q):
    """Updates the diagonal elements of Q in-place,
    i.e. Q will be changed by this function."""
    # Make sure that Q has square shape:
    if (not(Q.shape[0]==Q.shape[1])):
        print "Q doesn't have square shape in init_matrix(); aborting now."
        return
    for d in range(0,Q.shape[0]):
        Q[d,d]=0
        Q[d,d]=-np.sum(Q[d])

def p_inf(Q, debug=False):
    """Calculates p_inf for a Q matrix. Eq. 16, 17"""
    
    if (not(Q.shape[0]==Q.shape[1])):
        print "Q doesn't have square shape in init_matrix(); aborting now."
        return

    # add a unit vector:
    u = np.ones((Q.shape[0],1))
    S = np.concatenate((Q, u), axis=1)
    
    # Note that NumPy uses matrix multiplication for np.matrix,
    # but element-wise multiplication for np.array, so we have
    # to use np.dot here for matrix multiplication.
    p_mat = np.dot(np.transpose(u), np.linalg.inv((np.dot(S,np.transpose(S)))))
    return p_mat[0]

def mat_solve(Q, debug=False):
    """Returns the solutions to the equations given by
    Colquhoun & Hawkes, chapter 20, as a tuple:
    lambda, A. Eq. 88, 89."""

    lambda_r, X = np.linalg.eig(-Q)
    Y = np.linalg.inv(X)
    A = list()
    for i in range(0, Y.shape[0]):
        X_mat = np.empty((X.shape[0],1))
        Y_mat = np.empty((1,Y.shape[1]))
        X_mat[:,0] = X[:,i]
        Y_mat[0] = Y[i]
        A.append(np.dot(X_mat,Y_mat))

    # sort lambda and A. np.sort won't work here because
    # we need to synchronize lambda and A, so we use a
    # simple bubble sort.
    n_b = len(lambda_r)
    exchanged = True
    while exchanged and n_b >= 0:
        exchanged = False 
        for i in range(0, n_b-1):
            if (lambda_r[i] > lambda_r[i+1]):
                temp_l = lambda_r[i]
                temp_a = A[i]
                lambda_r[i] = lambda_r[i+1]
                A[i] = A[i+1]
                lambda_r[i+1] = temp_l
                A[i+1] = temp_a
                exchanged = True
        n_b -= 1

    if debug:
        for i in range(0, Y.shape[0]):
            print lambda_r[i]
            print A[i]
    return lambda_r, A

def p(t, p_0, p_inf, lambda_r, A):
    """Returns the probality of a channel being in a certain 
    state at time t. The state is characterised by p_inf (the
    probablity at equilibrium) and p_0 (the initial probability).
    The rates lambda_r and the amplitude terms A are the eigen-
    values and the eigenvectors of -Q, respectively. Eq. 24,26,27."""
    p_ret = np.empty((p_0.shape[0], len(t)))
    for j in range(0, p_0.shape[0]):
        sum_p = 0
        for i in range(1, len(A)):
            w_ij = 0
            for r in range(0, A[i].shape[1]):
                w_ij += p_0[r] * A[i][r,j]
            sum_p += w_ij*np.exp(-t*lambda_r[i])
 
        p_ret[j] = p_inf[j] + sum_p

    return p_ret

if __name__ == "__main__":

    np.set_printoptions(precision=3)

    # solution as in Ch. 20, p. 597, eq. 18
    p_inf1_ref = np.array([
            2.483e-5, 1.862e-3, 6.207e-5, 4.965e-3, 9.931e-1])
    p_inf1_hjc = np.array([
            2.483e-5, 1.862e-3, 6.207e-5, 4.965e-3, 9.931e-1])
    lambda_ref = np.array([
            0.000e0, 1.018e-1, 2.022e0, 3.094e0, 1.941e1])
    lambda_hjc = np.array([
            0.000e0, 1.018e-1, 2.022e0, 3.094e0, 1.941e1])
    Q_100 = init_matrix_colqu(100.0e-09)
    Q_10 = init_matrix_colqu(10.0e-09)
    t = np.arange(0,1000.0,0.01)

    # C++
    time0 = time.time()

    p_inf0_cpp = qmat.p_inf(Q_10)
    p_inf1_cpp = qmat.p_inf(Q_100)
    lambda_cpp, A_cpp = qmat.mat_solve(Q_100)
    y_cpp = qmat.p(t, p_inf0_cpp, p_inf1_cpp, lambda_cpp, A_cpp)

    print "\nC++ implementation took", 
    print np.round((time.time()-time0)*1e3, 3), "ms"

    # Pure Python
    time0 = time.time()

    p_inf0_py = p_inf(Q_10)
    p_inf1_py = p_inf(Q_100)
    lambda_py, A_py = mat_solve(Q_100)
    y_py = p(t, p_inf0_py, p_inf1_py, lambda_py, A_py)

    print "Python implementation took",
    print np.round((time.time()-time0)*1e3, 3), "ms"

    # Results
    sys.stdout.write("\np_inf (blue book) = %s\n" % np.array_str(p_inf1_ref, precision=3))
    sys.stdout.write("p_inf (hjcfit) =    %s\n" % np.array_str(p_inf1_hjc, precision=3))
    sys.stdout.write("p_inf (C++) =       %s\n" % np.array_str(p_inf1_cpp, precision=3))
    sys.stdout.write("p_inf (Python) =    %s\n\n" % np.array_str(p_inf1_py, precision=3))

    sys.stdout.write("lambda (blue book) = %s\n" % np.array_str(lambda_ref, precision=3))
    sys.stdout.write("lambda (hjcfit) =    %s\n"    % np.array_str(lambda_hjc, precision=3))
    sys.stdout.write("lambda (C++) =       %s\n" % np.array_str(lambda_cpp, precision=3))
    sys.stdout.write("lambda (Python) =    %s\n\n" % np.array_str(lambda_py, precision=3))
    
    try:
        import matplotlib.pyplot as plt
        plt.plot(t,y_cpp[0],t,y_py[0])
        plt.show()
    except ImportError:
        print "Matplotlib is required for plotting"
