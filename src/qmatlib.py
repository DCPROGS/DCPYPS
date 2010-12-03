#! /usr/bin/python

__author__="RemisLape"
__date__ ="$11-Oct-2010 10:33:07$"

import math
import numpy as np
from numpy import linalg as nplin

def getGs(Q, kA, kB, debug=False):
    """Calculate GBA and GAB matrices.
       GBA=-QBB^(-1)*QBA
       GAB=-QAA^(-1)*QAB
    """
    kE = kA + kB

    QBB = Q[kA:kE, kA:kE]
    QBA = Q[kA:kE, 0:kA]
    QAA = Q[0:kA, 0:kA]
    QAB = Q[0:kA, kA:kE]

    GAB = np.dot(nplin.inv(-1 * QAA), QAB)
    GBA = np.dot(nplin.inv(-1 * QBB), QBA)

    if debug: print 'GAB= ', GAB
    if debug: print 'GBA= ', GBA
    return GAB, GBA

def getEigs(Q, k=None, debug=False):
    """
    Calculate eigenvalues, eigenvectors and spectral matrices.
    Return eigenvalues and spectral matrices.
    
    Parameters
    ----------
    Q : array_like, shape (k, k)
    k : {non-zero int, None}
        size of `Q`

    Returns
    -------
    eigvals : ndarray, shape (k,)
    A : list of ndarrays with shape (k, k); length `k`
        Spectral matrices of Q.
    """

    eigvals, M = nplin.eig(Q)
    if debug: print 'eigenvalues=', eigvals
    if debug: print 'eigenvectors=', M
    N = nplin.inv(M)

    if k is None or k<0 or k>N.shape[0]:
        k = N.shape[0]

    A = []
    for i in range(0, k):
        Ai = np.zeros((k, k))
        for j in range(0, k):
            for m in range(0, k):
                Ai[j, m] = M[j, i] * N[i, m]
        A.append(Ai)
    if debug: print 'spectral matrices =', A

    return eigvals, A

def getExpQsub(Q, k, t, debug=False):
    """Calculate exponential of a matrix.
    Q = exp(Q * t)
    """

    eigvals, A = getEigs(-Q, k, debug)
    expQ = np.zeros((k, k))
    for i in range(0, k):
        for j in range(0, k):
            for m in range(0, k):
                temp = A[m][i, j] * math.exp(-eigvals[m] * t)
                expQ[i, j] = expQ[i, j] + temp

    return expQ

def getEGs(G12, G21, k1, k2, expQ22):
    """Calculate eGAF or eGFA for calculation of initial vectors.
    (I - GAB * (I - expQBB) * GBA)^-1 * GAB * expQBB
    """

    temp = np.eye(k1) - np.dot(np.dot(G12, np.eye(k2) - expQ22), G21)
    eG12 = np.dot(np.dot(nplin.inv(temp), G12), expQ22)

    return eG12

def phiHJC(eG12, eG21, k1, k2, debug=False):
    """Solves phi*(I-eG12*eG21)=0"""

    if k1 == 1:
        pinf = 1
        return pinf

    Qsub = np.eye(k1) - np.dot(eG12, eG21)
    u = np.ones((k1,1))
    S = np.concatenate((Qsub, u), 1)
    pinf = np.dot(u.transpose(), nplin.inv(np.dot(S, S.transpose())))

    return pinf

def dARSdS(tres, Q11, Q12, Q22, Q21, G12, G21, expQ22, expQ11, k1, k2):
    """
    Python implementation of DC's DARSDS subroutine inside HJCMEAN.FOR.
    Evaluates [-dAR(s)/ds at s=0 for means etc of HJC distributions.
    Output is DARS- kA x kA matrix.
    SFF = I - exp(QFF * tres)
    First evaluate [dVA(s) / ds] * s = 0.
    Q1 = -inv(QAA) * GAF*SFF*GFA - GAF*SFF* inv(QFF)* GFA +
    + tres*GAF* expQFF*GFA
    CORRECTION by DC (26-May-92): sign of last term should be (-)

    Then: DARS = inv(VA)*(QAA^-2) - inv(VA)*Q1*inv(VA)*inv(QAA =
    = inv(VA) * [inv(QAA) - Q1*inv(VA)] * inv(QAA)
    where VA = I - GAF * SFF * GFA
    """

    invQ11 = nplin.inv(Q11)
    invQ22 = nplin.inv(Q22)

    #SFF = I - EXPQF
    I = np.eye(k2)
    SFF = I - expQ22

    #Q1 = tres * GAF * exp(QFF*tres) * GFA
    Q1 = tres * np.dot(G12, np.dot(expQ22, G21))
    #Q2 = GAF * SFF * inv(QFF) * GFA
    Q2 = np.dot(G12, np.dot(SFF, np.dot(invQ22, G21)))
    #Q3 = -inv(QAA) * GAF * SFF * GFA
    Q3 = np.dot(np.dot(np.dot(-invQ11, G12), SFF), G21)
    Q1 = Q1 - Q2 + Q3

    # VA = I - GAF * SFF * GFA
    I = np.eye(k1)
    VA = I - np.dot(np.dot(G12, SFF), G21)

    # DARS = inv(VA) * (QAA**-2) - inv(VA) * Q1 * inv(VA) * inv(QAA) =
    #      = inv(VA) * [inv(QAA) - Q1 * inv(VA)] * inv(QAA)

    Q3 = invQ11 + - np.dot(Q1, nplin.inv(VA))
    DARS = np.dot(np.dot(nplin.inv(VA), Q3), invQ11)

    return DARS

def pinf(Q):
    """Calculate ecquilibrium occupancies by adding a column of ones
    to Q matrix. Pinf = uT*invert((S*transpos(S))).
    """

    u = np.ones((Q.shape[0],1))
    S = np.concatenate((Q, u), 1)
    pinf = np.dot(u.transpose(), nplin.inv((np.dot(S,S.transpose()))))
    return pinf[0]

def phiO(Q, kA, kB, kC, debug=False):
    """
    Calculate initial vector for openings.
    """

    k = kA + kB + kC
    uA = np.ones((kA,1))
    p = pinf(Q)
    pF = p[kA:k]
    QFA = Q[kA:k, 0:kA]
    nom = np.dot(pF, QFA)
    denom = np.dot(nom,uA)
    phiO = nom / denom

    if debug: print 'phiO=', phiO
    return phiO

def phiS(Q, kA, kB, kC, debug=False):
    """
    Calculate inital vector for shuttings.
    """

    kF = kB + kC
    phiOp = phiO(Q, kA, kB, kC)
    GAF, GFA = getGs(Q, kA, kF)
    phiS = np.dot(phiOp, GAF)

    if debug: print 'phiS=', phiS
    return phiS
