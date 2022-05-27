import codecs
import sys
import numpy as np
import random
import math

strmap = {
    "0": "₀", "1": "₁", "2": "₂", "3": "₃", "4": "₄", "5": "₅", "6": "₆",
    "7": "₇", "8": "₈", "9": "₉",
    "lambda": "\u03BB", "alpha": "\u03B1", "beta": "\u03B2", "rho": "\u03C1",
    "Qbar" : u'Q\u0305', 'bar' : u'\u0305',
    "in" : 	u"\u2208", "pi" : u"\u03C0", "notequal" : u"\u2260",
    "tau" : u"\u03C4",
    "epsilon" : u"\u03B5", "Delta" : u"\u0394", "gamma": u"\u03B3", "Gamma": u"\u0393"
}

sval = {}
rval = {}
lambdaval = {}
rhoval = {}
alphaval = {}
betaval = {}

def nbAsStrSub(nb):
    return ''.join([strmap[c] for c in str(nb)])

def printMatrix(matrix):
    print(stringMatrix(matrix))

def stringMatrix(matrix):
    return '\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix])

def printTwoMatrix(m1,m2):
    print('\n'.join(['\t'.join([c1 + "=" + str(c2) \
     if isinstance(c1,str) else str(c1) for c1,c2 in zip(r1,r2)]) \
      for r1,r2 in zip(m1,m2)]))

def printOneCol(matrix):
    for row in matrix:
        print(row)

def printOneColTwoMatrix(m1,m2):
    for r1,r2 in zip(m1,m2):
        print(str(r1) + "\t=\t" + str(r2))

def printOneColThreeMatrix(m1,m2,m3):
    for r1,r2,r3 in zip(m1,m2,m3):
        print(str(r1) + "\t=\t" + str(r2) + "\t=\t" + str(r3))

def buildGamma(matrix,n,i):
    if i == 1:
        matrix = [[strmap["lambda"] + nbAsStrSub(1), 1], [1, strmap["rho"] + nbAsStrSub(1)]]
    else:
        matrix1 = []
        matrix2 = []
        for x in matrix:
            l1 = []
            l2 = []
            for y in x:
                l1.append(y)
                l2.append(y)
            l1.extend([strmap["lambda"] + nbAsStrSub(i), 1])
            l2.extend([1,strmap["rho"] + nbAsStrSub(i)])
            matrix1.append(l1)
            matrix2.append(l2)
        matrix = matrix1 + matrix2
    if i == n:
        return matrix
    else:
        return buildGamma(matrix,n,i+1)

def buildGammaSymbolic(n):
    matrix = []
    for i in range(1,2**n + 1):
        matrix.append([strmap["gamma"] + nbAsStrSub(i) + nbAsStrSub(j) for j in range(1, 2*n + 1)])
    return matrix

def buildX(n, terminaison=""):
    matrix = []
    for i in range(1,n+1):
        if i == 1:
            matrix = ["s" + terminaison + nbAsStrSub(1), "r" + terminaison + nbAsStrSub(1)]
        else:
            matrix.append("s" + terminaison + nbAsStrSub(i))
            matrix.append("r" + terminaison + nbAsStrSub(i))
    return matrix

def buildmatrix(matrix,letter1,letter2,n,mode,i = 1):
    if i == 1:
        if mode == "uplet":
            matrix = [(letter1 + nbAsStrSub(1),),(letter2 + nbAsStrSub(1),)]
        if mode == "func":
            matrix = [letter1 + nbAsStrSub(1),letter2 + nbAsStrSub(1)]
    else:
        matrix1 = []
        matrix2 = []
        for x in matrix:
            l1 = x
            l2 = x
            if mode == "uplet":
                l1 = l1 + (letter1 + nbAsStrSub(i),)
                l2 = l2 + (letter2 + nbAsStrSub(i),)
            if mode == "func":
                l1 = l1 + letter1 + nbAsStrSub(i)
                l2 = l2 + letter2 + nbAsStrSub(i)

            matrix1.append(l1)
            matrix2.append(l2)
        matrix = matrix1 + matrix2

    if i == n:
        return matrix
    else:
        return buildmatrix(matrix,letter1,letter2,n,mode,i+1)

def buildMR(matrix,n,i):
    if i == 1:
        matrix = [["s" + nbAsStrSub(1), 1], [1, "r" + nbAsStrSub(1)]]
    else:
        matrix1 = []
        matrix2 = []
        for x in matrix:
            l1 = []
            l2 = []
            for y in x:
                l1.append(y)
                l2.append(y)
            l1.extend(["s" + nbAsStrSub(i), 1])
            l2.extend([1,"r" + nbAsStrSub(i)])
            matrix1.append(l1)
            matrix2.append(l2)
        matrix = matrix1 + matrix2
    if i == n:
        return matrix
    else:
        return buildMR(matrix,n,i+1)

def buildBi(n,x1,y1,x2,y2):
    m1 = buildmatrix([],x1,y1,n,"uplet")
    m2 = buildmatrix([],x2,y2,n,"uplet")

    Bi = []
    for x1,x2 in zip(m1,m2):
        l = []
        for y1,y2 in zip(x1,x2):
            r = ('max', y1, y2)
            l.append(r)
        Bi.append(l.copy())
    return Bi



def buildIV(n,terminaison=""):
    matrix = []
    for i in range(1,n+1):
        if i == 1:
            matrix = [strmap["lambda"] + terminaison + nbAsStrSub(1), strmap["rho"] + terminaison + nbAsStrSub(1)]
        else:
            matrix.append(strmap["lambda"] + terminaison + nbAsStrSub(i))
            matrix.append(strmap["rho"] + terminaison + nbAsStrSub(i))
    return matrix


def buildOV(n,x1,y1,x2,y2):
    Bi = buildBi(n,x1,y1,x2,y2)
    OV = []
    for x in Bi:
        l = ("min",x)
        OV.append(l)
    return OV


def buildOVAlphaBeta(n, terminaison=""):
    return setAlphaBeta(buildOV(n,"s",
            "r" ,
             strmap["lambda"],
              strmap["rho"]), terminaison)

def setAlphaBeta(matrix,terminaison):
    nmatrix = []
    for x in matrix:
        l = len(x[1])
        i = 0
        c = ("min",)
        for z in x[1]:
            if z[1][0] == "s" and z[2][0] == strmap["lambda"]:
                c = c + (strmap["alpha"]+terminaison + z[1][1:],)
            if z[1][0] == "r" and z[2][0] == strmap["rho"]:
                c = c + (strmap["beta"]+terminaison + z[1][1:],)
            i+=1
            if i == l:
                nmatrix.append(c)
    return nmatrix

def SecondPartitionIsjTop(s,j):
    return [i for i in range(((s-1)*2**(j))+1, ((s-1)*2**(j))+ (2**(j-1))+1 )]

def SecondPartitionIsjBottom(s,j):
    return [i for i in range(((s-1)*2**(j)) + (2**(j-1))+1,  (s*2**(j))+1)]


def FirstPartitionIsj(s,j):
    return [i for i in range(((s-1)*2**(j))+1, (s*2**(j))+1)]

def FirstPartitionIsminusj(n,j):
    return [i for i in range(1,2**(n-j)+1)]

def firstPartition(n,j):
    return [FirstPartitionIsj(s,j) for s in FirstPartitionIsminusj(n,j)]

def secondPartition(n,j):
    P = []
    P.append([SecondPartitionIsjTop(s,j) for s in FirstPartitionIsminusj(n,j)])
    P.append([SecondPartitionIsjBottom(s,j) for s in FirstPartitionIsminusj(n,j)])
    return P

def PartitionsHandler(n):
    print("Partitions definitions:")
    I = [i for i in range(1,2**n+1)]
    J = [i for i in range(1,n+1)]
    print("I:" +str(I))
    print("J:" +str(J))
    for j in J:
        print("for j="+str(j)+ ", the first partition is composed by:")
        print(str(firstPartition(n,j)))
        print("for j="+str(j)+ ", the second partition is composed by:")
        print(str(secondPartition(n,j)))
        print("therefore:")
        l1 = []
        for k in secondPartition(n,j)[0]:
            for u in k:
                l1.append(strmap["gamma"] + str(u) + str(2*j-1) + "=" + strmap["lambda"]+str(j))
        l2 = []
        for k in secondPartition(n,j)[1]:
            for u in k:
                l2.append(strmap["gamma"] + str(u) + str(2*j) + "=" + strmap["rho"]+str(j))
        print(l1)
        print(l2)

def PartitionsExplainer(n):
    Gamma = buildGamma([],n,1)
    GS = buildGammaSymbolic(n)
    print(strmap["Gamma"] + strmap[str(n)] +":")
    printTwoMatrix(GS,Gamma)


def epsilonproduct(x,y):
    if x < y:
        return y
    if x >= y:
        return 0

def setInputVector():
    IV = [0.1, 1, 1, 0.8, 1, 0.3]
    lambdaj = [0.1, 1, 1]
    rhoj = [1, 0.8, 0.3]
    return [IV, lambdaj, rhoj]

def setOutputVector():
    return [0.3, 1, 0.3, 0.8, 0.3, 0.7, 0.3, 0.7]

def setVars(n, min = True):
    IV, lambdaj, rhoj = setInputVector()
    O = setOutputVector()
    if min is True:
        sol = solvingMinimalSolution(n)
    else:
        sol = solvingMaximalSolution(n)

    for i in range(1,n+1):
        sval["s" + nbAsStrSub(i)] = sol[2*(i-1)]
        rval["r" + nbAsStrSub(i)] = sol[2*(i-1)+1]

        lambdaval[strmap["lambda"] + nbAsStrSub(i)] = lambdaj[i-1]
        rhoval[strmap["rho"] + nbAsStrSub(i)] = rhoj[i-1]

        alphaval[strmap["alpha"] + nbAsStrSub(i)] =\
            max(sval["s" + nbAsStrSub(i)],lambdaval[strmap["lambda"] + nbAsStrSub(i)])
        betaval[strmap["beta"] + nbAsStrSub(i)] =\
            max(rval["r" + nbAsStrSub(i)],rhoval[strmap["rho"] + nbAsStrSub(i)])
    MR = buildMR([],n,1)
    MRval = []
    for x in MR:
        l = []
        for y in x:
            if isinstance(y, int):
                l.append(y)
                continue
            if y.startswith("s"):
                l.append(sval[y])
            else:
                l.append(rval[y])
        MRval.append(l.copy())
    IV = buildIV(n)
    IVval = []
    for x in IV:
        if x.startswith(strmap["lambda"]):
            IVval.append(lambdaval[x])
        if x.startswith(strmap["rho"]):
            IVval.append(rhoval[x])

    OV = buildOVAlphaBeta(n)
    OVval = []
    for x in OV:
        l = ()
        for y in x:
            if y.startswith("min"):
                l = l + ("min",)
            if y.startswith(strmap["alpha"]):
                l = l + (alphaval[y],)
            if y.startswith(strmap["beta"]):
                l = l + (betaval[y],)
        OVval.append(l)
    printOneColThreeMatrix(OV,OVval,O)

def yjnTop(n,j):
    O = setOutputVector()
    return max(O[i-1] for i in secondPartition(n,j)[0][0])

def yjnBottom(n,j):
    O = setOutputVector()
    return max(O[i-1] for i in secondPartition(n,j)[1][0])

def solvingMinimalSolution(n):
    IV, lambdaj, rhoj = setInputVector()
    sol = []
    for j in range(1,n+1):
        sol.append(epsilonproduct(lambdaj[j-1], yjnTop(n,j)))
        sol.append(epsilonproduct(rhoj[j-1], yjnBottom(n,j)))
    return sol





def solvingMaximalSolution(n):
    IV, lambdaj, rhoj = setInputVector()
    sol = []
    for j in range(1,n+1):
        sol.append(yjnTop(n,j))
        sol.append(yjnBottom(n,j))
    return sol



def SolvingHandler(n):
    IV = buildIV(n)
    OV = buildOV(n,"s", "r", strmap["lambda"], strmap["rho"])
    X = buildX(n)
    sol = solvingMinimalSolution(n)
    solmax = solvingMaximalSolution(n)
    print("Solving, given input data:")
    printOneColTwoMatrix(IV, setInputVector()[0])
    print("and output data:")
    printOneColTwoMatrix(OV, setOutputVector())


    print("Minimal Solution:")
    printOneColTwoMatrix(X,sol)
    setVars(n)
    print("Maximal Solution:")
    printOneColTwoMatrix(X,solmax)
    setVars(n, False)


def main(n):
    Gamma = buildGamma([],n,1)
    X = buildX(n)
    print(strmap["Gamma"] + strmap[str(n)] +":")
    printMatrix(Gamma)
    print("X:")
    printOneCol(X)

    PartitionsHandler(n)
    PartitionsExplainer(n)
    SolvingHandler(n)





if __name__ == "__main__":
    print("Let n = 3")
    n = 3
    main(n)
