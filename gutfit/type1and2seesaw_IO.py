import numpy as np
import random
import math
from random import randint, choice
from gutfit import model, parameterlist

#Debug#
testUPMNS = False
testf = False
test1 = False
testyu = False
testYl = False

testYnu = False
testlepto = False
testneutrino = False
testbenchmark = False
testsign = False
testPred = False
#r1,r2 log scan, cecnu mo, log as well

#########################################################################################################################################################################################################################
#INVERTED ORDERING WHEN THIS VARIABLE IS TRUE

IO = False

##############################àààà##############################################à


def classify_number(x):
    return 1 if x < 0.5 else -1

def matrix_diag3(d1,d2,d3):
    return np.array([[d1, 0.0, 0.0], [0.0, d2, 0.0], [0.0, 0.0, d3]], dtype = np.complex64)

# Generic Rotations #
def matrix_rot23(th23):
    return np.array([[1.0,          0.0 , 0.0],
                    [0.0,  np.cos(th23), np.sin(th23)],
                     [0.0, -np.sin(th23), np.cos(th23)]],dtype = np.complex64)

def matrix_rot12(th12):
    return np.array([[ np.cos(th12), np.sin(th12), 0.0],
                    [-np.sin(th12), np.cos(th12), 0.0],
                     [          0.0,  0.0,         1.0]], dtype = np.complex64)

def matrix_rot13(th13, delta):
    return np.array([[ np.cos(th13), 0.0, np.sin(th13) * np.exp(-1j * delta)],
                    [ 0.0, 1.0, 0.0                               ],
                    [-np.sin(th13)* np.exp(1j * delta), 0.0, np.cos(th13)]],
                    dtype=np.complex64)

def matrix_vckm(th12, th13, th23, delta):
    return matrix_rot23(th23) @ matrix_rot13(th13, delta) @ matrix_rot12(th12)

# Phase Matrices #

def matrix_phase(a1, a2, a3):
    return np.array([[np.exp(1j * a1),             0.0,             0.0],
                    [            0.0, np.exp(1j * a2),             0.0],
                    [            0.0,             0.0, np.exp(1j * a3)]],
                    dtype=np.complex64)

def matrix_phase2(a1, a2):
    return np.array([[np.exp(1j * a1),             0.0,             0.0],
                    [            0.0, np.exp(1j * a2),             0.0],
                     [            0.0,             0.0,             1.0]],
                        dtype=np.complex64)


class Type1And2SeeSaw_v4(model.Model):
    """The summary line for a class docstring should fit on one line.

If the class has public attributes, they may be documented here
in an ``Attributes`` section and follow the same formatting as a
function's ``Args`` section. Alternatively, attributes may be documented
inline with the attribute's declaration (see __init__ method below).

Properties created with the ``@property`` decorator should be documented
in the property's getter method.

Attributes:
    attr1 (str): Description of `attr1`.
    attr2 (:obj:`int`, optional): Description of `attr2`.

"""
    
    
   
    

    
    def __init__(self,IO= False):
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            param1 (str): Description of `param1`.
            param2 (:obj:`int`, optional): Description of `param2`. Multiple
                lines are supported.
            param3 (:obj:`list` of :obj:`str`): Description of `param3`.

        """
        #self.attr1 = param1
        #self.attr2 = param2
        #self.attr3 = param3  #: Doc comment *inline* with attribute

        #: list of str: Doc comment *before* attribute, with type specified
        #self.attr4 = ['attr4']

        #self.attr5 = None
        """str: Docstring *after* attribute, with type specified."""

        params = ["sign_parameter1",
           "sign_parameter2",
           "sign_parameter3",
           "sign_parameter4",
           "sign_parameter5",
           "generic_quark_phase_a1",
           "generic_quark_phase_a2",
           "data_quark_th12",
           "data_quark_th13",
           "data_quark_th23",
           "data_quark_delta",
           "data_quark_yu",
           "data_quark_yc",
           "data_quark_yt",
           "data_quark_yd",
           "data_quark_ys",
           "data_quark_yb",
           "model4_mR",
           "model4_r1",
           "model4_r2",
           "model4_cnu",
           "model4_ce"
           ]
        super().__init__(params)
        self.IO = IO
        #self.randomphase(self.sign_parameter1, self.sign_parameter2, self.sign_parameter2, self.sign_parameter3, self.sign_parameter4, self.sign_parameter5)
        #self.randomsign(self.sign_parameter)
    
    
    @property
    def val(self):
        '''
        

        Returns:
            value (TYPE): DESCRIPTION.

        '''
        self.randomphase(self.sign_parameter1,self.sign_parameter2, self.sign_parameter3, self.sign_parameter4, self.sign_parameter5)
        value = np.abs(
                self.PRED(
                    self.sign_parameter1,
		            self.sign_parameter2, 
		            self.sign_parameter3,
		            self.sign_parameter4, 
		            self.sign_parameter5, 
                    self.generic_quark_phase_a1,
                    self.generic_quark_phase_a2,
                    self.data_quark_th12,
                    self.data_quark_th13,
                    self.data_quark_th23,
                    self.data_quark_delta,
                    self.data_quark_yu,
                    self.data_quark_yc,
                    self.data_quark_yt,
                    self.data_quark_yd,
                    self.data_quark_ys,
                    self.data_quark_yb,
                    self.model4_mR,
                    self.model4_r1,
                    self.model4_r2,
                    self.model4_cnu,
                    self.model4_ce
                    )
                )
        return value

   
    def randomphase(self,s1,s2,s3,s4,s5):  #reparametrize sign for quark yukawa couplings
        
    
        self.ydrand  = classify_number(s1)
        self.ysrand  = classify_number(s2)
        self.ybrand  = classify_number(s3)
        self.yurand  = classify_number(s4)
        self.ycrand  = classify_number(s5)
        self.ytrand  = 1
        
   
    
    def matrix_Yd(self,s1,s2,s3,s4,s5, a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce):
        Pa      = matrix_phase2(a1, a2)
        Vckm    = matrix_vckm(th12q, th13q, th23q, deltaq)
        ydrand  = self.ydrand*yd
        #print(self.ydrand)
        #print(s)
        ysrand  = self.ysrand*ys
        ybrand  = self.ybrand*yb
        Yddiag  = matrix_diag3(ydrand, ysrand, ybrand)
        Vckmc   = np.conjugate(Vckm)
        Yukd    = Pa @ Vckm @ Yddiag  @ np.transpose(Vckmc) @  np.conjugate(Pa)
        
      #  print(Yddiag)
       # print("Here Yd")
       # print(Yd)
        if (test1==True):
            print("Yd")
           # print()                                                                                                                                                                                       
            print(Yukd)
        return  Yukd

#Modify the function below for a different model, at the moment used for 2209.00021 model                                                                                                                  

    def MnuTheory(self,s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce):
        Yd        = self.matrix_Yd(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        Yu        = matrix_diag3(self.yurand*yu, self.ycrand*yc,  self.ytrand*yt)
        ReYd      = np.real(Yd)
        ImYd      = np.imag(Yd)
        cnulogged = cnu
        r2logged  = r2
        r1logged  = r1
        ydrand    = yd
        ysrand    = -ys
        ybrand    =  yb
        type1p1   = (8 * r2logged * (r2logged+1) * Yu)/(r2logged-1)
        type1p2   = -(16 * r2logged*r2logged * ReYd)/(r1logged * (r2logged-1))
        type1p3   = ((r2logged-1)/r1logged) * (r1logged * Yu + 1j * cnulogged * ImYd) @ np.linalg.inv(r1logged * Yu - ReYd) @ (r1logged * Yu - 1j * cnulogged * ImYd)
        type1     = (mR) * (type1p1 + type1p2 + type1p3)
        if (testyu == True):
            #print("Mnu")                                                                                                                                                                                  
            print(Yu)                                                                                                                                                                                     
            #print(self.sign_parameter)
        return  type1

    def YlTheory(self,s1,s2,s3,s4,s5,a1,a2,th12q,th13q,th23q,deltaq,yu,yc,yt,yd,ys,yb,mR, r1, r2, cnu,ce):
        Yd = self.matrix_Yd(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        Yu = matrix_diag3(self.yurand*yu,  self.ycrand*yc,  self.ytrand*yt)
        ReYd = np.real(Yd)
        ImYd = np.imag(Yd)
        celogged = ce
        r2logged  = r2
        r1logged  = r1
        ydrand    = yd
        ysrand    = -ys
        ybrand    =  yb
        part1 = -(4*r1logged*Yu)/(r2logged-1)
        part2 = ((r2logged+3)*ReYd)/(r2logged-1)
        part3 = 1j*celogged*ImYd
        Ylepton = part1+part2+part3
        #print(Ylepton)
        if (testYl==True):
            print("Ye")
            print(Ylepton)
        return Ylepton


    def MixingMatrix(self,  s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR,r1, r2, cnu, ce): #diagonalize both Mnu and Yl, use the method linalg.eigh, ordered eigenvalues             
        Mnu = self.MnuTheory(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)*(10**9)
        Ylepton = self.YlTheory(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)
       # leptonmasses, Ul = np.linalg.eig(Ylepton)                                                                                                                                                         
        leptomasssquared, Vl = np.linalg.eigh(Ylepton @ np.transpose(np.conjugate(Ylepton)))
        if (testlepto==True):
            print(np.transpose(np.conj(Vl)) @ Ylepton @ np.transpose(np.conj(Ylepton)) @ Vl)
            print(np.transpose(np.conj(Vl)) @ Vl)
        neutrinomasses, Vnu = np.linalg.eigh(Mnu @ np.transpose(np.conj(Mnu)))
        if (testneutrino == True):
            print(Mnu)
        return abs(leptomasssquared), abs(neutrinomasses), Vl,Vnu



    def UPMNS(self, s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce):
        lm, nm, Vl, Vnu = self.MixingMatrix(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        if (testUPMNS == True):
            print(UMPNS(PL()))
        return  np.transpose(np.conj(Vl)) @ Vnu

    def angles(self, s1,s2,s3,s4,s5,a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce):
        U = self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)
        sin13 = abs(self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)[0][2])
        tan12 = abs(self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)[0][1]/self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)[0][0])
        tan23 = np.abs(self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)[1][2]/self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)[2][2])
        delta = np.angle(U[0][1]*U[1][2]*np.conj(U[0][2])*np.conj(U)[1][1]+(np.sin(np.arctan(tan12))**2)*(np.cos(np.arcsin(sin13))**2)*(np.sin(np.arctan(tan23))**2)*(sin13**2))
        aM1 = np.angle(U[0][1])
        aM2 = np.angle(U[1][2])
        return sin13, tan12, tan23, delta, aM1, aM2

    def masses(self, s1,s2,s3,s4,s5,a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce):
        Mnu = self.MnuTheory(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR, r1, r2, cnu, ce)*(10**9)
        nm, Vnu = np.linalg.eigh(np.transpose(np.conjugate(Mnu))@Mnu)
        electron_ordering = np.argsort(np.abs(Vnu[:,0]))
        nm = nm[electron_ordering]
        if self.IO:
            print("MassesIO")
            print(nm[0]+ nm[1]+nm[2])
            return np.abs(nm[0]-nm[1]), -np.abs(nm[1]-nm[2])
        else:
            print(nm[0]+ nm[1]+nm[2])
            return np.abs(nm[0]-nm[1]), np.abs(nm[0]-nm[2])



    def RHNMatrix(self,s1,s2,s3,s4,s5, a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce):
        Yd = self.matrix_Yd(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        Yu = matrix_diag3(self.yurand*yu,  self.ycrand*yc,  self.ytrand*yt)
        ReYd = np.real(Yd)
        ImYd = np.imag(Yd)
        celogged = ce
        r2logged  = r2
        r1logged  = r1
        ydrand    = yd
        ysrand    = -ys
        ybrand    =  yb
        vSM = 246 #GeV                                                                                                                                                                                     
        f = Yu*(1/(r2logged-1))-ReYd*(1/(r1logged*(r2logged-1)))
        if (testf == True):
            print(f)
        MN = (f*(vSM**2))/(mR)  #Should be in GeV now                                                                                                                                                      
        if (testYnu == True):
            print(MN)
        return MN



    def PRED(self,s1,s2,s3,s4,s5,a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce): #return observables needed for the scan
        
        U = self.UPMNS(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR,r1, r2, cnu, ce)
        lm,nm, Vl, Vnu = self.MixingMatrix(s1,s2,s3,s4,s5,a1, a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb, mR,r1, r2, cnu, ce)
        msq12, msq13 = self.masses(s1,s2,s3,s4,s5,a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        yyye = np.sqrt(lm[0])
        yyymu = np.sqrt(lm[1])
        yyytau = np.sqrt(lm[2])
        s13, t12, t23,d, aM1,aM2 = self.angles(s1,s2,s3,s4,s5,a1 , a2,  th12q, th13q, th23q, deltaq, yu, yc, yt, yd, ys, yb,mR, r1, r2, cnu, ce)
        A = []
        A.append(np.arcsin(s13))
        A.append(np.arctan(t12))
        A.append(np.arctan(t23))
        A.append(msq12)
        A.append(msq13)
        A.append(yyye)
        A.append(yyymu)
        A.append(yyytau)
        A = np.asarray(A)
        return A
     
def readparameter(fname):  #reads the parameter card 
    d = {}
    with open(fname) as f:
        for line in f:
            l = line.strip()
            if len(l) == 0 or l.startswith("#"):
                continue
                print(l)

            k,v  = l.split(" ")
            d[k] = float(v)
    return d

def sample(x,E,PL,N,pnames): #Creates a distribution with N points for theory and data 
    D = []
    T = []
    trials=0
    for _ in range(N):
        p=PL()
        for num, pn in enumerate(pnames):
            p[pn] = x[num]
            print(p[pn])
            print(x[num])
        print(p)
        data = E(p)
        print(data)
        D.append(data)
        
    return D, True

if __name__=="__main__":
    import optparse, os, sys
    op = optparse.OptionParser(usage=__doc__)
    opts, args = op.parse_args()
    PL = parameterlist.ParameterList.fromConfigFile(args[0])
    N  = 1 # number of samples of PL
    pdict        =   readparameter(args[1])
    usethese = [] #it select free parameters, look at parameterlist.py for understand better what happens 

    bounds, pnames = PL.getBox(usethese)

    PMIN   = [b[0] for b in bounds]
    PMAX   = [b[1] for b in bounds]
    PLEN   = [PMAX[i] - PMIN[i] for i in range(len(pnames))]
    E =  Type1And2SeeSaw_v4()
    #print(sample([pdict[pn] for pn in pnames] ,E,PL,N, pnames))
    E.ydrand = -1
    E.ysrand = +1
    E.ybrand = +1

    E.yurand = +1
    E.ycrand = -1
    E.ytrand = -1
    
    
    '''
    PL =  parameterlist.ParameterList.fromConfigFile("../examples/param_card3sigma_newsmall.dat")
#    E.sign_parameter = 0.078125
    # from IPython import embed
   # embed()
    #print(E.OBS(PL()))

    if (testbenchmark):
        import time
        t0 = time.time()
        
      #  print("theta13_bench = 0.1489607684194921          th: %f" % E(PL())[0])
       # print("theta12_bench = 0.5712806657423788          th: %f" % E(PL())[1])
        #print("theta23_bench = 0.732017203612924         th: %f" % E(PL())[2])
       # print("delta_bench = -2.199033363780465         th: %f" % E(PL())[3])
       # print("a1 th: %f" % E(PL())[4])
       # print("a2 th: %f" % E(PL())[5])
       # print("m1_bench = 3.3594212690227954           th: %f" % E(PL())[6])
       # print("mbetabeta_bench = 5.830868469833373     th: %f" % E(PL())[7])
       # print("msq21 : %f" % E(PL())[8])
       # print("whatever : %f" % E(PL())[4])
       # print("me : %f" % E(PL())[10])
       # print("mmuon : %f" % E(PL())[11])
       # print("mtau : %f" % E(PL())[12])
        #print("m1_bench = th: %f" % E(PL())[6])
       # print("MN1_bench = 4.2292070742425476e11          th: %f" % E(PL())[13])
        #print("MN2_bench = 5.324270909148878e111          th: %f" % E(PL())[14])
        #print("MN3_bench = 1.6654004807249688e13          th: %f" % E(PL())[15])
        # print("m1_bench = 3.3594212690227954          th: %f" % E(PL())[11])
        # print("mbb_bench = 5.830868469833373          th: %f" % E(PL())[12])
        
       # print(time.time() - t0)

    if testsign:
        A = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
        for i in range(32):
            #print("comb")
            E.randomsign(float(i)/32)
            #E(PL())[4]
    '''
