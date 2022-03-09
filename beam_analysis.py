## Code adapted from Prof. Fredy Gabriel Ramírez Villanueva repository
## https://github.com/SirPrime/MatrixAnalysis-Beams.git

## Tutorial: YouTube Channel วิเคราะห์โครงสร้าง กับ อ.กิจ
## https://www.youtube.com/watch?v=hCmXwMQWafk&list=LL&index=6&t=3642s

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
#%matplotlib nbagg --> comment out if you use Jupyter Notebook
np.set_printoptions(precision=3)

#=========================================================================================
###Preprocess
#Bernoulli beam
class BeamB:
    '''We define a beam section.
     E: Modulus of elasticity
     I: Inertia of the cross section
     L: Span length'''
    def __init__(self, E, I, L):
        '''ATTRIBUTES:
             self.E: Modulus of elasticity
             self.I: Inertia of the cross section
             self.L: Span length
             self.k: stiffness matrix of the span'''
        self.E = E
        self.I = I
        self.L = L
        
        #Element stiffness matrix
        self.k = E * I / L**3 * np.array([
                [12., 6*L, -12, 6*L],
                [6*L, 4*L**2, -6*L, 2*L**2],
                [-12, -6*L, 12, -6*L],
                [6*L, 2*L**2, -6*L, 4*L**2]
            ])
#Loads Name
class Load:
    '''Clase Load'''
    def __init__(self, type):
        '''
        type = 0: Point Load
        type = 1: Distributed Load
        type = 2: Concentrated Momento
        '''
        self.type = type
    
    def type(self):
        if self.type == 0:
            print("Point Load")
        elif self.type == 1:
            print('Distributed Load')
        elif self.type == 2:
            print('Concentrated Momento')
        else:
            print('Undefined')

#Point Load
class PointLoad(Load):
    '''Point load class'''
    def __init__(self, P=0, a=0):
        '''Point load P.
         P: Load value. Positive down.
         a: Load position with respect to the left end of the section'''
        Load.__init__(self, 0)
        self.P = P
        self.a = a
    
    def __str__(self):
        return 'Point Load\n   Value= ' + str(self.P) + 'N' \
    + '\n   Position, x= ' + str(self.a) + 'm'
    
    #Qf = [Fy1, M1, Fy2, M2, ...]
    def Qf(self, L):
        '''Equivalent nodal reactions for a point Load.
         L: Beam length'''
        a = self.a
        b = L - a      
        return self.P / L**2 * np.array([
                [b**2 / L * (3*a + b)],
                [a * b**2],
                [a**2 / L * (a + 3*b)],
                [-a**2 * b]
            ])
    
    #Shear force in a section (beam without supports)
    def FQ(self, x, L):
        '''Contribution to the shear force in a section due to a point Load,
         x: position of the section considered with respect to the extreme left
         L: span length'''
        if self.a < x <= L:
            return -self.P
        else:
            return 0
         
    #Bending moment in a section (simply supported beam)
    def MF(self, x, L):
        '''Contribution to the bending moment in a section due to a punctual Load,
         x: position of the section considered with respect to the extreme left
         L: span length'''
        if 0 <= x < self.a:
            return (1 - self.a/L) * self.P * x
        elif x <= L:
            return self.a * self.P * (1 - x/L)
        else:
            return 0

#Distributed load
class DistributedLoad(Load):
    '''Distributed load class'''
    def __init__(self, q=0, a=0, l=0):
        '''distribbuted load q.
         q: load value. Positive down.
         a: distance between the left end of the span and the start of the load.
         l: length of distributed load'''
        Load.__init__(self, 1)
        self.q = q
        self.a = a
        self.l = l
    
    def __str__(self):
        return 'Load distribution\n   Value= ' + str(self.q) + 'N/m'\
    ', ' + '\n   Beginning= ' + str(self.a) + 'm' + '\n   Longitud= ' + str(self.l) + 'm'
    
    #Qf = [Fy1, M1, Fy2, M2,...]
    def Qf(self, L):
        '''Equivalent Nodal Reactions for a Load
         evenly distributed.
         L: beam length'''
        q = self.q
        a = self.a
        b = L - self.a - self.l
        return q * L / 2 * np.array([
                [1 - a/L**4*(2*L**3 - 2*a**2*L + a**3) - b**3/L**4*(2*L - b)],
                [L/6*(1 - a**2/L**4*(6*L**2 - 8*a*L + 3*a**2) - b**3/L**4*(4*L - 3*b))],
                [1 - a**3/L**4*(2*L - a) - b/L**4*(2*L**3 - 2*b**2*L + a**3)],
                [-L/6*(1 - a**3/L**4*(4*L - 3*a) - b**2/L**4*(6*L**2 - 8*b*L + 3*b**2))]
            ])
            
    #Shear force in a section (unsupported beam)
    def FQ(self, x, L):
        '''Contribution to the shear force in a section due to the distributed load.
         x: position of the section considered with respect to the extreme left
         L: Span length'''
        if self.a <= x < self.a + self.l:
            return -self.q * (x - self.a)
        elif x <= L:
            return -self.q * self.l
        else:
            return 0
    
    #Bending moment in a section (simply supported beam)
    def MF(self, x, L):
        '''Contribution to the shear force in a section due to the distributed load.
         x: position of the section considered with respect to the extreme left
         L: Span length'''
        V1 = self.q*self.l/L*(L - self.a - self.l/2)
        V2 = self.q*self.l - V1
        if 0 <= x < self.a:
            return V1 * x
        elif x <= self.a + self.l:
            return V1*x - 0.5*self.q*(x-self.a)**2
        elif x <= L:
            return V2 * (L - x)
        else:
            return 0

# Concentrated moment
class MomentConcentrated(Load):
    '''Clase momento concentrado'''
    def __init__(self, M=0, a=0):
        '''Concentrated moment M.
         M: value of the concentrated moment. Positive counterclockwise
         a: position of the moment with respect to the left end of the section'''
        Load.__init__(self, 2)
        self.M = M
        self.a = a
    
    def __str__(self):
        return 'Moment concentrate\n   Value= ' + str(self.M) + 'Nm' \
    + '\n   Posición, x= ' + str(self.a) + 'm'
    
    #[Fy1, M1, Fy2, M2,...]
    def Qf(self, L):
        '''Equivalent nodal reactions for a concentrated moment.
         L: beam length'''
        a = self.a
        b = L - a
        return self.M / L**2 * np.array([
                [-6*a*b/L],
                [b*(b - 2*a)],
                [6*a*b/L],
                [a*(a - 2*b)]
            ])
    
    #Shear force in a section (beam without supports)
    def FQ(self, x, L):
        '''Contribution to the shear force in a section due to the distributed load.
         x: position of the section considered with respect to the extreme left'''
        return 0
    
    #Bending moment in a section (simply supported beam)
    def MF(self, x, L):
        '''Contribution to the bending moment in a section due to a concentrated moment,
         These values correspond to that of a simply supported beam.
         x: position of the section considered with respect to the extreme left
         L: Span length'''
        if 0 <= x < self.a:
            return self.M / L * x
        elif self.a < x <= L:
            return self.M * (x/L - 1)
        else:
            return 0

#=========================================================================================
###Method
#displacement matrix : di = ['d1y', 'θ1', 'd2y', 'θ2', 'd3y', 'θ3',...]
def d_matrix(list_of_suport):
    d = []
    for i in range(0,len(list_of_suport)):
        dy, θ = 0, 0
        if list_of_suport[i] == 0:#displacement of fixed support
            d.append(dy)
            d.append(θ)
        elif list_of_suport[i] == 1:#displacement of pin support
            θ = 'θ'+str(i+1)
            d.append(dy)
            d.append(θ)
        elif list_of_suport[i] == 2:#displacement of free support
            dy = 'd'+str(i+1)
            θ = 'θ'+str(i+1)
            d.append(dy)
            d.append(θ)

    print(f"\nNodal displacement matrix \nd = {d}")
    return d
    
#Reaction matrix : R = ['F1y', 'M1', 'F2y', 'M2', 'F3y', 'M3',...]
def R_matrix(list_of_suport, R0):
    R = []
    for i in range(0,len(list_of_suport)):
        Fy, M = 0, 0
        if list_of_suport[i] == 0:#reaction of fixed support
            Fy = 'F'+str(i+1)
            M = 'M'+str(i+1)
            R.append(Fy)
            R.append(M)
        elif list_of_suport[i] == 1:#reaction of pin support
            Fy = 'F'+str(i+1)
            R.append(Fy)
            R.append(M)
        elif list_of_suport[i] == 2:#reaction of free support
            R.append(Fy)
            R.append(M)
    
    if len(R0) != 0:
        for i in range(0, len(R)):
            if R0[i] != 0:
                R[i] = R0[i]
      
    print(f"\nNodal reaction matrix \nR = {R}") 
    return R

#Stiffness matrix
#Assembly of the global stiffness matrix
def global_stiffness(nodes, spans, stretch):
    K = np.zeros((2*nodes, 2*nodes))
    for i in range(spans):
        K[2*i:2*i+4, 2*i:2*i+4] += stretch[i].k
    print(f"Stiffness matrix : K")
    print(f"{K}")
    
    return K

#Fixed End Force
#Local fixed end force
#Equivalent nodal reactions in each stretch
def local_FEF(spans, loads, stretch):
    QF = [0]*spans #to save the equivalent nodal reaction vectors of each stretch
    for i in range(spans): #go through all the stretches
        for j in range(len(loads[i])): #consider all the loads of each stretch
            QF[i] += loads[i][j].Qf(stretch[i].L)

    return QF

#Assembbly the global fixed end force
def global_FEF(nodes, spans, QF):
    Qf = np.zeros((2*nodes,1))
    for i in range(spans):
        Qf[2*i:2*i+4,:] += QF[i]
    print(f"\nGlobal Fixend Force, Qf :")
    print(f"{Qf}")
    return Qf

#Calculated unknown displacement
#If we know R, we don't know d.
#If we know d, we don't know R.
#[Ri] = [Ki][di]+[Qfi]
# R = ['0', '0', 'F2y', '0', 'F3y', 'M3', ...] --example
# d = ['d1', 'θ1', '0', 'θ2', '0', '0', ...] --example
# calculate d1, θ1, θ2,... exclude known == 0, 0, ...
def displacement(d, K, Qf, R):
    '''d : list of displacement vector
       K : np.array of global stiffness 
       Qf : np.array of global FEF
       R = list of nodal external force/reaction       
    '''
    #index of unknowm displacement
    J = np.where(np.array(d) != '0')[0].tolist()
    
    #index of reaction matched unknowm displacement
    I = J 

    R1 = np.zeros((len(J),1), dtype=float) #--->Matrix [Ri]

    #Assembly Ki matched index
    K1 = np.zeros((len(I), len(J)))
    for i in range(0,len(I)):
        for j in range(0,len(J)):
            K1[i][j] = K[I[i]][J[j]] #--->Matrix [Ki]

    #Assembly Qfi index
    Q01 = []
    for item in I:
        q01 = [Qf[item][0]]
        Q01.append(q01)    
    Q01 = np.array(Q01) #--->Matrix [Qfi] 
    
    #Assembly R index 
    for i in range(0, len(I)):
        R1[i][0] = (R[I[i]])

    # Calculate displacement
    # [di] = inv[Ki][Ri]+(-1*[Qfi])
    K1 = np.linalg.inv(K1) #inverse[K1]

    #unknown displacement
    di = np.dot(K1, R1+(-1*Q01)) #dot matrix

    return di #disp = m, θ = radian 

#Calculated Nodal Reaction
#[R] = [K][d] + [Qf]
def reaction(d, di, K, Qf):

    #index of d where di will be added
    ii=np.where(np.array(d) != '0')[0].tolist()

    #added di to displacement d-matrix
    for i in range(len(ii)):
        d[ii[i]] = di[i][0]

    #Convert to array(ix1)   
    dy = np.array(d).reshape(-1, 1)
    print(f"Nodal Displacement, [d] : d1, θ1, d2, θ2, ...:")
    print(f"{dy} m, radian, m, radian,...")
    
    
    #Calculated nodal reaction
    R = np.dot(K, dy) + Qf #dot matrix
    
    print(f"\nExternal Force/Nodal Reaction, [R] : F1, M1, F2, M2, ... :")
    print("[R] = [K][d] + [Qf]")
    print(f"{R/1000} kN, kN-m, kN, kN-m,...")
    return dy, R

#Calculated internal force
def internal_force(dy, b, QF, stretch):
    """
    dy : np.array of displacement 
    b : spans q'ty
    QF : local FEF
    """
    #Nodal displacements by stretch
    u = []
    for i in range(b):
        u.append(dy[2*i:2*i+4,:])
        # print(f"Local displacement: span {i+1} = {u[-1]} m, radian, m, radian,...")
    
    #Forces in each stretch
    # [Fi] = [ki][ui]+[QFi]
    F = []
    for i in range(b):
        F.append(stretch[i].k @ u[i] + QF[i])
        # print(f"Loacal force span {i+1} = {F[-1]} N, N-m, N, N-m,...")
        
    return u, F

#Shear force values
def xi_coordinate(spans, stretch):   
    #Number of sections to take for the graphics in each stretch
    numS = 1000
    Xt = [] #to save the x of each stretch
    for i in range(spans):
        Xt.append(np.linspace(0, stretch[i].L, numS)) #Sections location
    return numS, Xt

def shears(spans, stretch, loads, F):
    numS, Xt = xi_coordinate(spans, stretch)
    Shears = []
    for i in range(spans): #for each stretch

        #Shear like unsupported beams(Internal Shear)
        Q0 = np.zeros(numS) #y-y axis
        for j in range(len(loads[i])): #consider all the loads of each stretch
            m = 0 #para enumerar las secciones
            for x in Xt[i]: #to list the sections
                Q0[m] += loads[i][j].FQ(x, stretch[i].L) #Calculate Qi given xi
                m += 1

        #Shear at the extreme left, obtained from the calculation
        Q1 = F[i][0]

        #Total shear
        Shears.append(Q0+Q1)

    # Maximum and minimum shear force values (in each stretch)
    maxShear = [] #Maximum shear for each stretch
    minShear = [] #Minimal shear for each stretch
    XmaxQ= [] #locations of the maxima in each stretch
    XminQ = [] #locations of the minimum in each stretch

    print(f"\nSHEAR")
    for i in range(spans):
        maxQ = max(Shears[i]) #Máximo Shearnte
        minQ = min(Shears[i]) #Mínimo Shearnte
        print(f"Span {i+1} : maxQ = {maxQ/1000:.2f}, minQ = {minQ/1000:.2f} ,kN")

        maxShear.append(maxQ)
        minShear.append(minQ)
        indMaxQ = np.where(Shears[i] == maxQ )[0][0] #index of maximum shear
        indMinQ = np.where(Shears[i] == minQ )[0][0] #index of minimum shear
        XmaxQ.append(Xt[i][indMaxQ])#location of maximum shear
        XminQ.append(Xt[i][indMinQ])#location of minimum shear
        print(f"At location x = {Xt[i][indMaxQ]:.2f}, {Xt[i][indMinQ]:.2f} ,m")
        
    #Shear Force Values for Charts
    DFQ = []
    for i in range(spans):
        #Values for list type DFQ
        #Shear = (Shears[i]).tolist() #We go to kN and we convert to list, N
        Shear = (Shears[i]/1000).tolist() #We go to kN and we convert to list, kN
        DFQ += Shear
        
    return DFQ, maxShear, minShear, XmaxQ, XminQ

#Bending moment values
def moments(spans,  stretch, loads, F):    
    numS, Xt = xi_coordinate(spans, stretch)
    Moments = []
    for i in range(spans): #for each stretch

        #Moments like stretchs simply supported
        M0 = np.zeros(numS)
        for j in range(len(loads[i])): #consider all the loads of each stretch
            m = 0 #to list the sections
            for x in Xt[i]: #go through the sections
                M0[m] += loads[i][j].MF(x, stretch[i].L) 
                m += 1

        #Moments due to embedment or continuity of the beam
        M1 = -F[i][1] + (F[i][3] + F[i][1]) / stretch[i].L * Xt[i]

        #Total moment
        Moments.append(M0 + M1)

    #Maximum and minimum bending moment values (in each stretch)
    maxMoment = [] #Maximum moment in each stretch
    minMoment = [] #Minimum moment in each stretch
    XmaxF= [] #locations of maximum moments by stretch
    XminF = [] #locations of the minimum moments by stretch
    print(f"\nMOMENT")
    for i in range(spans):
        maxF = max(Moments[i]) #Máximo flector
        minF = min(Moments[i]) #Mínimo flector
        print(f"Span {i+1} : maxF = {maxF/1000:.2f}, minF = {minF/1000:.2f} ,kN-m")

        maxMoment.append(maxF)
        minMoment.append(minF)
        indMaxF = np.where(Moments[i] == maxF )[0][0] #index of maximum bending
        indMinF = np.where(Moments[i] == minF )[0][0] #index of minimum bending    
        XmaxF.append(Xt[i][indMaxF])#location of maximum bending
        XminF.append(Xt[i][indMinF])#location of minimum bending 
        print(f"At location x = {Xt[i][indMaxF]:.2f}, {Xt[i][indMinF]:.2f} ,m")
        
    #Bending moment values for graphs
    DMF = []
    for i in range(spans):
        #Values for list type DMF
        #Flex = (Moments[i]).tolist() #We go to kNm and convert to list, N-m
        Flex = (Moments[i]/1000).tolist() #We go to kNm and convert to list. kN-m
        DMF += Flex
        
    return DMF, maxMoment, minMoment, XmaxF, XminF

#=========================================================================================
###Diagram
#Values of x for graphs
def X_coordinate(spans, stretch, Xt):
    X = []
    Lacum = 0
    for i in range(spans):
        if i > 0:
            Lacum += stretch[i-1].L
        Xprov = Xt[i] + Lacum
        Xlist = Xprov.tolist()
        X += Xlist
    return X

def plot(lb, spans, Ltotal, stretch, DFQ, maxShear, minShear, XmaxQ, XminQ):
    numS, Xt = xi_coordinate(spans, stretch)
    X = X_coordinate(spans, stretch, Xt)
    
    #For N, N-m
#    lb_Q = ['SFD', 'Shear Force [N]', '$Q_{max} = $', '$Q_{min} = $', '$N, x= $']
#    lb_M = ['BMD', 'Moment [N-m]', '$M_{max} = $', '$M_{min} = $', '$N-m, x= $']
    
    #For kN, kN-m
    lb_Q = ['SFD', 'Shear Force [kN]', '$Q_{max} = $', '$Q_{min} = $', '$kN, x= $']
    lb_M = ['BMD', 'Moment [kN-m]', '$M_{max} = $', '$M_{min} = $', '$kN-m, x= $']
    
    if lb=='Shear':
        lb = lb_Q
        # plt.figure(1)
    else:
        lb = lb_M
        # plt.figure(2)
        plt.gca().invert_yaxis() #invert y-axis --> M+ == bottom
        
    #Shear Force Diagram    
    plt.plot(X, DFQ)
    plt.title(lb[0], fontsize = 12)
    plt.xlabel('x [m]')
    plt.ylabel(lb[1])
    plt.axhline(linewidth = 3)
    plt.xlim(0, Ltotal)
    plt.grid()

    #Texts for maximum and minimum values
    def colocarTextosQ(lb):
        LacumQ = 0
        for i in range(spans):
            if i > 0:
                LacumQ += stretch[i-1].L
            ubicMax = LacumQ + XmaxQ[i]
            ubicMin = LacumQ + XminQ[i]
            if ubicMax == Ltotal:
                ubicMax = Ltotal - stretch[i].L/2
            if ubicMin == Ltotal:
                ubicMin = Ltotal - stretch[i].L/2

            ##Plot For unit in N, N-m
#             plt.text(ubicMax, maxShear[i], lb[2] + \
#                      str(round(maxShear[i],2)) + lb[-1] + str(round(ubicMax,2)) \
#                      + '$m$')
#             plt.text(ubicMin, minShear[i], lb[3] + \
#                      str(round(minShear[i],2)) + lb[-1] + str(round(ubicMin,2)) \
#                      + '$m$')
                     
           #Plot For unit in kN, kN-m
            plt.text(ubicMax, maxShear[i]/1000*1.1, lb[2] + \
                     str(round(maxShear[i]/1000,2)) + lb[-1] + str(round(ubicMax,2)) \
                     + '$m$')
            plt.text(ubicMin, minShear[i]/1000*1, lb[3] + \
                     str(round(minShear[i]/1000,2)) + lb[-1] + str(round(ubicMin,2)) \
                     + '$m$')

    colocarTextosQ(lb)

    #To shade the graph.
    Xgraf = [0] + X
    Xgraf.append(Ltotal)

    DFQgraf = [0] + DFQ
    DFQgraf.append(0)

    plt.fill(Xgraf, DFQgraf, 'b', alpha=0.3)

    #Stretch dividers
    vertical = 0
    for i in range(spans - 1):
        vertical += stretch[i].L
        plt.axvline(vertical, color='black')

    # plt.show()

#=========================================================================================
####
def main(E, I, spans, s, loads, R0):
    print("BEAM ANALYSIS : METRIX STIFFNESS METHOD")
    print("Code adapted from Prof. Fredy Gabriel Ramírez Villanueva repository")
    print("https://github.com/SirPrime/MatrixAnalysis-Beams.git")
    print("")
    print("Tutorial: YouTube Channel 'วิเคราะห์โครงสร้าง กับ อ.กิจ'")
    print("https://www.youtube.com/watch?v=hCmXwMQWafk&list=LL&index=6&t=3642s")
    print("=========================================================================================")

    ##User input
    print("PROPERTIES :")
    #https://en.wikipedia.org/wiki/Young%27s_modulus
    # E = 200e9 # Es, GPa
    # https://optimalbeam.com/section-properties.php 
    # I = 252*10e-8 #m4
    print(f"E = {E} GPa, I = {I} m4")

    #----------------------------------------------------
    print(f"\nGEOMETRY :")
    #Define length of each span
    ##User input
    # spans = [6, 1]
    for i in range(len(spans)):
        print(f"Span {i+1} : {spans[i]} m")

    #----------------------------------------------------
    print(f"\nSUPPORT :")
    support = {
        "0" : "Fixed",
        "1" : "Pin",
        "2" : "Free"
    }
    #Define list of support type
    ##User input
    # s = [1, 1, 2] 
    leftSupport = s[0]
    rightSupport = s[-1]
    for i in range(len(spans)+1):
        print(f"Support {i+1} = {support.get(str(s[i]))}") 

    #----------------------------------------------------
    print(f"\nNODAL EXTERNAL FORCE, Ro :")
    #Define known/unknown vector of external force/reaction (+Up, -Down)
    #R0 = ['F1y', 'M1', 'F2y', 'M2', 'F3y', 'M3',...]
    ##User input
    # R0 = [0, 0, 0, 0, 0, 0] #N, N-m
    print(f"[Ro] = ['F1y', 'M1', 'F2y', 'M2', 'F3y', 'M3',...]")
    print(f"[Ro] = {R0}  N, N-m...")

    #----------------------------------------------------
    print(f"\nBERNOULLI BEAM :")
    #Define the stretchs of the continuous beam in a list
    #BeamB(Elasticity, Inertia, Length) for each stretch
    stretch = []
    for i in range(len(spans)):
        st = BeamB(E, I, spans[i]) 
        print(f"K{i+1}")
        print(f"{st.k}")
        stretch.append(st)

    #----------------------------------------------------
    print(f"\nLOAD:")
    #Define loads in each stretch
    #q = DistributedLoad (value, start, length), distance between the left end of the span and the start of the load, Down+ Up-
    #P = PointLoad(value, position), Load position with respect to the left end of the section, Down+ Up-
    #M = MomentConcentrated (value, position), position of the moment with respect to the left end of the section', counterclockwise+

    ##User input
    #provide load: unit in --> Newton, N
    
    #print load
    for i in range(0, len(loads)):
        print(f'Load in stretch {i+1} : ')
        print(*loads[i], sep = "\n")

    #Number of stretchs or bars
    b = len(stretch)

    #Number of nodes
    nodes = b + 1

    #Total length of the beam
    Ltotal = 0
    for i in range(b):
        Ltotal += stretch[i].L

    #----------------------------------------------------
    print(f"\nCALCULATION :")
    #Assembly Global Stiffness, K
    K = global_stiffness(nodes, b, stretch)

    #Calculate and assembly Fixed End Force (FEF)
    QF = local_FEF(b, loads, stretch)
    Qf = global_FEF(nodes, b, QF)

    #Create displacement matrix & reaction matrix
    Ro = R_matrix(s, R0) #list
    d = d_matrix(s) #list

    #Calculate unknown displacement
    di = displacement(d, K, Qf, Ro)

    print(f"\nMerge displacement and calculate reaction :")
    dy, R = reaction(d, di, K, Qf)

    #Calculate local displacement and local force
    u, F = internal_force(dy, b, QF, stretch) 
    
    #----------------------------------------------------
    #Calculate shears coordinate for plotting
    DFQ, maxShear, minShear, XmaxQ, XminQ = shears(b, stretch, loads, F)

    ###### Calculate moments coordinate for plotting
    DMF, maxMoment, minMoment, XmaxF, XminF = moments(b, stretch, loads, F)

    #plot diagram 
    plt.figure(1)
    #plot
    plt.subplot(211) #2 rows, 1 column, first position
    plt.margins(2)
    plot('Shear', b, Ltotal, stretch, DFQ, maxShear, minShear, XmaxQ, XminQ)
    
    #plot
    plt.subplot(212)  #2 rows, 1 column, second position
    plt.margins(2) 
    plot('Moment', b, Ltotal, stretch, DMF, maxMoment, minMoment, XmaxF, XminF)

    plt.subplots_adjust(hspace=0.5)
    plt.show()

    ask = input('Save image? Y|N : ').upper()
    if ask == 'Y':
        plt.figure().savefig('Data/test.png')

    print("=========================================================================================")
    print("If you have any comment, pls contact me at highwaynumber12@gmail.com")

#=========================================================================================
### PLAY GROUND ###
##  Example from URL : https://learnaboutstructures.com/sites/default/files/images/3-Frames/Det-Beam-Example-Moment.png 
##  Recheck :  https://platform.skyciv.com/beam
### comment out all below to run code ###
'''
E = 200e9 #Es, GPa 
I = 67500e-8 #m4

spans = [10.0, 5.0]
s = [1, 1, 2]
R0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Load in stretch 1
P1 = PointLoad(20000, 5)#, Down+ Up-
m1 = MomentConcentrated(-30000, 5)#counterclockwise +

# Load in stretch 2
q2 = DistributedLoad(10000, 0, 5)#, Down+ Up-

f1 = [P1, m1]
f2 = [q2]
loads = [
    f1,
    f2
    ]

if __name__ == '__main__':
    main(E, I, spans, s, loads, R0)
'''