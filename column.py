###COLUMN DESIGN : USD METHOD
#Credit Article: "https://structurepoint.org/pdfs/Interaction-Diagram-Tied-Reinforced-Concrete-Column-Symmetrical-ACI318-14.htm"
'''
What 's this program do?
-render interaction diagram associated with section and steel which user input
-render external force
-so you can see external force point(Mu, Pu) locate in IR-daigram 
-and adjust input section or reinforcment if you want for design.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
matplotlib.use('Agg')

from tabulate import tabulate

#Use TeX in graph
matplotlib.rcParams.update(
    {
        'text.usetex': False,
        'font.family': 'stixgeneral',
        'mathtext.fontset': 'stix',
    }
)
np.set_printoptions(precision=4)

#--------------------------------------------------------------
εc = 0.003 #concrete strain
𝜙 = {'6':6, '9':9, '12':12, '16':16, '20':20, '25':25, '28':28, '32': 32} #mm
A = {'6':28.27, '9':63.6, '12':113.1, '16':201, '20':314.6, '25':490.8, '28':615.7, '32': 631.3} #mm2


def beta_1():
    if fc <= 30: #N/mm2(MPa)
        β1 = 0.85
    elif 30 < fc < 55: #N/mm2
        β1 = 0.85 -0.05*(fc-30)/7
    else: 
        β1 = 0.65
    return β1

#Safety factor, 𝜙c
def 𝜙x(c,d):
    𝜙x = 0.65+0.25*((1/c/d) - 5/3) #tie
#     𝜙x = 0.75+0.15*((1/c/d) - 5/3) #spiral
    return 𝜙x

#--------------------------------------------------------------
#Try Reinforcement
def try_reinf():
    # 𝜙1  = input("Try tensile steel: diameter in mm :  ")
    # n1 = int(input("Quantity N = :  "))
    while True:
        𝜙1 = input('Try tensile steel: diameter in mm :  ')
        if (𝜙.get(𝜙1 ) == None):
            print("Wrong diameter! select again")
        else:
            n1= int(input('Quantities N = ? : '))           
            break

    while True:
        𝜙2 = input('Try compression steel: diameter in mm :  ')
        if (𝜙.get(𝜙1 ) == None):
            print("Wrong diameter! select again")
        else:
            n2= int(input('Quantities N = ? : '))           
            break
    
    while True:
        𝜙_traverse = input('Try traverse steel: diameter in mm :  ')
        if (𝜙.get(𝜙1 ) == None):
            print("Wrong diameter! select again")
        else:          
            break
              
    As1 = n1*A[𝜙1] #mm2
    As2 = n2*A[𝜙2] #mm2            
    Ag = b*h #mm2
    Ast = As1 + As2 #mm2
    
    print(f"As_tension = {As1:.2f}, As_compression = {As2:.2f}, As_total = {Ast:.2f} : mm2")
    print(f"Ag = {Ag:.2f} mm2")
    return 𝜙1, 𝜙2, 𝜙_traverse, As1, As2, Ast, Ag

#Percent Reinforcement
def percent_reinf(Ast, Ag):    
    ρg = Ast/Ag
    
    if 0.01 < ρg < 0.08:
        print(f"ρg = 0.01 < {ρg:.4f} < 0.08  OK ")
        return ρg
    else: 
        print(f"ρg = {ρg:.4f} out of range [0.01, 0.08]--> Used 0.01")
        ρg = 0.01
        return ρg
    
#Effective depth
def eff_depth(𝜙1, 𝜙2, 𝜙_traverse, h):    
    d2 = cv + 𝜙_traverse + 𝜙2/2 #Effective depth of Compression Steel, mm
    d = h - cv - 𝜙_traverse- 𝜙1/2 #Effective depth of Tension Steel, mm
    print (f"d = {d} mm, d' = {d2} mm")
    return d, d2

#--------------------------------------------------------------
#Forces in the concrete and steel
def force(a, d, d2, As1, As2, fs1, fs2, 𝜙c):
    '''
    Input unit
    a, b, h --> mm
    As --> mm2
    fc, fs --> MPa(N/mm2)
    '''
    Cc = 0.85*fc*a*b*1e-3 #kN, ACI318-14(22.2.2.4.1)
    if fs2 < 0:#Neutal axis < d2 
        Cs = As2*fs2*1e-3 #kN
        
    Cs = As2*(fs2-0.85*fc)*1e-3 #kN
    Ts = As1*fs1*1e-3 #kN
    Pn = Cc+Cs-Ts #kN
    𝜙Pn = 𝜙c*Pn #kN
    Mn = (Cs*(0.5*h-a/2) + Cs*(0.5*h-d2) + Ts*(d-h/2))*1e-3 #kN-m
    𝜙Mn = 𝜙c*Mn #kN-m

    e = Mn*10e3/Pn #mm  
    return Pn, 𝜙Pn, Mn, 𝜙Mn

##Pure Compression
def pure_compression(Ast, Ag, Pn, 𝜙Pn, Pu, p, table):
    #Nominal axial compressive strength at zero eccentricity
    P0 = (0.85*fc*(Ag-Ast)+fy*Ast)*1e-3 #kN, ACI 318-14 (22.4.2.2)

    #Factored axial compressive strength at zero eccentricity
    𝜙P0 = 0.65*P0 #ACI 318-14 (Table 21.2.2)
    𝜙Pn_max = 0.80*𝜙P0 #ACI318-14(Table22.4.2.1)
        
    p.append([0, P0, 0, 𝜙P0]) #save for plotting
    Pn.append(P0)
    𝜙Pn.append(𝜙P0)
    table.append(['Pure Compression', P0, 𝜙P0, 0, 0]) #save to render a table
    return 𝜙Pn_max

##Zero strain at concrete edge c = h
def zero_1(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    c = h #mm
    a = β1*c #mm
    #stress in T-steel
    εs1 = εc*(c-d)/c
    fs1 = -min(fy, Es*εs1) #N/mm2
    #stress in C-steel
    εs2 = εc*(c-d2)/c
    fs2 = min(fy, Es*εs2) #N/mm2
    
    𝜙c = 0.65
    Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2, As1, As2, fs1, fs2, 𝜙c)
    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx) 
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#save for plot
    table.append(['Zero strain at concrete edge ', Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])#save to render a table

## Zero strain at extreme steel reinforcement c = d
def zero_2(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    c = d #mm ACI 318-14 (22.2.2.4.2)
    a = β1*c #mm ACI 318-14 (22.2.2.4.1)
    #stress in T-steel <= 0 
    εs1 = 0
    fs1 = 0 #N/mm2
    #stress in C-steel
    εs2 = εc*(c-d2)/c
    fs2 = min(fy, Es*εs2) #N/mm2

    𝜙c = 0.65 #ACI 318-14 (Table 21.2.2)
    Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2,  As1, As2, fs1, fs2, 𝜙c)
    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx) 
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#save for plot
    table.append(['Zero strain at extreme steel reinforcement ',  Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])#save to render a table

##0.002 strain at extreme steel reinforcement
def t_1(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    εs1 = 0.002  
    c = d*εc/(εs1+εc) #mm
    a = β1*c #mm
    #stress in T-steel
    fs1 = min(fy, Es*εs1) #N/mm2
    #stress in C-steel
    εs2 = εc*(c-d2)/c
    fs2 = min(fy, Es*εs2) #N/mm2
    
    𝜙c = 0.65
    Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2,  As1, As2, fs1, fs2, 𝜙c)
    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx) 
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#save for plot
    table.append(['0.002 strain at extreme steel reinforcement ',  Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])#save to render a table

##Balance
def balance(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    εs1 = fy/Es
    fs1 = fy

    cb = d*εc/(εs1+εc) #mm
    a = β1*cb #mm

    εs2 = (cb-d2)*εc/cb
    fs2 = min(fy, Es*εs2)
    
    𝜙c = 𝜙x(cb, d)
    Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2,  As1, As2, fs1, fs2, 𝜙c)
    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx) 
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#save for plot
    table.append(['Balance',  Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])#save to render a table

## 0.005 strain at extreme steel reinforcement
def t_2(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    εs1 = 0.005 
    c = d*εc/(εs1+εc) #mm
    a = β1*c #mm
    #stress in T-steel     
    fs1 = min(fy, Es*εs1) #N/mm2
    #stress in C-steel
    εs2 = εc*(c-d2)/c
    fs2 = min(fy, Es*εs2) #N/mm2
    𝜙c = 𝜙x(c,d)
    
    Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2,  As1, As2, fs1, fs2, 𝜙c)
    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx) 
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#save for plotting
    table.append(['0.005 strain at extreme steel reinforcement', Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])#save to render a table 

## Pure Bending
def pure_bending(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    εsy = fy/Es
    #Try c
    c = d2 #mm
    while True:       
        εs1 = (d-c)*εc/c
        a = β1*c #mm
        fs1 = fy

        εs2 = (c-d2)*εc/c
        fs2 = min(fy, Es*εs2)
        𝜙b = 0.9
        temp = [0]
        Pnx, 𝜙Pnx, Mnx, 𝜙Mnx = force(a, d, d2,  As1, As2, fs1, fs2, 𝜙b)
        # print(c)
        # print(Pnx)
        if Pnx <= 0 :
            break  
        elif Pnx > temp[-1]: 
            c-=1 #Neutal axis < d2            
        else:
            c+=1 #Neutal axis > d2
        temp.append(Pnx)

    #clean axial   
    Pnx, 𝜙Pnx = 0, 0  

    Pn.append(Pnx)
    𝜙Pn.append(𝜙Pnx)
    Mn.append(Mnx)
    𝜙Mn.append(𝜙Mnx)
    p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])#append in list for plotting
    table.append(['Pure Bending',Pnx, 𝜙Pnx, Mnx,  𝜙Mnx])

## Pure Tension
def pure_tension(As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    Pnt = -fy*(As1+As2)*1e-3
    𝜙Pnt = 0.9*Pnt
    Mnt = 0
    𝜙Mnt = 0

    Pn.append(Pnt)
    𝜙Pn.append(𝜙Pnt)
    Mn.append(Mnt)
    𝜙Mn.append(𝜙Mnt)
    p.append([Mnt, Pnt, 𝜙Mnt, 𝜙Pnt])#save for plotting
    table.append(['Pure Tension', Pnt, 𝜙Pnt, Mnt, 𝜙Mnt])#save to render a table 

#--------------------------------------------------------------
def slender_check(L):
        '''
        Non-swayed frame --> K <= 1
        KL/r < 34 - 12M1b/M2b --> Short column

        Swayed frame --> K > 1
        KL/r < 22 --> Short Column

        rx = sqr(Ix/A), ry = sqr(Iy/A)
        '''
        I = b*(h**3)/12
        r = (I/(b*h))
        if frame == 'non-sway':
            r = 0.3*h
            if K*L/r < 22:
                print("OK. SHORT COLUMN")
                return "1"
            else:
                print(f"ratio = {K*L/r:.2f} > 22")
                return "0"
        #You can comment out for sway case and define M1b/M2b
        #TODO sway cal
        # else:
        #     if (K*L/r < 34 - 12*M1b/M2b):
        #         print("OK. SHORT COLUMN")
        #         return "1"
        #     else:
        #         print(f"ratio = {K*L/r:.2f} > 22")
        #         return "0"

##Calculate IR-diagram coordinate
def coor(β1, d, d2, As1, As2, Ast, Ag,Pn, 𝜙Pn, Mn, 𝜙Mn, p, table):
    '''
    I skip 2 case and skip label of plot for them too.
    You can comment out if you want and provide label at methot plot() too.
    '''
    print(f"\nCalculation coordinate of each case : ")
    𝜙Pn_max = pure_compression(Ast, Ag, Pn, 𝜙Pn, Pu, p, table)#Pure Compression
    zero_1(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#Zero strain at concrete edge c = h
    zero_2(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#Zero strain at extreme steel reinforcement c = d
    #t_1(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#0.002 strain at extreme steel reinforcement --> skip label at plot() too!
    balance(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#Balance
    #t_2(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#0.005 strain at extreme steel reinforcement --> skip label@ at lot() too!
    pure_bending(β1, d, d2, As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#Pure Bending
    pure_tension(As1, As2, Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)#Pure Tension

    return p, Pn, 𝜙Pn, Mn, 𝜙Mn, 𝜙Pn_max

##PLOT
def plot(Ag, p, Pn, 𝜙Pn, Mn, 𝜙Mn, 𝜙Pn_max, Pu, Mu):
    #label of plot point(skip 2 case)
    lb = ['Pure C', 'ε0', 'Zero T', 'Balance', 'Pure M', 'Pure T'] #skip 0.002ε, '0.005ε

    #If you want all, comment out below and comment out at method coor() too.
    # lb = ['Pure C', 'ε0', 'Zero T',  '0.002ε','Balance', '0.005ε', 'Pure M', 'Pure T']

    for i in range(1, len(p)):
        #plot Pn, Mn
        plt.plot(p[i][0], p[i][1], 'bo')
        #plot 𝜙Pn, 𝜙Mn,
        plt.plot(p[i][2], p[i][3], 'go')      
        #plot label
        if lb[i]:
            plt.text(p[i][0], p[i][1], lb[i]+'(' +str(round(p[i][0],2)) +','+' '+ str(round(p[i][1],2))+')')
            plt.text(p[i][2], p[i][3], lb[i]+'(' +str(round(p[i][2],2)) +','+' '+ str(round(p[i][3],2))+')')

    #plot cases
    plt.plot(Mn, Pn)
    plt.plot(𝜙Mn, 𝜙Pn)

    #plot external load
    plt.plot(Mu, Pu, 'rx')
    plt.text(Mu, Pu, (Mu, Pu))

    #plot horizontal line 0f 80%Pure compression
    plt.axhline(y=𝜙Pn_max, color='r', linestyle='-')

    #plot horizontal line 0f 10%Pure compression
    plt.axhline(y=0.65*0.10*fc*Ag/1000, color='g', linestyle='-')    
    plt.title('M-P Diagram', fontsize = 12)
    plt.xlabel("Mn, 𝜙Mn: kN-m") 
    plt.ylabel("Pn, 𝜙Pn: kN")
    plt.axhline(linewidth = 3)

#--------------------------------------------------------------
###DESIGN###
def main():
    print("COLUMN DESIGN : USD METHOD")
    print("Credit Article: https://structurepoint.org/pdfs/Interaction-Diagram-Tied-Reinforced-Concrete-Column-Symmetrical-ACI318-14.htm")
    print("#========================================================================================================")

    while True:
        ##from user input
        print("MATERIAL PROPERTIES")      
        print(f"f'c = {fc} MPa, fy ={fy} MPa, fv = {fv} MPa, Es = {Es} MPa")
        β1 = beta_1()

        ##from user input
        print(f"\nLOAD :")      
        print(f"Pu = {Pu} kN, Mux = {Mux} kN-m, Muy = {Muy} kN-m")

        ##from user input
        global b, h, L
        print(f"\nGEOMTRY :")     
        print(f"Try Geometry : b = {b} mm, h = {h} mm, L = {L} mm")
        check = slender_check(L)
        if check == "0":
            ask = input("Column is slender. Are you continue? y/n : ")
            if ask == "n":
                print("Consider bracing or edit section then run again.")
                break

        ##from user input
        #Try reinf.
        print(f"\nTry Reinforcement :")
        𝜙1, 𝜙2, 𝜙_traverse, As1, As2, Ast, Ag = try_reinf()
        s = min(16*int(𝜙1), 48*int(𝜙_traverse), min(b, h), 300)
        print(f"traverse spacing = {s} mm")

        ##Design
        ##X_X axis
        print(f"\nX-X axis------")
        d, d2 = eff_depth(int(𝜙1), int(𝜙2), int(𝜙_traverse), h)
        ρg = percent_reinf(Ast, Ag)

        #place holder
        Pn = []
        𝜙Pn = []
        Mn = [0]#for Pure Compression Mn = 0
        𝜙Mn = [0]#for Pure Compression 𝜙MnMn = 0
        p = []#place holder --> p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])
        table = []

        #calculate coordinate
        px, Pnx, 𝜙Pnx, Mnx, 𝜙Mnx, 𝜙Pn_maxx = coor(β1, d, d2, As1, As2, Ast, Ag,Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)
        print(tabulate(table, headers=['App', 'Pn, kN', '𝜙Pn, kN', 'Mn, kN-m', '𝜙Mn, kN-m'], tablefmt="psql"))

        #--------------------------------------------------------------
        ##Y-Y axis
        print(f"\nY-Y axis------")
        #swap 
        b, h = h, b

        check = slender_check(L)
        if check == "0":
            ask = input("Column is slender. Are you continue? y/n : ")
            if ask == "n":
                print("Consider bracing or edit section then run again.")
                break

        d, d2 = eff_depth(int(𝜙1), int(𝜙2), int(𝜙_traverse), h)
        # ρg = percent_reinf()

        #place holder list
        Pn = []
        𝜙Pn = []
        Mn = [0]#for Pure Compression Mn = 0
        𝜙Mn = [0]#for Pure Compression 𝜙MnMn = 0
        p = [] #place holder --> p.append([Mnx, Pnx, 𝜙Mnx, 𝜙Pnx])
        table = []

        #calculate coordinate
        py, Pny, 𝜙Pny, Mny, 𝜙Mny, 𝜙Pn_maxy = coor(β1, d, d2, As1, As2, Ast, Ag,Pn, 𝜙Pn, Mn, 𝜙Mn, p, table)
        print(tabulate(table, headers=['App', 'Pn, kN', '𝜙Pn, kN', 'Mn, kN-m', '𝜙Mn, kN-m'], tablefmt="psql"))

        #--------------------------------------------------------------
        #plot IR diagram --> Matplotlib
        plt.figure(1)

        #plot x-x axis, Mx
        plt.subplot(1, 2, 1) #1 rows, 2 column, first position
        plt.margins(0)
        plot(Ag, px, Pnx, 𝜙Pnx, Mnx, 𝜙Mnx, 𝜙Pn_maxx, Pu, Mux)

        #plot y-y axis, My
        plt.subplot(1, 2, 2)
        plt.margins(0) 
        plot(Ag, py, Pny, 𝜙Pny, Mny, 𝜙Mny, 𝜙Pn_maxy, Pu, Muy)
        plt.grid()  
        plt.show()

        ask = input('Save image? Y|N').upper()
        if ask == 'Y':
            plt.figure().savefig('IR/test.png')
        
        
        print("#========================================================================================================")
        print("If you have any comment, pls contact me at highwaynumber12@gmail.com")
        break

### USER PLAYGROUD ###
fc = 32 #N/mm2 (MPa)
fy = 500 #N/mm2 SD40
fv = 235 #N/mm2 SR24
Es = 200000 #N/mm2
cv = 45 #covering, mm

Pu = 2464 #kN
Mux = 94 #kN-m
Muy = 50 #kN-m

frame = 'non-sway' #swayed frame or non-sway
K = 1
b = 350 #mm
h = 400 #mm
L = 2500 #mm, clear height
#L = L/2 if we brace column at middle --> L = L/2  or any!!

if __name__ == '__main__':
    main()
