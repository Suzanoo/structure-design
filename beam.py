### Typical RC beam design : USD method ###
import numpy as np
import pandas as pd
from tabulate import tabulate

#=============================================================================================
###PREPARATION
𝜙 = {'6':6, '9':9, '12':12, '16':16, '20':20, '25':25, '28':28, '32': 32} #mm
A = {'6':0.2827, '9':0.636, '12':1.131, '16':2.01, '20':3.146, '25':4.908, '28':6.157, '32': 6.313} #cm2

def beta(fc): #β1
    if fc <= 30: #MPa
        β1 = 0.85
    elif fc > 30 and fc < 55: #MPa
        β1 = 0.85 -0.05*(fc-30)/7
    else: 
        β1 = 0.65
    return β1

#Effective depth
def eff_depth(dia_main, dia_traverse, c, h):    
    d1 = c + dia_traverse/10 + dia_main/10/2 #Effective depth of Compression Steel
    d = h - d1 #Effective depth of Tension Steel
    print (f"Effective Depth : \nd = {d} cm, d1 = {d1} cm")  
    return d, d1

#percent Reinforcement
def percent_reinf(fc, fy):
    β1 = beta(fc)
    ρmin = max(np.sqrt(fc)/(4*fy), 1.4/fy)
    ρb = (0.85*fc/fy)*β1*(600/(600+fy))
    ρmax1 = 0.75*ρb
    ρ = 0.50*ρb  #conservative!!!
    print(f"\n% Reinforcement : \nρmin = {ρmin:.4f}, ρmax = {ρmax1:.4f}")  
    return ρmin, ρb, ρmax1, ρ

#Section Capacity
def capacity(𝜙b, ρ, b, d, fc, fy):    
    As = ρ*b*d #cm2
    a = As*fy/(0.85*fc*b) #cm.
    𝜙Mn1 = 𝜙b*As*fy*(d-a/2)/1000 #kN-m  
    return 𝜙Mn1

#Classification
def classification(Mu, 𝜙Mn1):
    𝜙Mn2 = Mu-𝜙Mn1 #kN-m  
    print(f"\nClassification")
    if 𝜙Mn2 > 0:
        class_ = 'double_reinforcement'
        print('double_reinforcement')
        print(f'Mu = {Mu:.2f}, 𝜙Mn1 = {𝜙Mn1:.2f}, 𝜙Mn2 = {𝜙Mn2:.2f} :kg.m')
    else:
        class_ = 'singly_reinforcement'
        𝜙Mn2 = 0
        print('singly_reinforcement')
        print(f'Mu = {Mu:.2f}, 𝜙Mn1 = {𝜙Mn1:.2f}, 𝜙Mn2 = {𝜙Mn2:.2f} :kg.m')      
    return class_

###If it's single reinforcement
def singly_reinf(𝜙b, fc, fy, ρmin, ρmax1, b, d, Mu):   
    Ru = abs(Mu)*1000/(b*d**2) #MPa
    ρ_req = 0.85*(fc/fy)*(1-(1-2*(Ru/𝜙b)/(0.85*fc))**0.5)

    if ρ_req > ρmin:
        As_major = ρ_req*b*d #As_minor = 0
    else:
        As_major = ρmin*b*d #As_minor = 0

    print(f'Ru = {Ru:.2f}, ρ_req = {ρ_req:.4f}, ρmin = {ρmin:.4f}, ρmax = {ρmax1:.4f}')
    print(f'As_major = {As_major:.2f} cm2, As_minor = 0 cm2')
    return As_major

###If it's double reinforcement
def double_reinf(𝜙b, fc, fy, ρ, ρmax1, b, d, d1, 𝜙Mn1, Mu): 
    β1 = beta(fc)  
    𝜙Mn2 = abs(Mu)-𝜙Mn1
    As1 = ρ*b*d #cm2
    As2 = (𝜙Mn2*1000/𝜙b) /(fy*(d-d1)) #cm2
    As_major = As1 + As2 #cm2
       
    ρ1 = (As1+As2)/b*d
    ρ2 = As2/b*d
    
    if ρ1-ρ2 > 0.85*fc*d1*β1*(600/(600-fy))/(fy*d):
        fs = fy
        As_minor = As2
        print(f"fs' = {fs:.2f} --> yeild : OK")
        
        #check ρ
        if ρmax1+ρ2 > 1.4/fy and ρmax1+ρ2 < ρ1:
            print("ρmin < ρ < ρmax ---> OK")
        else:
            print("ρ ---< Out of range")
    else:
        fs = 600*(1-(d1/d)*(600+fy)/600)
        print(f"fs' = {fs:.2f} --> fs' not yeild ")
        
        a = β1*d*(600/(600-fy))
        As_minor = (As1*fy-0.85*fc*a*b)/fs
    
    print(f'As_major = {As_major:.2f} cm2, As_minor = {As_minor:.2f} cm2')
    return fs, As_major, As_minor   

#=============================================================================================
###METHOD
#design stirrup
def shear_design(𝜙v, fc, fv, b, d, Av, Vu):    
    𝜙Vc = 𝜙v*np.sqrt(fc)*b*d/60 #kN --> f'c=MPa(N/mm2),
    𝜙Vs = np.abs(Vu - 𝜙Vc )#kN
    s_req = 𝜙v*Av*fv*d/10*𝜙Vs #cm
    
    lg = (1/3)*𝜙v*np.sqrt(fc)*b*d/10 #light shear, kN
    hv = (2/3)*𝜙v*np.sqrt(fc)*b*d/10 #heavy shear, kN
         
    if Vu <= 𝜙Vc:
        s_max = min(3*Av*fv/b, d/2, 60)#cm, ACI 11.5.5.3;11-13
        print(f'Vu = {Vu:.2f}, 𝜙Vc = {𝜙Vc:.2f}, 𝜙Vs = 0')
        print(f's_max = {s_max:.2f} cm')
        return s_req, s_max

    elif 𝜙Vc < Vu <= 𝜙Vc+ lg:
        s_max = min(𝜙v*Av*fv*d/(𝜙Vs*10), d/2, 60) #cm, ACI 11.5.6.4;11-16
        print(f'Vu = {Vu:.2f}, 𝜙Vc = {𝜙Vc:.2f}, 𝜙Vs = {𝜙Vs:.2f}')
        print(f's_max = {s_max} cm')
        return s_req, s_max

    elif 𝜙Vc + lg  < Vu <= 𝜙Vc+ hv:
        s_max = min(𝜙v*Av*fv*d/(𝜙Vs*10), d/4, 30) #cm, ACI 11.5.6.4;11-16
        print(f'Vu = {Vu:.2f}, 𝜙Vc = {𝜙Vc:.2f}, 𝜙Vs = {𝜙Vs:.2f}')
        print(f's_max = {s_max:.2f} cm')
        return s_req, s_max
        
    else:
        print(f'Vu = {Vu:.2f}, 𝜙Vc = {𝜙Vc:.2f}, 𝜙Vs = {𝜙Vs:.2f}')
        print('Heavy shear --> Revised cross section') 

#design reinforcment
def mainbar_req(𝜙b, fc, fy, c, b, h, Mu, dia_main, dia_traverse, d, d1):   
    
    #percent reinforcement
    ρmin, ρb, ρmax1, ρ = percent_reinf(fc, fy)
    
    #checked capacity of section
    𝜙Mn1 = capacity(𝜙b, ρ, b, d, fc, fy)

    #classification
    class_ = classification(Mu, 𝜙Mn1)

    if class_ == "singly_reinforcement":
        As_major = singly_reinf(𝜙b, fc, fy, ρmin, ρmax1, b, d, Mu)
        data = [As_major]
        return data
    else:
        fs, As_major, As_minor = double_reinf(𝜙b, fc, fy, ρ, ρmax1, b, d, d1, 𝜙Mn1, Mu)
        data = [fs, As_major, As_minor]
        return data

def traverse(𝜙v, fc, fv, b, d, Av, Vu): 
    # d, d1 = eff_depth(𝜙1, 𝜙2, c, h)
    #Calculation traverse spacing
    s_req, s_max = shear_design(𝜙v, fc, fv, b, d, Av, Vu)
    print(f'Traverse design : s_req = {s_req:.2f} cm, s_max = {s_max:.2f} cm')
    input('Select spacing : ')

#=============================================================================================
###DESIGN
def design(fy, data):          
    while True:
        # Double Reinforcement
        # Decision for compression steel if it's not yeild
        if (len(data) != 1) and (data[0] < fy): # data = [[fs, As_major, As_minor]]
            ask = input("fs' not yeild --> Break : Y/N : ").upper()
            if ask == 'Y':
                break
            else:
                pass

            t = ['As_major', 'As_minor']
            for i in range(1, len(data)-1): # data = [[fs, As_major, As_minor]]
                print(f'{t[i-1]} reinf : As-req = {data[i]:.2f} cm2 : Please provide As : ')
                while True:
                    dia = input('Diameter = ? : ')
                    if (𝜙.get(dia) == None):
                        print("Wrong diameter! select again")
                    else:
                        N = int(input('Quantities N = ? : '))
                        As_assign = N*A[dia]
                        if (As_assign < data[i]):
                            print("As_provide < As_required --> select again")               
                        else:
                            break

                if dia in (6, 9 ):
                    result = (f'{N} - RB{dia} : As = {As_assign:.2f} cm2')
                    print(f'{result} --> OK ')
                else:
                    result = (f'{N} - DB{dia} : As = {As_assign:.2f} cm2')
                    print(f'{result} --> OK ')

        ##Singly reinforcement
        else:
            print(f'As-req = {data[0]:.2f} cm2 : Please provide As :') #data = [As_major]
            while True:
                dia = input('Diameter = ? : ')
                if (𝜙.get(dia) == None):
                    print("Wrong diameter! select again")
                else:
                    N = int(input('Quantities N = ? : '))
                    As_assign = N*A[dia]
                    if (As_assign < data[0]):
                        print("As_provide < As_required --> select again")               
                    else:
                        break

            if dia in (6, 9 ):
                result = (f'{N} - RB{dia} : As = {As_assign:.2f} cm2')
                print(f'{result} --> OK ')
            else:
                result = (f'{N} - DB{dia} : As = {As_assign:.2f} cm2')
                print(f'{result} --> OK ') 
        break
#=============================================================================================
###Design and make Report
def main(𝜙b, 𝜙v, fc, fy, fv, Es, b, h, l, c, 𝜙1, 𝜙2, Mu1, Mu2, Vu, df):
    print("TYPICAL RC BEAM DESIGN : USD METHOD")
    print("========================================================================================================")
    print("PROPERTIES")
    print(f"f'c = {fc} Mpa, fy = {fy} Mpa, fv = {fv} MPa, Es = {Es} MPa")
    β1 = beta(fc)
    print(f"𝜙b = {𝜙b}, 𝜙v = {𝜙v}, β1 = {β1}")

    print(f"\nGEOMETRY")
    print(f"b = {b} cm, h = {h} cm, l = {l} m")

    print(f"\nLOAD")
    print(f"Moment at midspan Mu1 = {Mu1} kN-m, \nMoment at support Mu2 = {Mu2} kN-m, \nMax shear Vu = {Vu} kN")

    print(f"\nDATABASE")
    print(tabulate(
        df,
        headers=df.columns,
        floatfmt=".2f",
        showindex=True,
        tablefmt="psql",
        )) 
    #effective depth
    d, d1 = eff_depth(𝜙1, 𝜙2, c, h)

    #design at mid-span 
    print(f"AT MID-SPAN")
    data = mainbar_req(𝜙b, fc, fy, c, b, h, Mu1, 𝜙1, 𝜙2, d, d1)
    design(data)

    print(f"\nAT SUPPORT")
    data = mainbar_req(𝜙b, fc, fy, c, b, h, Mu2, 𝜙1, 𝜙2, d, d1)
    design(data)

    print(f"\nTRAVERSE")
    traverse()

    print("========================================================================================================")
### Comment out below for test
'''
### PLAYGROUD ###
#Factor and Constants:
##User input
𝜙b = 0.90
𝜙v = 0.85

#Material:
##User input
fc = 18 #Mpa
fy = 295 #MPa SD30
fv = 235 #Mpa SR24
Es = 200000 #Mpa
c = 3 #covering(cm)

#Database
df = pd.read_csv('Data/Deform_Bar.csv')

#Cross section of beam and beam length
##User input
b, h, l = 20, 40, 4.05 , #cm cm m

#define initial rebar
##User input
𝜙1 = '12' #main
𝜙2 = '6' #traverse
𝜙1, As = 𝜙[𝜙1], A[𝜙1] #main & main area
𝜙2, Av = 𝜙[𝜙2], 2*A[𝜙2] #traverse & traverse area

#Moment at midspan, moment at support, shear
##User input
Mu1, Mu2, Vu = 12, 23.5, 28 #kN-m, kN-m, kN

if __name__ == '__main__':
    main()
'''