###ISOLATED FOOTING DESIGN : USD METHOD
###Article Credit : https://www.scribd.com/document/400417851/191061 by รศ.อมร พิมานมาศ และคณะ'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objs as go
# %matplotlib inline
from tabulate import tabulate

import beam as beam
##===============================================================================
#CONSTANCE
𝜙 = {'6':6, '9':9, '12':12, '16':16, '20':20, '25':25, '28':28, '32': 32} #mm
A = {'6':0.2827, '9':0.636, '12':1.131, '16':2.01, '20':3.146, '25':4.908, '28':6.157, '32': 6.313} #cm2

##Method
#Cal.critical length
def initial(B, b, hc, d):
    xv = B/2+b+d/100 #Vu critical d from column edge
    x1 = B/2-b/2-d/200 #d/2 from the periphery of the column left, m
    x2 = B/2+b/2+d/200 #d/2 from the periphery of the column right, m
    Ap = (x2-x1)*(hc+d/100) #punching bearing area, m2
    p = 2*(b+d/100)+2*(hc+d/100) #punching perimeter, m
    xm = B/2+b/2 #Mu critical at column edge, m
    return xv, x1, x2, Ap, p, xm

#Cahculate soil capacity response footing area
#case axial load only --> qu = uniform load
#with external moment --> qu = trapazoid load
#SLS qu
def q(d, B, L, P, M, qa, γ):
    q1 = P/(B*L) - M*(B/2)/((1/12)*L*B**3)+γ*B*L*d/100 #kN/m2
    q2 = P/(B*L) + M*(B/2)/((1/12)*L*B**3)+γ*B*L*d/100 #kN/m2
    if q1 < qa and q2 < qa:
        print(f"q1 = {q1:.2f} kN/m2, q2 = {q2:.2f} kN/m2 < q_allow = {qa:.2f} kN/m2---> OK")
    else:
        print(f"q1 = {q1:.2f} kN/m2, q2 = {q2:.2f} kN/m2 > q_allow = {qa:.2f} kN/m2 ---> NOT OK")
        
#ULS qu
def qu(d, B, L, Pu, Mu, γ):
    qu1 = Pu/(B*L) - Mu*(B/2)/((1/12)*L*B**3)+γ*B*L*d/100 #kN/m2
    qu2 = Pu/(B*L) + Mu*(B/2)/((1/12)*L*B**3)+γ*B*L*d/100 #kN/m2
    print(f"qu1 = {qu1:.2f} kN/m2, qu2 = {qu2:.2f} kN/m2")
    return qu1, qu1

#Create qu equation at x point for any critical case   
def ωux(qu1, qu2, B, x):      
    ωux = qu1+ (qu2-qu1)*x/B #linear equation 
    return ωux #kN/m2

#check shear
def shear(x, B, L, d, qu1, qu2,):
    ωu = ωux(qu1, qu2, B, x) #qu at Vu-critical plane, kN/m2   
    Vu = (1/2)*(ωu+qu2)*(B-x)*L #kN    
    𝜙Vn = 𝜙v*(1/6)*np.sqrt(fc)*(B-x)*L*1000 #kN --> f'c = N/mm2, BL = mm2
    
#     vn = (Vu/𝜙v)/(1000*(B-x)*L) #N/mm2
#     v_allow = np.sqrt(fc)/6 
    Vu = np.abs(Vu)
    𝜙Vn = np.abs(𝜙Vn)

    if 𝜙Vn > Vu:
        print(f"𝜙Vn = {𝜙Vn:.2f} N/mm2 > Vu ={Vu:.2f} N/mm2 --> Shear capacity OK")
    else:
        print(f"𝜙Vn = {𝜙Vn:.2f} N/mm2 < Vu ={Vu:.2f} N/mm2 --> Shear capacity NOT OK")

#check punching shear        
def punching(x1, x2, B, L, Ap, p, d, qu1, qu2):
    ''' x1 --> d/2 from the periphery of the column left, m
        x2 --> d/2 from the periphery of the column right, m
        Ap --> punching bearing area, m2
        p --> punching perimeter, m
    '''
    ωu1 = ωux(qu1, qu2, B, x1) #qu at center plane, kN/m2
    ωu2 = ωux(qu1, qu2, B, x2) #qu at Vu'-critical plane, kN/m2
    
    Vup = (1/2)*(qu1+qu2)*B*L - (1/2)*(ωu2+ωu1)*Ap #kN  
    𝜙Vp = 𝜙v*0.25*np.sqrt(fc)*p*d*10 #kN --> f'c = N/mm2, pd = mm2
    
#     vnp = (Vup/𝜙v)/(p*d)/10 #N/mm2 
#     v_allow = min(1+np.sqrt(fc)/(6*B/L), (2+40*d*p/100)*np.sqrt(fc)/12) #N/mm2
    if 𝜙Vp > Vup:
        print(f"𝜙Vp = {𝜙Vp:.2f} N/mm2 > Vup ={Vup:.2f} N/mm2 --> Punching shear capacity OK")
    else:
        print(f"𝜙Vp = {𝜙Vp:.2f} N/mm2 < Vup ={Vup:.2f} N/mm2 --> Punching shear capacity NOT OK")
        
def moment(x, B, L, qu1, qu2 ): #Mu critical at column edge, m
    wu = ωux(qu1, qu2, B, x) #qu at Mu-critical plane, kN/m2
    
    #case axial load only --> qu = uniform load
    if wu == qu1:
        Mu = (wu*(B-x)*L*(B-x)/2)
        print(f"Mu = {Mu:.2f} kN-m")
        return(Mu) 
    
    #case with external moment --> qu = trapazoid load
    else:
        Mu = (wu*(B-x)*L*(B-x)/2) + (1/2)*(wu+qu2)*(B-x)*L*(2/3)*(B-x)
        print(f"Mu = {Mu:.2f} kN-m")
        return(Mu)  

def calculate(d, d1):
    qu1, qu2 = qu(d, B, L, Pu, Mux, γ)
    print('--------------------------------------------------------------')

    xv, x1, x2, Ap, p, xm = initial(B, bc, hc, d)
    
    shear(xv, B, L, d, qu1, qu2)#criticald at d from column edge
    print('--------------------------------------------------------------')

    punching(x1, x2, B, L, Ap, p, d, qu1, qu2) #critical at d/2 from column edge
    print('--------------------------------------------------------------')

    Mu = moment(xm, B, L, qu1, qu2)#critical at column edge, m

    data = beam.mainbar_req(𝜙b, fc, fy, c, B*100, L*100, Mu, 𝜙1, 𝜙2, d, d1) #mainbar_req(𝜙b, fc, fy, c, b, h, Mu, dia_main, dia_traverse)
    beam.design(fy, data) #design(fy, data)

##=====================================================================================================
###Design && make report
def main():
    print("ISOLATED FOOTING DESIGN : USD METHOD")
    print('Article Credit : https://www.scribd.com/document/400417851/191061 by รศ.อมร พิมานมาศ และคณะ')
    print("========================================================================================================")

    global B, L, bc, hc,Mux, Muy
    print(f"DATABASE")
    df = pd.read_csv('Data/Deform_Bar.csv')
    print(tabulate(
        df,
        headers=df.columns,
        floatfmt=".2f",
        showindex=False,
        tablefmt="psql",
        )) 

    #TODO
    #Prepare for drawing
    𝜙𝜙 = []
    n = []

    print(f"Footing dimension: B = {B} m, L = {L} m, depth = {t} m")
    print(f"Pu = {Pu:.2f} kN/m2, Mux = {Mux:.2f} kN-m, Muy = {Muy:.2f} kN-m")

    print(f"\nCALCULATED")
    #Check qu
    d, d1 = beam.eff_depth(𝜙1, 𝜙2, c, t*100) #eff_depth(dia_main, dia_traverse, depth)
    q(d, B, L, P, M, qa, γ)
    #-------------------------------------------------------------------------------------
    print(f"\nX-X AXIS")
    calculate(d, d1)
    # 𝜙𝜙.append(dia)
    # n.append(N)
    #-------------------------------------------------------------------------------------
    print(f"\nY-Y AXIS")
    B, L = L, B #swap variable
    b = B*100  #cm, define for calling method mainbar_req(b, h, Mu, dia_main, dia_traverse) in beam.py
    bc, hc = hc, bc #swap variable
    Mux, Muy = Muy, Mux #swap variable
    calculate(d, d1)
    # 𝜙𝜙.append(dia)
    # n.append(N)
    print("========================================================================================================")
    print("If you have any comment, pls contact me at highwaynumber12@gmail.com")

#-------------------------------------------------------------------------------------
### USER PLAYGROUD ###
#Material:
##User input
fc = 20 #Mpa
fy = 415 #MPa SD30
fv = 235 #Mpa SR24
Es = 200000 #Mpa

#Factor and Constants:
##User input
𝜙b = 0.90
𝜙v = 0.85

# Geometry
##User input
c = 7.5 #covering, cm
B, L, t,  = 1, .5, .3 #Footing dimension, m
b = B*100 #cm, define for calling method mainbar_req(b, h, Mu, dia_main, dia_traverse) in beam.py
h = t*100 #cm, define for calling method mainbar_req(b, h, Mu, dia_main, dia_traverse) in beam.py
bc, hc = 0.15, 0.15 #column cross section, m
qa = 300 #Bearing capacity of soil(kN/m2)
γ = 18 #Unit wt of soil(kN/m3) 

#try initial rebar
##User input
𝜙1 = '9' #main
𝜙2 = '9' #traverse
𝜙1, As = 𝜙[𝜙1], A[𝜙1]
𝜙2, Av = 𝜙[𝜙2], 2*A[𝜙2]

#Load
##User input
Pu = 2 #kN
Mux = 0 #kN-m
Muy = 0 #kN-m

wt = B*L*t*24 # self wt., kN
Pu = Pu+1.4*wt #kN

P = Pu/1.4
M = max(Mux, Muy)/1.4

if __name__ == '__main__':
    main()
