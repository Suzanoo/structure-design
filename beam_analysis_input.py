import beam_analysis as analysis

#https://en.wikipedia.org/wiki/Young%27s_modulus
E = 200e9 # Es, GPa 
# https://optimalbeam.com/section-properties.php 
# Rec 30x40 cm
I = 67500e-8 #m4

##Span of each stretch
spans = []
i = 1 
print('Define span of each stretch.')
while True:   
    try:   
        s = float(input(f'Enter span for stretch {i} in meters : '))
        if s > 0:
            spans.append(s)
            print(f"span = {spans}")
            ask = input("Finish? Y|N : ").upper()      
            if ask == 'Y':
                break
            else:
                i +=1
        else:
            print('Badly input.Try again')
    except Exception as e:
        print('Badly input.Try again')

#--------------------------------------------------------------------
##Support type
s = []
print(f'\nDefine support type for each node. You have {len(spans)+1} nodes')
for i in range(1, len(spans)+2):
    while True:
        try:
            x = int(input(f'Define support type for node {i} --> fixd:0, pin:1, free:2 : '))
            if x in (0, 1, 2):
                s.append(x)
                print(f"support type = {s}") 
                break  
            else:
                print('Badly input.Try again')
        except Exception as e:
            print('Badly input.Try again')

#--------------------------------------------------------------------
##Nodal external loads
'''
Define external loads for each node.
For each node first define is Fy(kN), next is M(kN-m)
Finaly we have [R0] = ['F1y', 'M1', 'F2y', 'M2', 'F3y', 'M3',...]
'''
R0 = []
print(f'\nDefine external loads for each node. You have {len(s)} nodes')

for i in range(1, len(s)+1):
    while True:
        try:    
            f = float(input(f'Define Fy(kN) for node {i} Up-, Down+ : '))     
            m = float(input(f'Define moment(kN-m) for node {i} counterclockwise + : '))
            R0.append(f)
            R0.append(m)
            print(f'R0 = {R0}')          
            break
        except Exception as e:
            print('Badly input.Try again')

#--------------------------------------------------------------------
#Define loads in each stretch : unit in --> Newton, N
'''
q = DistributedLoad (value, start, length), distance between the left end of the span and the start of the load
P = PointLoad(value, position), Load position with respect to the left end of the section
M = MomentConcentrated (value, position),  position of the moment with respect to the left end of the section'
'''
print(f'\nDefine loads in each stretch : unit in --> Newton, Newton-meters')
print(f'You have {len(spans)} stretch')

loads = [[] for i in range(0, len(spans))] #[[], [], [],...]
for i in range(0, len(loads)):
    print(f'Define load for stretch {i+1} :')
    while True:
        try:
            type = input('Choose load type(P, q , M) or other if none : ').lower()
            if type in ('p', 'q', 'm'):
                if type == 'p':
                    value = float(input('Enter point load P(N) , Down+ Up- : '))
                    x = float(input('Enter position x(m) with respect to the left end of the section : '))
                    f = analysis.PointLoad(value, x)
                    loads[i].append(f)
                elif type == 'q':
                    value = float(input('Enter line load value q(N/m) , Down+ Up- : '))
                    start = float(input('Enter start point x(m) distance between the left end of the span and the start of the load : '))
                    length = float(input('Enter length of line load l(m) : '))
                    f = analysis.DistributedLoad(value, start, length)
                    loads[i].append(f)
                else:
                    value = float(input('Enter moment m(N-m) : '))
                    x = float(input('Enter position x(m) relative to the left node of the stretch, counterclockwise + : '))
                    f = analysis.MomentConcentrated(value, x)
                    loads[i].append(f)
            else:
                print(f'None for stretch {i+1}')
                break

            ask = input(f"Finish for stretch {i+1} Y|N : ").upper()      
            if ask == 'Y':
                break
        except:
            print('Badly input.Try again')
    print('#---------------------------------------------')
    
#print load
for i in range(0, len(loads)):
    print(f'Load in stretch {i+1} : ')
    print(*loads[i], sep = "\n")

##Call 
analysis.main(E, I, spans, s, loads, R0)

#--------------------------------------------------------------------
### PLAY GROUND ###
##  Example from URL : https://learnaboutstructures.com/sites/default/files/images/3-Frames/Det-Beam-Example-Moment.png 
##  Recheck :  https://platform.skyciv.com/beam
### comment out all below to run code ###
'''
spans = [10.0, 5.0]
s = [1, 1, 2]
R0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

Load in stretch 1
P1 = PointLoad(20000, 5)
m1 = MomentConcentrated(-30000, 5)

Load in stretch 2
q2 = DistributedLoad(10000, 0, 5)

f1 = [P1, m1]
f2 = [q2]
loads = [
    f1,
    f2
    ]

'''