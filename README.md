Structure design use python:

Typical RC beam design, USD method, ACI Thailand code:
    User provide Mu, Vu. 
    Program will check it's singly or double reinforcement and calculated As required and traverse spacing.
    Then user provide As and spacing.
    
Beam analysis (beam_analysis.py)
This program adapted from Prof.Fredy Gabriel Ramírez Villanueva repository : https://github.com/SirPrime/MatrixAnalysis-Beams.git
This program render shear force diagram and bending moment diagram of contineous beam, use stiffness metrix method.
User can directly input in beam_analysis.py or from command promp or terminal by use beam_analysis_input.py(make comment all input in  beam_analysis.py before use  beam_analysis_input.py) 
User defined:
    -Es, I of section
    -Length for each stretch in meters
    -Support type for each node indicated by digit 0:fixed, 1:pin, 2:free end
    -External force metrix for each node: for beam each node have Fy(+down, -up) and M(+counterclockwise) ex.if you have 2 stretch you must define R0 = [F1y, M1y, F2y, M2y, F3y, M3y], if none defined 0
    -List of list of loads in each stretch: There are 3 type of loads
        1.Point load(P), +down, -up in Newton N
        2.Line load(q). +down, -up in Newton N
        3.Moment concentrate(M), +counterclockwise in N-m
        ex.
        f1 = [P1, m1] list of load in stretch 1
        f2 = [q2] list of load in stretch 2
        loads = [
            f1,
            f2
            ]

Column design(column.py)
User try to provided column section, reinforcement and Pu, Mux, Muy
Program will render IR diagram of section provided in x-x axix and  y-y axis
User can see coordinate of (Mux, Pu) and (Muy, Pu) in IR graph