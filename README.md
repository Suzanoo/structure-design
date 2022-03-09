Structure design use python:__

Typical RC beam design, USD method, ACI Thailand code:__
    User provide Mu, Vu. __
    Program will check it's singly or double reinforcement and calculated As required and traverse spacing.__
    Then user provide As and spacing.__
    
Beam analysis (beam_analysis.py)__
This program adapted from Prof.Fredy Gabriel Ramírez Villanueva repository : https://github.com/SirPrime/MatrixAnalysis-Beams.git__
This program render shear force diagram and bending moment diagram of contineous beam, use stiffness metrix method.__
User can directly input in beam_analysis.py or from command promp or terminal by use beam_analysis_input.py(make comment all input in  beam_analysis.py before use  beam_analysis_input.py)__ 
User defined:__
    -Es, I of section__
    -Length for each stretch in meters__
    -Support type for each node indicated by digit 0:fixed, 1:pin, 2:free end__
    -External force metrix for each node: for beam each node have Fy(+down, -up) and M(+counterclockwise) ex.if you have 2 stretch you must define R0 = [F1y, M1y, F2y, M2y, F3y, M3y], if none defined 0__
    -List of list of loads in each stretch: There are 3 type of loads__
        1.Point load(P), +down, -up in Newton N__
        2.Line load(q). +down, -up in Newton N__
        3.Moment concentrate(M), +counterclockwise in N-m__
        ex.__
        f1 = [P1, m1] list of load in stretch 1__
        f2 = [q2] list of load in stretch 2__
        loads = [f1, f2]__

Column design(column.py)__
User try to provided column section, reinforcement and Pu, Mux, Muy__
Program will render IR diagram of section provided in x-x axix and  y-y axis__
User can see coordinate of (Mux, Pu) and (Muy, Pu) in IR graph__
