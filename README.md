Structure design use python:<br />

Typical RC beam design, USD method, ACI Thailand code:<br />
&emsp;     User provide Mu, Vu.<br />
&emsp;     Program will check it's singly or double reinforcement and calculated As required and traverse spacing.<br />
&emsp;     Then user provide As and spacing.<br />
    
Beam analysis (beam_analysis.py)<br />
This program adapted from Prof.Fredy Gabriel Ramírez Villanueva repository : https://github.com/SirPrime/MatrixAnalysis-Beams.git<br />
This program render shear force diagram and bending moment diagram of contineous beam, use stiffness metrix method<br />
User can directly input in beam_analysis.py or from command promp or terminal by use beam_analysis_input.py(make comment all input in beam_analysis.py before use  beam_analysis_input.py)<br />
User defined:<br />
&emsp;     -Es, I of section<br />
&emsp;     -Length for each stretch in meters<br />
&emsp;     -Support type for each node indicated by digit 0:fixed, 1:pin, 2:free end<br />
&emsp;     -External force metrix for each node: for beam each node have Fy(+down, -up) and M(+counterclockwise) ex.if you have 2 stretch you must define R0 = [F1y, M1y, F2y, M2y, F3y, M3y], if none defined 0<br />
&emsp;     -List of list of loads in each stretch: There are 3 type of loads<br />
&emsp;         1.Point load(P), +down, -up in Newton N<br />
&emsp;         2.Line load(q). +down, -up in Newton N<br />
&emsp;         3.Moment concentrate(M), +counterclockwise in N-m<br />
&emsp;         ex.<br />
&emsp;         f1 = [P1, m1] list of load in stretch 1<br />
&emsp;         f2 = [q2] list of load in stretch 2<br />
&emsp;         loads = [f1, f2]<br />

Column design(column.py)<br />
User try to provided column section, reinforcement and Pu, Mux, Muy<br />
Program will render IR diagram of section provided in x-x axis and  y-y axis<br />
User can see coordinate of (Mux, Pu) and (Muy, Pu) in IR graph<br />
