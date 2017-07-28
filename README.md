# RGraph_minimum_dose_route
## Methods of minimizing doses incurred by external exposure while moving in radiation hazardous areas.
<p>
<b>Abstract.</b>
Radiation doses received by workers during their movement within areas contaminated as a result of events and activities, leading to emergency or existing exposure situations, may provide a substantial contribution to total external exposure during remediation work. This paper describes an approach to minimize worker external exposure in these circumstances, based on graph theory. The paper describes several tasks, including: searching for a route with the lowest dose, searching for an optimal bypass with a given set of control points and searching for the optimal road network coverage. Classical graph theory algorithms have been used (Dijkstra's algorithm, Chinese postman problem and travelling salesman problem). Algorithms for solving the above mentioned problems have been developed and were included in the information-analytical system for radiation safety. This software has been applied for optimization of protection during remediation work at the Andreeva Bay site of temporary storage for spent fuel and radioactive waste in the Kola Peninsula, both in the context of existing exposure situations and improving the preparedness for emergency exposure situations.
<p>
<h3> 1. interpolation.R </h3>
Use this code to create a grid from a data file with measurements of ambient dose rate equivalent.
<br>
<h3>2. Route_finder_from_A_to_B.R </h3>
Use this code to find a route with a minimum dose between two points on the map
<br>
<h3>3. Hamiltonian_path_(optimal_bypass_of_selected_points).R</h3>
Use this code to find an optimal bypass of selected points on the map (by default for 7 points)
<br>
<h3>4. Chinese_postman_problem.R</h3>
Use this code to find a solution of a Chinese postman problem. 
The Chinese postman problem is to find a shortest closed path or circuit that visits every edge of a (connected) undirected graph.
By default route starts in vertex 1.
<p><br>
© Chizhov Konstantin. <br>
Contact email: <a href="mailto:kchizhov@fmbcfmba.ru">kchizhov@fmbcfmba.ru</a>

Paper in JRP: http://iopscience.iop.org/article/10.1088/1361-6498/aa7c4f <br>
Published 28 July 2017 • © 2017 IOP Publishing Ltd  <br>
Journal of Radiological Protection, Volume 37, Number 3 <br> <p>
K. Chizhov et al., Methods of minimising doses incurred by external exposure while moving in radiation hazardous areas. Journal of Radiological Protection. 37, 697–714 (2017).
