##simulation of membrane tube formation for spontaneous and deviatoric curvature
0. Copy all the .m file in a matlab accessible folder
1. Give input in the file of new_sim_oct.m
2. change the directory in line 121 and 132 of tubedBAR_solver_new_par.m
3. do not change anything in the function of endoinit.m and plotMemProfileArea.m
4. Reading results file
	a. shapeD_0aC_0.txt  (D_0 and C_0 are the value of deviatoric and spontaneous curvature)
		stores data columnwise for r, z, C, W, H, lam at each location in the domain of (-A,A), where A is the total area
	b. energyD_0aC_0.txt  (D_0 and C_0 are the value of deviatoric and spontaneous curvature)
		store for lambda_0,k_0, bending energy density, total energy density, length of the tube

