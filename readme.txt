
----------------------------------------------
This is the readme for the model associated with the paper: 
"Bioinspired figure-ground discrimination via visual motion smoothing". Submitted to Plos Computational Biology. 
(Nov 12, 2022).

This matlab code was contributed by Zhihua Wu.

The main program is EMD_Lobula_network.m, which is self explanatory. The main program needs 5 functions as listed below.
1) stimulus.m
2) LobulaUnits.m
3) Rect.m
4) RK4_Lc.m
5) F_measure.m

When EMD_Lobula_network is typed into the matlab command line, you should see five figures. Changing the values of the variables that define the stimulus and the lobula units (e.g., motion_type, v_b, v_o, dot_size, Lc_tau, and receptive_s in EMD_Lobula_network.m) (variables controlling the stimulus contrast are in stimulus.m), you should obtain results shown in figures in the paper. 
=============================================