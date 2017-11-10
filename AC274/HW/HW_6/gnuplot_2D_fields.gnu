set pm3d
unset surface
set pm3d map
set size square
i = 10                                                  		  # change this i value to see different time steps  
splot 'BGK_2.ruv2d' index i u 1:2:3 title 'density'         		  # this is for the DENSITY profile, i corresponds to the various output times
#splot 'BGK_2.ruv2d' index i u 1:2:($4*$4+$5*$5)**0.5 title 'Vel Mag'    #! this is for the VELOCITY MAGNITUDE (de-comment this 
           								 #! line and comment the previous)
