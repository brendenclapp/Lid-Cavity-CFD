import time
to = time.time()

import numpy as np 
import Data ; import u_solver ; import v_solver
import Pressure ; import correction ; import plotting

geo = Data.Geometry()
var = Data.Variables(geo)
Fields = Data.Fields(geo)
Faces = Data.Faces(geo)
Tri = Data.TDMA(geo)


print("Starting SIMPLE Algorithm Solver...")
print(f"Grid: {geo.Nx}x{geo.Ny} | Re: {var.rho * var.u_inlet * geo.ly / var.mu:.1f}")

var.itt = 0

#BCs
Fields.u[0,:] = var.u_inlet
Fields.v[:,0] = 0; Fields.v[:geo.Ny-1] = 0

while ((var.res_u > var.utol) or (var.res_v > var.vtol) or (var.res_P > var.ptol)) and var.itt < 80:
    #Step 1 Momentum coeff
    u_solver.Coeff_u(geo,var,Fields,Faces)
    v_solver.Coeff_v(geo, var, Fields, Faces) 

    print(f"-------{var.itt}----------")
    print(f"resu:{var.res_u}")
    print(f"resv:{var.res_v}")
    print(f"resP:{var.res_P}")

    #Step 2 Momentum solve
    u_solver.TDMA_u(var, geo, Fields, Faces, Tri)
    v_solver.TDMA_v(geo, Fields, Faces, Tri)

    #Step 3 Pressure Coupling
    Pressure.Coeff_P(geo, var, Faces, Fields)
    Pressure.TDMA_P(geo, Faces, Fields, Tri)

    #Step 4 Correction
    correction.correction(geo, Faces, Fields)


    if var.itt % 5 == 0: print(f"Itteration {var.itt} Completed")
    var.itt += 1

#Step 5 Plotting
tf = time.time()
print(f"Simulation Finished at {tf-to:.2f}s with {var.itt} Itterations")
plotting.P_plot(geo,Fields)
