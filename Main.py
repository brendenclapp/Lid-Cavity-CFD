import numpy as np
import Data
import u_solver
import v_solver
import Pressure
import correction

geo = Data.Geometry()
var = Data.Variables(geo)
Fields = Data.Fields(geo)
Faces = Data.Faces(geo)
Tri = Data.TDMA(geo)
Coupler = Data.Coupler(geo)

print("Starting SIMPLE Algorithm Solver...")
print(f"Grid: {geo.Nx}x{geo.Ny} | Re: {var.rho * var.u_lid * geo.lx / var.mu:.1f}")



#BCs
Fields.u[0,:] = 0; Fields.u[geo.Nx, :] = 0 
Fields.v[:,0] = 0; Fields.v[:geo.Ny-1] = 0
print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))

for mit in range (1):

    #Step 1 Momentum coeff
    u_solver.Coeff_u(geo,var,Fields,Faces)
    v_solver.Coeff_v(geo, var, Fields, Faces) 

    #Step 2 Momentum solve
    u_solver.TDMA_u(var, geo, Fields, Faces, Tri)
    v_solver.TDMA_v(var, geo, Fields, Faces, Tri)

    #Step 3 Pressure Coupling
    #Pressure.Coeff_P(geo, var, Faces, Fields, Coupler)
    #Pressure.TDMA_P(geo, var, Faces, Fields, Coupler, Tri)

    #Step 4 Correction
    #correction.correction(geo, var, Faces, Fields, Coupler, Tri)