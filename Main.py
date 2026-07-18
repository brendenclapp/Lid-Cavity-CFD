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
print(f"Grid: {geo.Nx}x{geo.Ny} | Re: {var.rho * var.u_inlet * geo.lx / var.mu:.1f}")


#BCs
mit = 0

for mit in range (50):
    Fields.u[0,:] = var.u_inlet
    Fields.v[:, 0] = 0
    Fields.v[:, -1] = 0

    #---------------------- convective flux  -------------------------------------------------------------

    for i in range (geo.Nx+1):                                                                  # generates west/east faces, F = puA
        for j in range (geo.Ny):    
            Faces.F_we[i,j] = var.rho  * Fields.u[i,j] * geo.dy
        
    for i in range (geo.Nx):                                                                    # generates north/south faces, F = pvA
        for j in range (geo.Ny+1):
            Faces.F_ns[i,j] = var.rho * Fields.v[i,j] * geo.dx

    print('F_we final results')
    print(np.array2string(np.flipud(Faces.F_we.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('F_ns final results')
    print(np.array2string(np.flipud(Faces.F_ns.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))



    u_solver.Coeff_u(geo,var,Fields,Faces)
    u_solver.TDMA_u(var, geo, Fields, Faces, Tri)


    #v_solver.Coeff_v(geo, var, Fields, Faces)
    #v_solver.TDMA_v(var, geo, Fields, Faces, Tri)

    Pressure.Coeff_P(geo, var, Faces, Fields, Coupler)
    Pressure.TDMA_P(geo, var, Faces, Fields, Coupler, Tri)

    correction.correction(geo, var, Faces, Fields, Coupler, Tri)