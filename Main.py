import time
to = time.time()

import numpy as np ; import os ; from PIL import Image ; import glob
import Data ; import u_solver ; import v_solver
import Pressure ; import correction ; import plotting

os.makedirs("frames", exist_ok=True)

geo = Data.Geometry()
var = Data.Variables(geo)
Fields = Data.Fields(geo)
Faces = Data.Faces(geo)
Tri = Data.TDMA(geo)


print("Starting SIMPLE Algorithm Solver...")
print(f"Grid: {geo.Nx}x{geo.Ny} | Re: {var.rho * var.u_lid * geo.lx / var.mu:.1f}")


#BCs
Fields.u[0,:] = 0; Fields.u[geo.Nx, :] = 0 
Fields.v[:,0] = 0; Fields.v[:geo.Ny-1] = 0

t = 0
t_final = 3
dt = 0.05
pitt = 1
var.aPo = (var.rho*geo.dy*geo.dx)/dt

while t < t_final:

    print(f"Time step {f"{t:.3f}"}s of {t_final}s")
    Fields.u_old = np.copy(Fields.u)
    Fields.v_old = np.copy(Fields.v)
    Fields.P_old = np.copy(Fields.P)
    var.res_u = 1 ; var.res_v = 1 ; var.res_P = 1
    var.itt = 1

    while ((var.res_u > var.utol) or (var.res_v > var.vtol) or (var.res_P > var.ptol)) and var.itt < 5 :

        #Step 1 Momentum coeff
        u_solver.Coeff_u(geo,var,Fields,Faces)
        v_solver.Coeff_v(geo, var, Fields, Faces) 

        #Step 2 Momentum solve
        u_solver.TDMA_u(var, geo, Fields, Faces, Tri)
        v_solver.TDMA_v(geo, var, Fields, Faces, Tri)

        #Step 3 Pressure Coupling
        Pressure.Coeff_P(geo, var, Faces, Fields)
        Pressure.TDMA_P(geo, Faces, Fields, Tri)

        #Step 4 Correction
        correction.correction(geo, Faces, Fields)

        var.itt += 1

    #Step 5 Plotting
    plotting.P_plot(geo,Fields, pitt,t)

    #Time / itteration step
    t += dt
    pitt += 1


# Completion
tf = time.time()
print(f"Simulation Finished at elapsed time of {tf-to:.2f}s" )


# Gif Maker
frames = []

for file in sorted(glob.glob("frames/*.png")):
    frames.append(Image.open(file))


frames[0].save(
    "cavity.gif",
    save_all=True,
    append_images=frames[1:] + [frames[-1]]*30,
    duration=150,
    loop=5
)