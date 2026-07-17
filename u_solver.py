    
import numpy as np



def Coeff_u(geo,var,Fields,Faces):
   
   Fields.u_old = np.copy(Fields.u) 

   for i in range (1, geo.Nx):
        for j in range (geo.Ny):
            
            # (1, Nx) no need to generate coefficients for i == 0 or i == Nx because the BC already tell us what we need to know
            # UF_we/ns already take into account the Mesh GUI via its reliance on u field

            # Staggered face values for u nodes formed from averaging the scalar centered faces
        
            Faces.Fw_u[i,j] = 0.5 * (Faces.F_we[i,j] + Faces.F_we[i-1,j])
            Faces.Fe_u[i,j] = 0.5 * (Faces.F_we[i,j] + Faces.F_we[i+1,j])
            Faces.Fs_u[i,j] = 0.5 * (Faces.F_ns[i,j] + Faces.F_ns[i-1,j])
            Faces.Fn_u[i,j] = 0.5 * (Faces.F_ns[i,j+1] + Faces.F_ns[i-1,j+1])

        
            Faces.a_w_u[i,j] = var.Dw + max(Faces.Fw_u[i,j], 0)
            Faces.a_e_u[i,j] = var.De + max(-Faces.Fe_u[i,j], 0)

            #South wall diffusion added
            if j == 0:
                Faces.a_s_u[i,j] = var.D_wall + max(Faces.Fs_u[i,j], 0)
            else:
                Faces.a_s_u[i,j] = var.Ds + max(Faces.Fs_u[i,j], 0)

            #North wall diffusion added
            if j == geo.Ny-1:
                Faces.a_n_u[i,j] = var.D_wall + max(-Faces.Fn_u[i,j],0)
            else:
                Faces.a_n_u[i,j] = var.Dn + max(-Faces.Fn_u[i,j],0)

            Faces.a_P_u[i,j] = Faces.a_w_u[i,j] + Faces.a_e_u[i,j] + Faces.a_s_u[i,j] + Faces.a_n_u[i,j]
            #Faces.a_P_u[i,j] += (Faces.Fe_u[i,j]-Faces.Fw_u[i,j]) + ((Faces.Fn_u[i,j]-Faces.Fs_u[i,j]))

            Faces.dPdx[i,j] =  (Fields.P[i-1,j]-Fields.P[i,j])*geo.dy                


def TDMA_u(var, geo, Fields, Faces, Tri):

    from scipy.linalg import solve_banded
  
    iu = 0
    ju = 0
    ittu = 0

    for ittu in range (5):
        for iu in range (1, geo.Nx):
            for ju in range (geo.Ny):

                Tri.upper[ju] = -Faces.a_n_u[iu,ju]
                Tri.diag[ju] = Faces.a_P_u[iu,ju]
                Tri.lower[ju] = -Faces.a_s_u[iu,ju]
                Tri.RHS[ju] = (Faces.dPdx[iu,ju])

                #West / East, Inlet BC, Outlet BC
                if iu == 1:

                    Tri.RHS[ju] += (Faces.a_w_u[iu,ju]*var.u_inlet) 
                    Tri.RHS[ju] += (Faces.a_e_u[iu,ju]*Fields.u[iu+1,ju])
                
                elif iu == geo.Nx-1:
                    
                    Tri.diag[ju] -= Faces.a_e_u[iu,ju]
                    Tri.RHS[ju] += (Faces.a_w_u[iu,ju]*Fields.u[iu-1,ju]) 

                else:

                    Tri.RHS[ju] += (Faces.a_w_u[iu,ju]*Fields.u[iu-1,ju]) 
                    Tri.RHS[ju] += (Faces.a_e_u[iu,ju]*Fields.u[iu+1,ju])

                

                # South wall BC
                if ju == 0:

                    Tri.lower[ju] = 0
                    Tri.RHS[ju] += Faces.a_s_u[iu,ju] * 0

                # North wall BC
                if ju == geo.Ny-1:

                    Tri.upper[ju] = 0
                    Tri.RHS[ju] += Faces.a_n_u[iu,ju] * 0


            ab = np.zeros((3, geo.Ny))
            ab[0, 1:] = Tri.upper[:-1]
            ab[1,:] = Tri.diag
            ab[2,:-1] = Tri.lower[1:]
            sol = solve_banded((1,1), ab, Tri.RHS)
            Fields.u[iu,:] = sol
                
    
    Fields.u_psu = np.copy(Fields.u) 
    Fields.u_psu[-1, :] = Fields.u_psu[-2, :]

    print('u_TDMA results')
    print(np.array2string(np.flipud(Fields.u_psu.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))