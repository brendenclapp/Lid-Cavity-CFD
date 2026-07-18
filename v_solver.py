
import numpy as np
from scipy.linalg import solve_banded

def Coeff_v(geo, var, Fields, Faces):

    Fields.v_old = np.copy(Fields.v) 

    for i in range (geo.Nx):
            for j in range (1, geo.Ny):
                
                # (1, Ny) no need to generate coefficients for j == 0 or j == Ny because the BC already tell us what we need to know
                # VF_we/ns already take into account the Mesh GUI via its reliance on v field

                # Staggered face values for v nodes formed from averaging the scalar centered faces

                Faces.Fn_v[i,j] = 0.5*(Faces.F_ns[i,j]+Faces.F_ns[i,j+1])
                Faces.Fs_v[i,j] = 0.5*(Faces.F_ns[i,j]+Faces.F_ns[i,j-1])
                Faces.Fe_v[i,j] = 0.5*(Faces.F_we[i,j]+Faces.F_we[i,j-1])
                Faces.Fw_v[i,j] = 0.5*(Faces.F_we[i+1,j]+Faces.F_we[i+1,j-1])

                # TVD neighbor coefficients formed from v CV faces and v CV diffusion terms

                Faces.a_w_v[i,j] = var.Dw + max(Faces.Fw_v[i,j], 0)
                Faces.a_e_v[i,j] = var.De + max(-Faces.Fe_v[i,j], 0)
                Faces.a_s_v[i,j] = var.Dn + max(Faces.Fs_v[i,j], 0)
                Faces.a_n_v[i,j] = var.Ds + max(-Faces.Fn_v[i,j],0)

              
                Faces.a_P_v[i,j] = Faces.a_w_v[i,j] + Faces.a_e_v[i,j] + Faces.a_s_v[i,j] + Faces.a_n_v[i,j] #+ (Faces.Fe_v[i,j]-Faces.Fw_v[i,j]) + (Faces.Fn_v[i,j]-Faces.Fs_v[i,j])

                Faces.dPdy[i,j] = (Fields.P[i,j]-Fields.P[i,j-1]) * geo.dy                     

    print('a_w_v')
    print(np.array2string(np.flipud(Faces.a_w_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_e_v')
    print(np.array2string(np.flipud(Faces.a_e_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_n_v')
    print(np.array2string(np.flipud(Faces.a_n_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_s_v')
    print(np.array2string(np.flipud(Faces.a_s_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_P_v')
    print(np.array2string(np.flipud(Faces.a_P_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))


def TDMA_v(var, geo, Fields, Faces, Tri):
     
    iv = 0
    jv = 0
    ittv = 0

    for iv in range (geo.Nx):
        for jv in range (1, geo.Ny):

            k = jv - 1
             

            Tri.upper_v[k] = -Faces.a_n_v[iv,jv]
            Tri.diag_v[k] = Faces.a_P_v[iv,jv]
            Tri.lower_v[k] = -Faces.a_s_v[iv,jv]
            Tri.RHS_v[k] = (-Faces.dPdy[iv,jv])


            # West Inlet, dirchlett
            if iv == 0:         
                 
                Tri.RHS_v[k] += (Faces.a_w_v[iv,jv]*var.v_inlet) 
                Tri.RHS_v[k] += (Faces.a_e_v[iv,jv]*Fields.v[iv+1,jv])


            # East outlet, Zero-gradiant
            elif iv == geo.Nx-1:
                 
                 Tri.diag_v[k] -= Faces.a_e_v[iv,jv]
                 Tri.RHS_v[k] += (Faces.a_e_v[iv,jv]*Fields.v[iv,jv])

            # Interior
            else:
                 
                Tri.RHS_v[k] += (Faces.a_w_v[iv,jv]*Fields.v[iv-1,jv]) 
                Tri.RHS_v[k] += (Faces.a_e_v[iv,jv]*Fields.v[iv+1,jv])

            # South wall
            if jv == 1:
                 
              
                Tri.lower_v[k] = 0
                Tri.RHS_v[k] += Faces.a_s_v[iv,jv] * 0

            # North wall
            elif jv == geo.Ny-1:
                    

                Tri.upper_v[k] = 0
                Tri.RHS_v[k] += Faces.a_n_v[iv,jv] * 0

        print('-------------------------------------------------------------------------')
        print('column', iv)
        print('lower')
        print(np.array2string(Tri.lower_v, precision=6))
        print('diag')
        print(np.array2string(Tri.diag_v, precision=6))
        print('upper')
        print(np.array2string(Tri.upper_v, precision=6))
        print('RHS')
        print(np.array2string(Tri.RHS_v, precision=6))

        ab = np.zeros((3, geo.Ny-1))
        ab[0, 1:] = Tri.upper_v[:-1]
        ab[1,:] = Tri.diag_v
        ab[2,:-1] = Tri.lower_v[1:]
        sol = solve_banded((1,1), ab, Tri.RHS_v)
        Fields.v[iv,1:geo.Ny] = sol
    

    Fields.v_psu = np.copy(Fields.v) 

    print('v_TDMA results')
    print(np.array2string(np.flipud(Fields.v_psu.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))