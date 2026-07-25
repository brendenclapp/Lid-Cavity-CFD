
import numpy as np
from scipy.linalg import solve_banded

def Coeff_v(geo, var, Fields, Faces):


    for i in range (geo.Nx):
        for j in range (1, geo.Ny):


            #===================================== EAST ======================================================== 
                # Left Wall
                if i == 0:

                    ce = var.rho*(Fields.u[i,j] + Fields.u[i,j-1])*0.5
                    de = var.mu * geo.dy / geo.dx
                    we = (var.mu * geo.dx)/ 3*geo.dy
                    
                    Faces.a_e_v[i,j] = -((abs(ce)-ce)/2)*geo.dx - de - we

                elif i == geo.Nx-1:

                    Faces.a_e_v[i,j] = 0

                # East A Interior
                else:
                    ce = var.rho*(Fields.u[i,j] + Fields.u[i,j-1])*0.5
                    de = var.mu*geo.dy/geo.dx
        
                    Faces.a_e_v[i,j] = -((abs(ce)-ce)/2)*geo.dy - de
    
            #===================================== WEST ======================================================== 
                # Left Wall
                if i == 0:

                    Faces.a_w_v[i,j] = 0

                elif i == geo.Nx-1:

                    cw = var.rho*(Fields.u[i+1,j] + Fields.u[i+1,j-1])*0.5
                    dw = var.mu * geo.dy / geo.dx
                    ww = (var.mu * geo.dx)/ 3*geo.dy
        
                    Faces.a_w_v[i,j] = -((abs(cw)+cw)/2)*geo.dy - dw - ww


                # West A Interior
                else:
                    cw = var.rho*(Fields.u[i+1,j] + Fields.u[i+1,j-1])*0.5
                    dw = var.mu * geo.dy / geo.dx
        
                    Faces.a_w_v[i,j] = -((abs(cw)+cw)/2)*geo.dy - dw
    
            #===================================== NORTH ======================================================== 
            #North A Interior
                cn = var.rho*(Fields.v[i,j] + Fields.v[i,j+1])*0.5
                dn = var.mu * geo.dx / geo.dy

                Faces.a_n_v[i,j] = -((abs(cn)-cn)/2)*geo.dx - dn

            #===================================== SOUTH ======================================================== 

            #South A Interior
                cs = var.rho*(Fields.v[i,j] + Fields.v[i,j-1])*0.5
                ds = var.mu * geo.dx / geo.dy

                Faces.a_s_v[i,j] = -((abs(cs)+cs)/2)*geo.dx - ds
            
            #===================================== CENTER ======================================================== 
            #Left Wall
                if i == 0:

                    Faces.a_P_v[i,j] = (((abs(cn)-cn)/2)*geo.dx + dn) + (((abs(cs)+cs)/2)*geo.dx + ds) \
                                        + ((((abs(ce)-ce)/2)*geo.dy + de) + (3*var.mu*geo.dx)/geo.dy)
                #Top Wall
                elif j == geo.Ny-1:

                    Faces.a_P_v[i,j] = (((abs(cn)-cn)/2)*geo.dx + dn) + (((abs(cs)+cs)/2)*geo.dx + ds) \
                                        + ((((abs(cw)+cw)/2)*geo.dy + de) + (3*var.mu*geo.dx)/geo.dy)

                # Central A Interior
                else:
                    Faces.a_P_v[i,j] = (((abs(ce)-ce)/2)*geo.dy + de) + (((abs(cw)+cw)/2)*geo.dy + dw) \
                                    + (((abs(cn)-cn)/2)*geo.dx + dn) + ((abs(cs)+cs)/2)*geo.dx + ds

            #===================================== PRESSURE ========================================================

                Faces.dPdy[i,j] = (Fields.P[i,j-1] - Fields.P[i,j]) * geo.dy
                
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