import numpy as np
from scipy.linalg import solve_banded

def Coeff_P(geo, var, Faces, Fields, Coupler):

    Coupler.b[:] = 0.0
    Fields.P[0,0] = 0

    for i in range(geo.Nx):
        for j in range(geo.Ny):

            #================= EAST ==========================
            if i == geo.Nx-1:
            
                Faces.a_e_P[i,j] = 0

            else:
                Faces.a_e_P[i,j] = -((var.rho*geo.dy*geo.dy)/ Faces.a_P_u[i+1,j])


            #================== WEST ==================================
            if i == 0:

                Faces.a_w_P[i,j] = 0

            else:
                Faces.a_w_P[i,j] = -((var.rho*geo.dy*geo.dy)/ Faces.a_P_u[i,j])

            #======================== NORTH====================================

            if j == geo.Ny-1:

                Faces.a_n_P[i,j] = 0

            elif Faces.a_P_v[i,j+1] < 1e-12:

                Faces.a_n_P[i,j] = 0

            else:
                Faces.a_n_P[i,j] = -((var.rho*geo.dx*geo.dx)/ Faces.a_P_v[i,j+1])

            #======================= SOUTH ====================================
            if j == 0:
            
                Faces.a_s_P[i,j] = 0

            elif Faces.a_P_v[i,j] < 1e-12:

                Faces.a_s_P[i,j] = 0

            else:

                Faces.a_s_P[i,j] = -((var.rho*geo.dx*geo.dx)/ Faces.a_P_v[i,j])

            # ====================== CENTER ===============================
            Faces.a_P_P[i,j] = (-Faces.a_e_P[i,j]) + (-Faces.a_w_P[i,j]) + (-Faces.a_n_P[i,j]) + (-Faces.a_s_P[i,j])

            #======================== MASS IMBALANCE ===================================
            if i == 0:

                Coupler.b[i,j] += (var.rho*Fields.u[i+1,j])*geo.dy

            elif i == geo.Nx-1:

                Coupler.b[i,j] += -(var.rho*Fields.u[i,j])*geo.dy

            else:

                Coupler.b[i,j] += (((var.rho*Fields.u[i+1,j])-(var.rho*Fields.u[i,j]))*geo.dy)

            if j == 0:

                 Coupler.b[i,j] += var.rho*Fields.v[i,j+1]*geo.dx

            elif j == geo.Ny-1:

                Coupler.b[i,j] += -(var.rho*Fields.v[i,j])*geo.dx

            else:

                Coupler.b[i,j] += ((var.rho*Fields.v[i,j+1])-(var.rho*Fields.v[i,j]))*geo.dx

                 
    print('a_w_P')
    print(np.array2string(np.flipud(Faces.a_w_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_e_P')
    print(np.array2string(np.flipud(Faces.a_e_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_n_P')
    print(np.array2string(np.flipud(Faces.a_n_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_s_P')
    print(np.array2string(np.flipud(Faces.a_s_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_P_P')
    print(np.array2string(np.flipud(Faces.a_P_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))

            

def TDMA_P(geo, var, Faces, Fields, Coupler, Tri):


    #========================================== COLUMNS ============================
    for ittp in range (1):
        
        for i in range (geo.Nx):
            Tri.upper_P[:] = 0
            Tri.diag_P[:]  = 0
            Tri.lower_P[:] = 0
            Tri.RHS_P[:]   = 0
            for j in range (geo.Ny):

                Tri.upper_P[j] = Faces.a_n_P[i,j]
                Tri.diag_P[j]  = Faces.a_P_P[i,j]
                Tri.lower_P[j] = Faces.a_s_P[i,j]
                Tri.RHS_P[j]   = -Coupler.b[i,j]

                if i == 0:
                    
                    Tri.RHS_P[j] += (-Faces.a_e_P[i,j]*Fields.Pp[i+1,j])

                elif i == geo.Nx-1:

                    Tri.RHS_P[j] += (-Faces.a_w_P[i,j]*Fields.Pp[i-1,j])

                else:
                            
                    Tri.RHS_P[j] += (-Faces.a_w_P[i,j]*Fields.Pp[i-1,j] + -Faces.a_e_P[i,j]*Fields.Pp[i+1,j])


                if j == 0:

                    Tri.lower_P[j] = 0
                    
                    #Residual
                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]

                elif j == geo.Ny-1:

                    Tri.upper_P[j] = 0

                    #Residual
                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] -Faces.a_s_P[i,j]*Fields.Pp[i,j-1]

                else:

                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] -Faces.a_s_P[i,j]*Fields.Pp[i,j-1] -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]
            

            ab = np.zeros((3, geo.Ny))
            
            ab[0,1:] = Tri.upper_P[:-1]
            ab[1,:]  = Tri.diag_P
            ab[2,:-1] = Tri.lower_P[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_P)
            Fields.Pp[i,:] = sol

            print('-------------------------------------------------------------------------')
            print('column', i)
            print('lower')
            print(np.array2string(np.flipud(Tri.upper_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('diag')
            print(np.array2string(np.flipud(Tri.diag_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('upper')
            print(np.array2string(np.flipud(Tri.lower_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('RHS')
            print(np.array2string(np.flipud(Tri.RHS_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))

        print('P_TDMA results')
        print(np.array2string(np.flipud(Fields.Pp.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
                
                    

        #========================================== ROWS ============================
        for j in range (geo.Ny):
            Tri.upper_P[:] = 0
            Tri.diag_P[:]  = 0
            Tri.lower_P[:] = 0
            Tri.RHS_P[:]   = 0
            for i in range (geo.Nx):

                Tri.upper_P[i] = Faces.a_e_P[i,j]
                Tri.diag_P[i]  = Faces.a_P_P[i,j]
                Tri.lower_P[i] = Faces.a_w_P[i,j]
                Tri.RHS_P[i]   = -Coupler.b[i,j]

                if i == 0:

                    Tri.lower_P[i] = 0

                    #Residiual
                    Tri.RHS_P[i] += -Faces.a_e_P[i,j]*Fields.Pp[i+1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j]

                elif i == geo.Nx-1:
                
                    Tri.upper_P[i] = 0

                    #Residiual
                    Tri.RHS_P[i] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j]

                else:
                
                    #Residiual
                    Tri.RHS_P[i] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j] -Faces.a_e_P[i,j]*Fields.Pp[i+1,j]

                if j == 0:
                
                    Tri.RHS_P[i] += -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]

                elif j == geo.Ny-1:

                    Tri.RHS_P[i] += -Faces.a_s_P[i,j]*Fields.Pp[i,j-1] 

                else:

                    Tri.RHS_P[i] += -Faces.a_n_P[i,j]*Fields.Pp[i,j+1] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1]


            ab = np.zeros((3, geo.Nx))
            
            ab[0,1:] = Tri.upper_P[:-1]
            ab[1,:]  = Tri.diag_P
            ab[2,:-1] = Tri.lower_P[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_P)
            Fields.Pp[:,j] = sol


            print('-------------------------------------------------------------------------')
            print('column', i)
            print('lower')
            print(np.array2string(np.flipud(Tri.upper_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('diag')
            print(np.array2string(np.flipud(Tri.diag_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('upper')
            print(np.array2string(np.flipud(Tri.lower_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
            print('RHS')
            print(np.array2string(np.flipud(Tri.RHS_P.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))

    
        
    print('P_TDMA results')
    print(np.array2string(np.flipud(Fields.Pp.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))
    
    
        