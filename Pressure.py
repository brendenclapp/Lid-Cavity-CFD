import numpy as np
from scipy.linalg import solve_banded

def Coeff_P(geo, var, Faces, Fields):

    Faces.b[:] = 0.0
    Fields.P[geo.Nx-1,:] = 0
    #print('auP results')
    #print(np.array2string(np.flipud(Faces.a_P_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))

    #print('avP results')
    #print(np.array2string(np.flipud(Faces.a_P_v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
   # print('a_v_P results')
   # print(np.array2string(np.flipud(Faces.a_P_v.T),formatter={'float_kind': lambda x: f"{x:8.9f}"}, max_line_width= 1000000000))
   # print('u results')
    #print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:8.9f}"}, max_line_width= 1000000000))
                         

    for i in range(geo.Nx-1):
        for j in range(geo.Ny):

            #================= EAST ==========================
            if i == geo.Nx-2:
            
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

            elif Faces.a_P_v[i,j+1] == 0:

                Faces.a_n_P[i,j] = 0

            else:
                Faces.a_n_P[i,j] = -((var.rho*geo.dx*geo.dx)/ Faces.a_P_v[i,j+1])

            #======================= SOUTH ====================================
            if j == 0:
            
                Faces.a_s_P[i,j] = 0

            elif Faces.a_P_v[i,j] == 0:

                Faces.a_s_P[i,j] = 0

            else:

                Faces.a_s_P[i,j] = -((var.rho*geo.dx*geo.dx)/ Faces.a_P_v[i,j])

            # ====================== CENTER ===============================
            Faces.a_P_P[i,j] = (-Faces.a_e_P[i,j]) + (-Faces.a_w_P[i,j]) + (-Faces.a_n_P[i,j]) + (-Faces.a_s_P[i,j])

            #======================== MASS IMBALANCE ===================================
            if i == 0:

                Faces.b[i,j] += (((var.rho*Fields.u[i+1,j])-(var.rho*var.u_inlet))*geo.dy)

            else:

                Faces.b[i,j] += (((var.rho*Fields.u[i+1,j])-(var.rho*Fields.u[i,j]))*geo.dy)

            if j == 0:

                Faces.b[i,j] += var.rho*Fields.v[i,j+1]*geo.dx

            elif j == geo.Ny-1:

                Faces.b[i,j] += -(var.rho*Fields.v[i,j])*geo.dx

            else:

                Faces.b[i,j] += ((var.rho*Fields.v[i,j+1])-(var.rho*Fields.v[i,j]))*geo.dx  

    var.res_P = 0.0
    for i in range (1, geo.Nx-1):
        for j in range (1, geo.Ny-1):

            var.res_P = var.res_P + ((Faces.a_P_P[i,j]*Fields.P[i,j]) + (Faces.a_w_P[i,j]*Fields.P[i-1,j]) \
                    + (Faces.a_e_P[i,j]*Fields.P[i+1,j]) + (Faces.a_s_P[i,j]*Fields.P[i,j-1]) \
                    + (Faces.a_n_P[i,j]*Fields.P[i,j+1]) - Faces.b[i,j])**2

    var.res_P = np.sqrt(var.res_P)       


    """
    print('apw tot results')
    print(np.array2string(np.flipud(Faces.a_w_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('ape tot results')
    print(np.array2string(np.flipud(Faces.a_e_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('aps tot results')
    print(np.array2string(np.flipud(Faces.a_s_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('ans tot results')
    print(np.array2string(np.flipud(Faces.a_n_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('b tot results')
    print(np.array2string(np.flipud(Faces.b.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('app tot results')
    print(np.array2string(np.flipud(Faces.a_P_P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
         """       

def TDMA_P(geo, Faces, Fields, Tri):

    Fields.Pp[:] = 0.0

    for ittp in range (20):
        
        #========================================== COLUMNS ============================
        for i in range (geo.Nx-1):
            Tri.upper_Pcol[:] = 0
            Tri.diag_Pcol[:]  = 0
            Tri.lower_Pcol[:] = 0
            Tri.RHS_Pcol[:]   = 0
            for j in range (geo.Ny):

                Tri.upper_Pcol[j] = Faces.a_n_P[i,j]
                Tri.diag_Pcol[j]  = Faces.a_P_P[i,j]
                Tri.lower_Pcol[j] = Faces.a_s_P[i,j]
                Tri.RHS_Pcol[j]   = -Faces.b[i,j]

                if i == 0:
                    
                    Tri.RHS_Pcol[j] += (-Faces.a_e_P[i,j]*Fields.Pp[i+1,j])

                elif i == geo.Nx-2:

                    Tri.RHS_Pcol[j] += (-Faces.a_w_P[i,j]*Fields.Pp[i-1,j])

                else:
                            
                    Tri.RHS_Pcol[j] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_e_P[i,j]*Fields.Pp[i+1,j]


                if j == 0:

                    Tri.lower_Pcol[j] = 0
                    
                    #Residual
                    Tri.RHS_Pcol[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_n_P[i,j]*Fields.Pp[i,j+1]

                elif j == geo.Ny-1:

                    Tri.upper_Pcol[j] = 0

                    #Residual
                    Tri.RHS_Pcol[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1]

                else:

                    Tri.RHS_Pcol[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1] -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]
            

            ab = np.zeros((3, geo.Ny))
            
            ab[0,1:] = Tri.upper_Pcol[:-1]
            ab[1,:]  = Tri.diag_Pcol
            ab[2,:-1] = Tri.lower_Pcol[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_Pcol)
            Fields.Pp[i,:] += sol


        #print('p col results')
        #print(np.array2string(np.flipud(Fields.Pp.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))


        #========================================== ROWS ============================
        for j in range (geo.Ny):
            Tri.upper_Prow[:] = 0
            Tri.diag_Prow[:]  = 0
            Tri.lower_Prow[:] = 0
            Tri.RHS_Prow[:]   = 0
            for i in range (geo.Nx-1):

                Tri.upper_Prow[i] = Faces.a_e_P[i,j]
                Tri.diag_Prow[i]  = Faces.a_P_P[i,j]
                Tri.lower_Prow[i] = Faces.a_w_P[i,j]
                Tri.RHS_Prow[i]   = -Faces.b[i,j]

                # Inlet
                if i == 0:

                    Tri.lower_Prow[i] = 0

                    #Residiual
                    Tri.RHS_Prow[i] += -Faces.a_e_P[i,j]*Fields.Pp[i+1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j]

                elif i == geo.Nx-2:

                    Tri.upper_Prow[i] = 0

                    #Residiual
                    Tri.RHS_Prow[i] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j]

                else:
                
                    #Residiual
                    Tri.RHS_Prow[i] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_P_P[i,j]*Fields.Pp[i,j] -Faces.a_e_P[i,j]*Fields.Pp[i+1,j]

                if j == 0:
                
                    Tri.RHS_Prow[i] += -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]

                elif j == geo.Ny-1:

                    Tri.RHS_Prow[i] += -Faces.a_s_P[i,j]*Fields.Pp[i,j-1] 

                else:

                    Tri.RHS_Prow[i] += -Faces.a_n_P[i,j]*Fields.Pp[i,j+1] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1]

            """
            print(f'----------- row {j}---------------')
            print('rhs results')
            print(Tri.RHS_Prow)
            print('lower results')
            print(Tri.lower_Prow)
            print('diag results')
            print(Tri.diag_Prow)
            print('upper results')
            print(Tri.upper_Prow)
            #print('rhs results')
            #print(np.array2string(np.flipud(Tri.RHS_Prow.T),formatter={'float_kind': lambda x: f"{x:8.6f}"}, max_line_width= 1000000000))   
            """

            ab = np.zeros((3, geo.Nx-1))
            
            ab[0,1:] = Tri.upper_Prow[:-1]
            ab[1,:]  = Tri.diag_Prow
            ab[2,:-1] = Tri.lower_Prow[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_Prow)
            Fields.Pp[:,j] += sol

        #print('Pp row results')
        #print(np.array2string(np.flipud(Fields.Pp.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
                        

    #print('pp tot results')
    #print(np.array2string(np.flipud(Fields.Pp.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
                    
                            