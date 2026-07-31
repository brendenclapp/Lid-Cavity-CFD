import numpy as np
from scipy.linalg import solve_banded

def Coeff_P(geo, var, Faces, Fields):

    Faces.b[:] = 0.0
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

                Faces.b[i,j] += (var.rho*Fields.u[i+1,j])*geo.dy

            elif i == geo.Nx-1:

                Faces.b[i,j] += -(var.rho*Fields.u[i,j])*geo.dy

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



def TDMA_P(geo, Faces, Fields, Tri):

    Fields.Pp[:] = 0.0

    for ittp in range (20):
       
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
                Tri.RHS_P[i]   = -Faces.b[i,j]

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
            Fields.Pp[:,j] += sol
            
        #========================================== COLUMNS ============================
        for i in range (geo.Nx):
            Tri.upper_P[:] = 0
            Tri.diag_P[:]  = 0
            Tri.lower_P[:] = 0
            Tri.RHS_P[:]   = 0
            for j in range (geo.Ny):

                Tri.upper_P[j] = Faces.a_n_P[i,j]
                Tri.diag_P[j]  = Faces.a_P_P[i,j]
                Tri.lower_P[j] = Faces.a_s_P[i,j]
                Tri.RHS_P[j]   = -Faces.b[i,j]

                if i == 0:
                    
                    Tri.RHS_P[j] += (-Faces.a_e_P[i,j]*Fields.Pp[i+1,j])

                elif i == geo.Nx-1:

                    Tri.RHS_P[j] += (-Faces.a_w_P[i,j]*Fields.Pp[i-1,j])

                else:
                            
                    Tri.RHS_P[j] += -Faces.a_w_P[i,j]*Fields.Pp[i-1,j] - Faces.a_e_P[i,j]*Fields.Pp[i+1,j]


                if j == 0:

                    Tri.lower_P[j] = 0
                    
                    #Residual
                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_n_P[i,j]*Fields.Pp[i,j+1]

                elif j == geo.Ny-1:

                    Tri.upper_P[j] = 0

                    #Residual
                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1]

                else:

                    Tri.RHS_P[j] += -Faces.a_P_P[i,j]*Fields.Pp[i,j] - Faces.a_s_P[i,j]*Fields.Pp[i,j-1] -Faces.a_n_P[i,j]*Fields.Pp[i,j+1]
            

            ab = np.zeros((3, geo.Ny))
            
            ab[0,1:] = Tri.upper_P[:-1]
            ab[1,:]  = Tri.diag_P
            ab[2,:-1] = Tri.lower_P[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_P)
            Fields.Pp[i,:] += sol
