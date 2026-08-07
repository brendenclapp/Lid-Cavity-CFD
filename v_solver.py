
import numpy as np
from scipy.linalg import solve_banded

def Coeff_v(geo, var, Fields, Faces):


    for i in range (geo.Nx):
        for j in range (1, geo.Ny):


            #===================================== EAST ======================================================== 
        
                # East face at Outlet
                if i == geo.Nx-1:

                    Faces.a_e_v[i,j] = 0

                # East A Interior
                else:
                    ce = var.rho*(Fields.u[i,j] + Fields.u[i,j-1])*0.5
                    de = var.mu*geo.dy/geo.dx
        
                    Faces.a_e_v[i,j] = -((abs(ce)-ce)/2)*geo.dy - de
    
            #===================================== WEST ======================================================== 
                # West Inlet
                if i == 0:

                    Faces.a_w_v[i,j] = 0

                # West A Interior
                else:
                    cw = var.rho*(Fields.u[i,j] + Fields.u[i,j-1])*0.5
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
         
                Faces.a_P_v[i,j] =  -(Faces.a_e_v[i,j] + Faces.a_w_v[i,j] + Faces.a_n_v[i,j] + Faces.a_s_v[i,j])

            #===================================== PRESSURE ========================================================

                Faces.dPdy[i,j] = (Fields.P[i,j-1] - Fields.P[i,j]) * geo.dy

    var.res_v = 0.0
    for i in range (1, geo.Nx-1):
        for j in range (1, geo.Ny):

            var.res_v = var.res_v + ((Faces.a_P_v[i,j]*Fields.v[i,j]) + (Faces.a_w_v[i,j]*Fields.v[i-1,j]) \
                    + (Faces.a_e_v[i,j]*Fields.v[i+1,j]) + (Faces.a_s_v[i,j]*Fields.v[i,j-1]) \
                    + (Faces.a_n_v[i,j]*Fields.v[i,j+1]) - Faces.dPdy[i,j])**2

    var.res_v = np.sqrt(var.res_v)       


def TDMA_v(geo, Fields, Faces, Tri):
     
        #========================================== COLUMNS ============================
        for uitt in range (2):
            for i in range(geo.Nx):
    
                Tri.upper_colv[:] = 0
                Tri.diag_colv[:]  = 0
                Tri.lower_colv[:] = 0
                Tri.RHS_colv[:]   = 0
    
                for j in range(1,geo.Ny):
    
                    k = j - 1
    
                    Tri.upper_colv[k] = Faces.a_n_v[i,j]
                    Tri.diag_colv[k]  = Faces.a_P_v[i,j]*(1.2)
                    Tri.lower_colv[k] = Faces.a_s_v[i,j]
                    Tri.RHS_colv[k] = Faces.dPdy[i,j]
    
    
                    if i == 0:
    
                       
                        Tri.RHS_colv[k] += (-Faces.a_e_v[i,j]*Fields.v[i+1,j]) 
    
                    elif i == geo.Nx-1:
    
                        
                        Tri.RHS_colv[k] += (-Faces.a_w_v[i,j]*Fields.v[i-1,j])
    
                    else:
    
                       
                        Tri.RHS_colv[k] += (-Faces.a_w_v[i,j]*Fields.v[i-1,j] + -Faces.a_e_v[i,j]*Fields.v[i+1,j])
    
    
        
                    if j == 1:
    
                        Tri.lower_colv[k] = 0
    
                        #Residual
                        Tri.RHS_colv[k] += -Faces.a_P_v[i,j]*Fields.v[i,j] -Faces.a_n_v[i,j]*Fields.v[i,j+1]
    
                    elif j == geo.Ny-1:
    
                        Tri.upper_colv[k] = 0
    
                        Tri.RHS_colv[k] += Faces.a_n_v[i,j]*0
    
                        #Residual
                        Tri.RHS_colv[k] += -Faces.a_P_v[i,j]*Fields.v[i,j] -Faces.a_s_v[i,j]*Fields.v[i,j-1]
    
                    else:
    
                        #Residual
                        Tri.RHS_colv[k] += -Faces.a_P_v[i,j]*Fields.v[i,j] -Faces.a_s_v[i,j]*Fields.v[i,j-1] -Faces.a_n_v[i,j]*Fields.v[i,j+1]
    
    
    
                ab = np.zeros((3, geo.Ny-1))
    
                ab[0,1:] = Tri.upper_colv[:-1]
                ab[1,:]  = Tri.diag_colv
                ab[2,:-1] = Tri.lower_colv[1:]
                sol = solve_banded((1,1), ab, Tri.RHS_colv)
                Fields.v_tilde[i,1:-1] = sol
    
            for i in range(geo.Nx):
                for j in range(1, geo.Ny):
    
                    Fields.v[i,j] = Fields.v[i,j] + Fields.v_tilde[i,j]
                    
            #print('v col results')
            #print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
                                
            #================================================= ROWS =========================================================

            for j in range(1, geo.Ny):
                Tri.upper_rowv[:] = 0
                Tri.diag_rowv[:]  = 0
                Tri.lower_rowv[:] = 0
                Tri.RHS_rowv[:]   = 0
                for i in range(geo.Nx):
    
                    k = i
    
                    Tri.upper_rowv[k] = Faces.a_e_v[i,j]
                    Tri.diag_rowv[k]  = Faces.a_P_v[i,j]*(1.2)
                    Tri.lower_rowv[k] = Faces.a_w_v[i,j]
                    Tri.RHS_rowv[k]   = Faces.dPdy[i,j]
    
                    if i == 0:
    
                        Tri.lower_rowv[k] = 0
    
                        #Residiual
                        Tri.RHS_rowv[k] += -Faces.a_e_v[i,j]*Fields.v[i+1,j] - Faces.a_P_v[i,j]*Fields.v[i,j]
    
                    elif i == geo.Nx-1:
    
                        Tri.upper_rowv[k] = 0
    
                        #Residiual
                        Tri.RHS_rowv[k] += -Faces.a_w_v[i,j]*Fields.v[i-1,j] - Faces.a_P_v[i,j]*Fields.v[i,j]
    
                    else:
    
                        #Residiual
                        Tri.RHS_rowv[k] += -Faces.a_w_v[i,j]*Fields.v[i-1,j] - Faces.a_P_v[i,j]*Fields.v[i,j] -Faces.a_e_v[i,j]*Fields.v[i+1,j]
    
                    if j == 1:
    
                        Tri.RHS_rowv[k] += -Faces.a_n_v[i,j]*Fields.v[i,j+1]
    
                    elif j == geo.Ny-1:
    
                        Tri.RHS_rowv[k] += -Faces.a_s_v[i,j]*Fields.v[i,j-1] + -Faces.a_n_v[i,j]*0
    
                    else:
    
                        Tri.RHS_rowv[k] += -Faces.a_n_v[i,j]*Fields.v[i,j+1] - Faces.a_s_v[i,j]*Fields.v[i,j-1]
    
                ab = np.zeros((3, geo.Nx))
    
                ab[0,1:] = Tri.upper_rowv[:-1]
                ab[1,:]  = Tri.diag_rowv
                ab[2,:-1] = Tri.lower_rowv[1:]
    
                sol = solve_banded((1,1), ab, Tri.RHS_rowv)
    
                Fields.v_tilde[:, j] = sol
    
            for i in range(geo.Nx):
                    for j in range(1, geo.Ny):
                        Fields.v[i,j] = Fields.v[i,j] + Fields.v_tilde[i,j]

            #print('v row results')
            #print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
                    
        #print('v tot results')
        #print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
        
    
      
        