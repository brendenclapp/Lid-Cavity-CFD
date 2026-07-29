    
import numpy as np



def Coeff_u(geo,var,Fields,Faces):
   
   
    for i in range (1,geo.Nx):
        for j in range (geo.Ny):

            #===================================== EAST ======================================================== 
            # East A Interior
        
            ce = var.rho*(Fields.u[i+1,j] + Fields.u[i,j])*0.5
            de = var.mu*geo.dy/geo.dx

            Faces.a_e_u[i,j] = -((abs(ce)-ce)/2)*geo.dy - de

            #===================================== WEST ======================================================== 
            # West A Interior

            cw = var.rho*(Fields.u[i-1,j] + Fields.u[i,j])*0.5
            dw = var.mu * geo.dy / geo.dx

            Faces.a_w_u[i,j] = -((abs(cw)+cw)/2)*geo.dy - dw

            #===================================== NORTH ======================================================== 
            #Bottom Wall
            if j == 0:

                cn = var.rho*(Fields.v[i,j+1] + Fields.v[i-1,j+1])*0.5
                dn = var.mu * geo.dx / geo.dy
                wn = (var.mu * geo.dx)/ 3*geo.dy
                
                Faces.a_n_u[i,j] = -((abs(cn)-cn)/2)*geo.dx - dn - wn

            #Top Wall
            elif j == geo.Ny-1:

                Faces.a_n_u[i,j] = 0

            # North A Interior
            else:
            
                cn = var.rho*(Fields.v[i,j+1] + Fields.v[i-1,j+1])*0.5
                dn = var.mu * geo.dx / geo.dy

                Faces.a_n_u[i,j] = -((abs(cn)-cn)/2)*geo.dx - dn

            #===================================== SOUTH ======================================================== 
            #Bottom Wall
            if j == 0:

                Faces.a_s_u[i,j] = 0

            #Top Wall
            elif j == geo.Ny-1:

                cs = var.rho*(Fields.v[i,j] + Fields.v[i-1,j])*0.5
                ds = var.mu * geo.dx / geo.dy
                ws = (var.mu * geo.dx)/ 3*geo.dy
                
                Faces.a_s_u[i,j] = -((abs(cs)+cs)/2)*geo.dx - ds - ws

            
            # South A Interior
            else:
                cs = var.rho*(Fields.v[i,j] + Fields.v[i-1,j])*0.5
                ds = var.mu * geo.dx / geo.dy

                Faces.a_s_u[i,j] = -((abs(cs)+cs)/2)*geo.dx - ds

            #===================================== CENTER ======================================================== 
            #Bottom Wall
            if j == 0:

                Faces.a_P_u[i,j] = (((abs(ce)-ce)/2)*geo.dy + de) + (((abs(cw)+cw)/2)*geo.dy + dw) \
                                + (((abs(cn)-cn)/2)*geo.dx + dn + (3*var.mu*geo.dx)/geo.dy)

            #Top Wall
            elif j == geo.Ny-1:

                Faces.a_P_u[i,j] = (((abs(ce)-ce)/ 2)*geo.dy + de) + (((abs(cw)+cw)/2)*geo.dy + dw) \
                                    + (((abs(cs)+cs)/2)*geo.dx + ds + (3*var.mu*geo.dx)/geo.dy)

            # Central A Interior
            else:
                Faces.a_P_u[i,j] = (((abs(ce)-ce)/2)*geo.dy + de) + (((abs(cw)+cw)/2)*geo.dy + dw) \
                                + (((abs(cn)-cn)/2)*geo.dx + dn) + ((abs(cs)+cs)/2)*geo.dx + ds

            #===================================== PRESSURE ========================================================

            if j == geo.Ny - 1:

                Faces.dPdx[i,j] = ((Fields.P[i-1,j] - Fields.P[i,j]) * geo.dy) + ((8/(3*geo.dy)*var.mu*var.u_lid*geo.dx))

            else:

                Faces.dPdx[i,j] = (Fields.P[i-1,j] - Fields.P[i,j]) * geo.dy

    """
    print('a_w_u')
    print(np.array2string(np.flipud(Faces.a_w_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_e_u')
    print(np.array2string(np.flipud(Faces.a_e_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_n_u')
    print(np.array2string(np.flipud(Faces.a_n_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_s_u')
    print(np.array2string(np.flipud(Faces.a_s_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('a_P_u')
    print(np.array2string(np.flipud(Faces.a_P_u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    """

    
def TDMA_u(var, geo, Fields, Faces, Tri):

    from scipy.linalg import solve_banded

    #========================================== COLUMNS ============================
    for uitt in range (2):
        for i in range(1, geo.Nx):

            Tri.upper_col[:] = 0
            Tri.diag_col[:]  = 0
            Tri.lower_col[:] = 0
            Tri.RHS_col[:]   = 0

            for j in range(geo.Ny):

                k = j

                Tri.upper_col[k] = Faces.a_n_u[i,j]
                Tri.diag_col[k]  = Faces.a_P_u[i,j]*(1.2)
                Tri.lower_col[k] = Faces.a_s_u[i,j]
                Tri.RHS_col[k] = Faces.dPdx[i,j]


                if i == 1:

                 
                    Tri.RHS_col[k] += (-Faces.a_e_u[i,j]*Fields.u[i+1,j])

                elif i == geo.Nx-1:

                  
                    Tri.RHS_col[k] += (-Faces.a_w_u[i,j]*Fields.u[i-1,j])

                else:

                    
                    Tri.RHS_col[k] += (-Faces.a_w_u[i,j]*Fields.u[i-1,j] + -Faces.a_e_u[i,j]*Fields.u[i+1,j])


    
                if j == 0:

                    Tri.lower_col[k] = 0

                    #Residual
                    Tri.RHS_col[k] += -Faces.a_P_u[i,j]*Fields.u[i,j] -Faces.a_n_u[i,j]*Fields.u[i,j+1]

                elif j == geo.Ny-1:

                    Tri.upper_col[k] = 0

                    Tri.RHS_col[k] += Faces.a_n_u[i,j]*var.u_lid

                    #Residual
                    Tri.RHS_col[k] += -Faces.a_P_u[i,j]*Fields.u[i,j] -Faces.a_s_u[i,j]*Fields.u[i,j-1]

                else:

                    #Residual
                    Tri.RHS_col[k] += -Faces.a_P_u[i,j]*Fields.u[i,j] -Faces.a_s_u[i,j]*Fields.u[i,j-1] -Faces.a_n_u[i,j]*Fields.u[i,j+1]



            ab = np.zeros((3, geo.Ny))

            ab[0,1:] = Tri.upper_col[:-1]
            ab[1,:]  = Tri.diag_col
            ab[2,:-1] = Tri.lower_col[1:]
            sol = solve_banded((1,1), ab, Tri.RHS_col)
            Fields.u_tilde[i,:] = sol

        for i in range(1, geo.Nx):
            for j in range(geo.Ny):

                Fields.u[i,j] = Fields.u[i,j] + Fields.u_tilde[i,j]
        
                        
    # ================================================= ROWS =========================================================
        for j in range(geo.Ny):
            Tri.upper_row[:] = 0
            Tri.diag_row[:]  = 0
            Tri.lower_row[:] = 0
            Tri.RHS_row[:]   = 0
            for i in range(1, geo.Nx):

                k = i-1

                Tri.upper_row[k] = Faces.a_e_u[i,j]
                Tri.diag_row[k]  = Faces.a_P_u[i,j]*(1.2)
                Tri.lower_row[k] = Faces.a_w_u[i,j]
                Tri.RHS_row[k]   = Faces.dPdx[i,j]

                if i == 1:

                    Tri.lower_row[k] = 0

                    #Residiual
                    Tri.RHS_row[k] += -Faces.a_e_u[i,j]*Fields.u[i+1,j] - Faces.a_P_u[i,j]*Fields.u[i,j]

                elif i == geo.Nx-1:

                    Tri.upper_row[k] = 0

                    #Residiual
                    Tri.RHS_row[k] += -Faces.a_w_u[i,j]*Fields.u[i-1,j] - Faces.a_P_u[i,j]*Fields.u[i,j]

                else:

                    #Residiual
                    Tri.RHS_row[k] += -Faces.a_w_u[i,j]*Fields.u[i-1,j] - Faces.a_P_u[i,j]*Fields.u[i,j] -Faces.a_e_u[i,j]*Fields.u[i+1,j]

                if j == 0:

                    Tri.RHS_row[k] += -Faces.a_n_u[i,j]*Fields.u[i,j+1]

                elif j == geo.Ny-1:

                    Tri.RHS_row[k] += -Faces.a_s_u[i,j]*Fields.u[i,j-1] + -Faces.a_n_u[i,j]*var.u_lid

                else:

                    Tri.RHS_row[k] += -Faces.a_n_u[i,j]*Fields.u[i,j+1] - Faces.a_s_u[i,j]*Fields.u[i,j-1]

            ab = np.zeros((3, geo.Nx-1))

            ab[0,1:] = Tri.upper_row[:-1]
            ab[1,:]  = Tri.diag_row
            ab[2,:-1] = Tri.lower_row[1:geo.Nx-1]

            sol = solve_banded((1,1), ab, Tri.RHS_row)

            Fields.u_tilde[1:geo.Nx,j] = sol

        for i in range (1, geo.Nx):
                for j in range (geo.Ny):
                    Fields.u[i,j] = Fields.u[i,j] + Fields.u_tilde[i,j]

  
    #print('u_TDMA results')
    #print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:8.5f}"}, max_line_width= 1000000000))