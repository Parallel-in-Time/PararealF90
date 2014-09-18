import math
def generate_q0(Nx, Ny, Nz):
    dx = 1.0/Nx
    dy = 1.0/Ny
    dz = 1.0/Nz
    f = open('q0.dat','w')    
    for i in range(-2,Nx+4):
        for j in range(-2,Ny+4):
            for k in range(-2,Nz+4):
                x = 0.5*dx + (i-1)*dx
                y = 0.5*dy + (j-1)*dy
                z = 0.5*dy + (k-1)*dz
                q = math.sin(2.0*math.pi*x)*math.sin(2.0*math.pi*y)*math.sin(2.0*math.pi*y)
                data = "%.35f" % q
                f.write(data+"\n")
    f.close