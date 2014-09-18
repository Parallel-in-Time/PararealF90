def build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose):
    with open("parameter.in", "w") as myfile:
        myfile.write("&param \n")
        myfile.write("nu="+str(nu)+"\n")
        myfile.write("Nx="+str(Nx)+"\n")
        myfile.write("Ny="+str(Ny)+"\n")
        myfile.write("Nz="+str(Nz)+"\n")
        myfile.write("N_fine="+str(N_fine)+"\n")
        myfile.write("N_coarse="+str(N_coarse)+"\n")
        myfile.write("Niter="+str(Niter)+"\n")
        myfile.write("Tend="+str(Tend)+"\n")
        myfile.write("do_io="+str(do_io)+"\n")
        myfile.write("be_verbose="+str(be_verbose)+"\n")
        myfile.write("/ \n")