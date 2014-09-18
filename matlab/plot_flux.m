clear all;
close all;

Nx = 60;
Ny = 15;
Nz = 15;

Nthreads = 1;

addpath('~/Programs/Matlab/ExampleCodes/WENO5');
BC = 'periodic';

x.xleft     = 0;
x.xright    = 1;
[xaxis, dx] = GenerateMesh(x, Nx);
xaxis_int   = linspace(x.xleft, x.xright, Nx+1);
[yaxis, dy] = GenerateMesh(x, Ny);
[zaxis, dz] = GenerateMesh(x, Nz);

q    = load('q.txt', 'ascii');
fq   = load('fluxi.txt', 'ascii');
rq   = load('rq.txt', 'ascii');

assert(length(q)  == Nx*Ny*Nz*Nthreads, 'Mismatch in length of q...');
assert(length(fq) == (Nx+1)*Ny*Nz*Nthreads, 'Mismatch in length of fq...');
assert(length(rq) == Nx*Ny*Nz*Nthreads, 'Mismatch in length of rq ...');

Q    = zeros(Nx,   Ny, Nz);
FQ   = zeros(Nx+1, Ny, Nz);
RQ   = zeros(Nx,   Ny, Nz);

[X,Y,Z] = ndgrid(xaxis, yaxis, zaxis);
[X_int, Y_int, Z_int] = ndgrid(xaxis_int, yaxis, zaxis);

ufh  = @(x,y,z) sin(2*pi*x);% .* sin(2*pi*y) .* sin(2*pi*z);
fufh = @(x,y,z) 0.5*sin(2*pi*x).^2;
rqfh = @(x,y,z) -2*pi*sin(2*pi*x).*cos(2*pi*x);

% FORTRAN uses column major order to store, so the loop ordering has to be
% kk, jj, ii
counter = 1;
for nn=1:Nthreads
    for kk=1:Nz
        for jj=1:Ny
            for ii=1:Nx
                Q(ii,jj,kk) = q(counter);
                RQ(ii,jj,kk) = rq(counter);
                counter     = counter+1;
            end
        end
    end
end

Qh = ufh(X,Y,Z);
Qdiff = Q - Qh;
assert( max(max(max(abs(Qdiff))))<1e-13, 'Computed Q and retrieved Q do not coincide...');

counter = 1;
for nn=1:Nthreads
    for kk=1:Nz
        for jj=1:Ny
            for ii=1:Nx+1
                FQ(ii,jj,kk) = fq(counter);
                counter     = counter+1;
            end
        end
    end
end


FQh = fufh(X_int,Y_int,Z_int);

figure(1);
mesh(yaxis, xaxis_int, FQ(:,:,5)); view(-90,0);
title('Retrieved FQ');

figure(2);
mesh(yaxis, xaxis_int, FQh(:,:,5)); view(-90,0);
title('Computed FQ');

figure(3);
mesh(yaxis, xaxis_int, FQ(:,:,5) - FQh(:,:,5)); view(-90,0);
title('Difference');

RQh = rqfh(X,Y,Z);
figure(4);
mesh(yaxis, xaxis, RQ(:,:,5)); view(-90,0);
title('Retrieved RQ');

figure(5);
mesh(yaxis, xaxis, RQh(:,:,5)); view(-90,0);
title('Computed RQ');