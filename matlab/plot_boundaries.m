clear all;
close all;

Nx = 32;
Ny = 32;
Nz = 32;

Nthreads = 1;

addpath('~/Programs/Matlab/ExampleCodes/WENO5');
BC = 'periodic';

x.xleft     = 0;
x.xright    = 1;
[xaxis, dx] = GenerateMesh(x, Nx);
yaxis       = xaxis;
zaxis       = xaxis;
dy          = dx;
dz          = dx;

q    = load('q.txt', 'ascii');
qall = load('qall.txt', 'ascii');


assert(length(q)    ==(Nx+6)*(Ny+6)*(Nz+6), 'Mismatch in length of q...');
assert(length(qall) ==(Nx+6)*(Ny+6)*(Nz+6), 'Mismatch in length of qall...');

Q    = zeros(Nx+6, Ny+6, Nz+6);
Qall = zeros(Nx+6, Ny+6, Nz+6);

[X,Y,Z] = ndgrid(xaxis, yaxis, zaxis);

ufh  = @(x,y,z) sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);

%
% ufh  = @(x,y,z) sin(2*pi*x).*sin(2*pi*y);
% rufh = @(x,y,z) -2*pi*cos(2*pi*x).*sin(2*pi*y)-2*pi*cos(2*pi*y).*sin(2*pi*x);

% FORTRAN uses column major order to store, so the loop ordering has to be
% kk, jj, ii
counter = 1;

for kk=1:Nz+6
    for jj=1:Ny+6
        for ii=1:Nx+6
            Q(ii,jj,kk) = q(counter);
            Qall(ii,jj,kk) = qall(counter);
            counter     = counter+1;
        end
    end
end

ind = 1;
Qplot = reshape(Q(:,ind,:), Nx+6, Ny+6);
mesh(Qplot)