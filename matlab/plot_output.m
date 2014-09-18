clear all;
close all;

Nx = 28;
Ny = 27;
Nz = 25;

Nthreads = 1;

addpath('~/Programs/Matlab/ExampleCodes/WENO5');
BC = 'periodic';

x.xleft     = 0;
x.xright    = 1;
[xaxis, dx] = GenerateMesh(x, Nx);
[yaxis, dy] = GenerateMesh(x, Ny);
[zaxis, dz] = GenerateMesh(x, Nz);
dy          = dx;
dz          = dx;

q    = load('qfield.txt', 'ascii');
rq   = load('rqfield.txt', 'ascii');
fqxy = load('fqxyfield.txt','ascii');

assert(length(rq)  ==Nx*Ny*Nz*Nthreads, 'Mismatch in length of rq...');
assert(length(q)   ==Nx*Ny*Nz*Nthreads, 'Mismatch in length of q...');
assert(length(fqxy)==Nx*Ny*Nz*Nthreads, 'Mismatch in length of fqxz...');

Q    = zeros(Nx, Ny, Nz);
RQ   = zeros(Nx, Ny, Nz);
FQXY = zeros(Nx, Ny, Nz);

[X,Y,Z] = ndgrid(xaxis, yaxis, zaxis);

ufh = @(x,y,z) sin(2*pi*x) .* sin(2*pi*y) .* sin(2*pi*z);
%rufh = @(x,y,z) -4*pi^2*sin(2*pi*x);
%rufh = @(x,y,z) -2*pi*cos(2*pi*x).*sin(2*pi*x);

% --- linear ---
rufh = @(x,y,z) -2*pi*cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z) ...
                -2*pi*sin(2*pi*x).*cos(2*pi*y).*sin(2*pi*z) ...
                -2*pi*sin(2*pi*x).*sin(2*pi*y).*cos(2*pi*z);

% --- nonlinear ---
% rufh = @(x,y,z) -2*pi*sin(2*pi*x).*cos(2*pi*x).*sin(2*pi*y).^2.*sin(2*pi*z).^2 ...
%                 -2*pi*sin(2*pi*y).*cos(2*pi*y).*sin(2*pi*x).^2.*sin(2*pi*z).^2 ...
%                 -2*pi*sin(2*pi*z).*cos(2*pi*z).*sin(2*pi*x).^2.*sin(2*pi*y).^2;


% FORTRAN uses column major order to store, so the loop ordering has to be
% kk, jj, ii
counter = 1;
for nn=1:Nthreads
    for kk=1:Nz
        for jj=1:Ny
            for ii=1:Nx
                Q(ii,jj,kk) = q(counter);
                counter     = counter+1;
            end
        end
    end
end

Qh = ufh(X,Y,Z);

figure(1);
mesh(yaxis, xaxis, Q(:,:,5)); view(0,90);
title('Retrieved Q');

figure(2);
mesh(yaxis, xaxis, Qh(:,:,5)); view(0,90);
title('Computed Q');

Qdiff = Q - Qh;
assert( max(max(max(abs(Qdiff))))<1e-13, 'Computed Q and retrieved Q do not coincide...');

counter = 1;
% see comment above
for kk=1:Nz
    for jj=1:Ny
        for ii=1:Nx
            
            RQ(ii,jj,kk)    = rq(counter);
            FQXY(ii,jj,kk)  = fqxy(counter);
            counter         = counter+1;
        end
    end
end

figure(3)
mesh(yaxis, xaxis, RQ(:,:,5)); view(0,90);
title('Retrieved discrete RQ');

figure(4)
mesh(yaxis, xaxis, FQXY(:,:,5)); view(0,90);
title('Retrieved analytical RQ');

RQfh = rufh(X,Y,Z);

figure(5)
mesh(yaxis, xaxis, RQfh(:,:,5)); view(0,90);
title('Computed analytical RQ');

RQdiff = FQXY - RQfh;
assert( max(max(max(abs(RQdiff))))<1e-13,'Computed exact RQ and retrieved FQXY do not coincide...');

% IMPORTANT: Note that now the x-axis corresponds to the first dimension
% of Q and is therefore plotted along the vertical axis... when plotting
% with mesh etc, the roles of x and y are thus interchanged.



