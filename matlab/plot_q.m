clear all;
close all;

Nx = 20;
Ny = 20;
Nz = 20;

Nthreads = 1;

addpath('~/Programs/Matlab/ExampleCodes/WENO5');

x.xleft     = 0;
x.xright    = 1;
[xaxis, dx] = GenerateMesh(x, Nx);
[yaxis, dy] = GenerateMesh(x, Ny);
[zaxis, dz] = GenerateMesh(x, Nz);

q    = load('q.txt','ascii');
qref = load('qref.txt','ascii');

assert( length(q)==Nx*Ny*Nz, 'Mismatch in length of q ...');
assert( length(qref)==Nx*Ny*Nz, 'Mismatch in length of qref ...');

Q    = zeros(Nx, Ny, Nz);
Qref = zeros(Nx, Ny, Nz);

counter = 1;
for nn=1:Nthreads
    for kk=1:Nz
        for jj=1:Ny
            for ii=1:Nx
                Q(ii,jj,kk)    = q(counter);
                Qref(ii,jj,kk) = qref(counter);
                counter     = counter+1;
            end
        end
    end
end

max_def_j = 0;
for jj=1:Ny-1
    max_def_j = max(max_def_j, max(max(abs(Q(:,jj+1,:)-Q(:,jj,:)))));
end 
fprintf('Maximum variation in j direction: %5.2f \n', max_def_j);

max_def_k = 0;
for kk=1:Nz-1
   max_def_k = max(max_def_k, max(max(abs(Q(:,:,kk+1)-Q(:,:,kk)))));
end
fprintf('Maximum variation in k direction: %5.2f \n', max_def_k);

figure(1)
mesh(yaxis, xaxis, Q(:,:,1)); view(-90,0)
title('Q');

figure(2)
mesh(yaxis, xaxis, Qref(:,:,1)); view(-90,0)
title('Qref');

Qdiff = Q(:,:,:,1) - Qref(:,:,:,1);

figure(3)
mesh(yaxis, xaxis, Qdiff(:,:,1)); view(-90,0);
title('Qdiff');