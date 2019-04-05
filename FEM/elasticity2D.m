%% STATIONARY ELASTICITY, SQUARES 2D, VECTORIZED %%
nelx = 120;
nely = 30;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS
%F = sparse([2*(nely+1)*nelx+2,2*(nely+1)*(nelx+1)],[1 2],[1 -1],2*(nely+1)*(nelx+1),2);
%F = sparse(floor(0.5*(nelx+1)*(nely+1))+1,1,1,(nely+1)*(nelx+1),1);
%F = 10*ones((nely+1)*(nelx+1),1);
%fixeddofs = (nelx+1)*(nely+1)-nely:1:(nelx+1)*(nely+1);
%fixeddofs = (nely+1)*(ceil(2*(nelx+1)/5)-1)+1:nely+1:(nely+1)*(ceil(3*(nelx+1)/5)-1)+1;
%fixeddofs = 2*(nelx)*(nely+1)+1:2*(nely+1)*(nelx+1);
%F = sparse(nelx*(nely+1)+2,1,-1,2*(nely+1)*(nelx+1),1);
%fixeddofs = union(2*nely+2,2*(nelx+1)*(nely+1));
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
fixeddofs = union(1:2:2*(nely+1),2*(nelx+1)*(nely+1));
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
U = zeros(2*(nelx+1)*(nely+1),1);
%% FE-ANALYSIS
sK = repmat(KE(:),nelx*nely,1);
K = sparse(iK,jK,sK);
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
Uh = reshape(U(1:2:2*(nely+1)*(nelx+1)),nely+1,nelx+1);
Uv = reshape(U(2:2:2*(nely+1)*(nelx+1)),nely+1,nelx+1);
%% PLOT
[xx,yy] = meshgrid(1:nelx+1,1:nely+1);
%set(gcf,'Position',get(0,'Screensize'));
figure(1);
quiver(xx,yy,Uh,Uv);
%set(gcf,'Position',get(0,'Screensize'));
figure(2);
surf(Uh);
%set(gcf,'Position',get(0,'Screensize'));
figure(3);
surf(Uv);