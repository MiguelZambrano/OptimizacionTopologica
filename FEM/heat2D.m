%% STATIONARY HEAT CONDUCTION, SQUARES 2D, VECTORIZED %%
nelx = 40;
nely = 30;
%% PREPARE FINITE ELEMENT ANALYSIS
KE = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4];
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+1 nely -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS
%F = 10*ones((nely+1)*(nelx+1),1);
%fixeddofs = (nelx+1)*(nely+1)-nely:1:(nelx+1)*(nely+1);
F = sparse(floor(0.5*(nelx+1)*(nely+1))+1,1,1,(nely+1)*(nelx+1),1);
fixeddofs = (nely+1)*(ceil(2*(nelx+1)/5)-1)+1:nely+1:(nely+1)*(ceil(3*(nelx+1)/5)-1)+1;
alldofs = 1:(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
U = zeros((nelx+1)*(nely+1),1);
%% FE-ANALYSIS
sK = repmat(KE(:),nelx*nely,1);
K = sparse(iK,jK,sK);
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
U = reshape(U,nely+1,nelx+1);
%% PLOT
%set(gcf,'Position',get(0,'Screensize'));
surf(U);