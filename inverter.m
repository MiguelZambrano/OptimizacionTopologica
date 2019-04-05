%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%function top88(nelx,nely,volfrac,penal,rmin,ft)
nelx=80; nely=40; volfrac=0.3; penal=3.0; rmin=3.0; ft=1; %rmin=width*0.04
beta = 100;
eta = 0.5;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;    %Minimum Young modulus
nu = 0.3;   %Poisson's ratio
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);  %Stiffness matrix
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% DEFINE LOADS AND SUPPORTS (HALF FORCE INVERTER)
din = 1;
dout = 2*nelx*(nely+1)+1;
F = sparse(2*(nely+1)*(nelx+1),2);
F(din,1) = 1;
F(dout,2) = -1;
fixeddofs = union(2:2*(nely+1):2*(nely+1)*(nelx+1), 2*(nely+1):-1:2*(nely+1)-3);
U = zeros(2*(nely+1)*(nelx+1),2);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
loop = 0;
loopbeta = 0;
change = 1;
%% HEAVISIDE FILTER
if ft == 1 || ft == 2
    xPhys = x;
elseif ft == 3
    xTilde = x;
    %xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
end
%% INITIALIZE MMA OPTIMIZER
m = 1;                % The number of general constraints.
n = nelx*nely;             % The number of design variables x_j.
xmin = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
low = ones(n,1);        % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0 = 1;                % The constants a_0 in the term a_0*z.
a = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
obj = NaN(200,1);
changeplot = NaN(200,1);
volume = NaN(200,1);
%% START ITERATION
while change > 0.01
    loop = loop + 1;
    loopbeta = loopbeta + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    K(din,din) = K(din,din) + 1;
    K(dout,dout) = K(dout,dout) + 0.001;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    U1 = U(:,1); U2 = U(:,2);
    ce = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),[nely,nelx]);
    c = U(dout,1);
    dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);
    obj(loop+1) = c;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    elseif ft == 3
        %dx = beta*exp(-beta*xTilde)+exp(-beta);
        %%dx = (1-beta*(tanh(beta*(xTilde-eta)))*(tanh(beta*(xTilde-eta))))/(tanh(beta*eta)+tanh(beta*(1-eta)));
        dx = (beta*sech(beta*(x-eta)).^2)/(tanh(beta*eta)+tanh(beta*(1-eta)));
        dc(:) = H*(dc(:).*dx(:)./Hs);
        dv(:) = H*(dv(:).*dx(:)./Hs);
    end
    %% METHOD OF MOVING ASYMPTOTES
    fval  = sum(xPhys(:)) / (volfrac*nelx*nely) - 1;
    dfdx  = dv(:)' / (volfrac*nelx*nely);
    dfdx2 = 0*dv(:)';
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low,upp] = ...
        mmasub(m, n, loop, x(:), xmin, xmax, xold1, xold2, c, dc(:), 0*dc(:), fval, dfdx, dfdx2, low, upp, a0, a, c_MMA, d);
    % Update MMA Variables
    xnew = reshape(xmma,nely,nelx);
    change = max(abs(xnew(:)-x(:)));
    if ft == 1
        xPhys(:) = xnew(:);
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    elseif ft == 3
      xTilde(:) = (H*xnew(:))./Hs;
      %xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
      xPhys = (tanh(beta*eta)+tanh(beta*(xTilde-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
    end
    xold2 = xold1(:);
    xold1 = x(:);
    x = xnew;
    changeplot(loop+1) = change;
    volume(loop+1) = mean(xPhys(:));
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
  %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
  if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 0.01)
    beta = 2*beta;
    loopbeta = 0;
    change = 1;
    fprintf('Parameter beta increased to %g.\n',beta);
  end
end
%% EXTRA PLOTS
obj(isnan(obj)) = [];
changeplot(isnan(changeplot)) = [];
volume(isnan(volume)) = [];

figure;
ax1 = subplot(3,1,1);
plot(obj,'r');
title(ax1,'Objective function')

ax2 = subplot(3,1,2);
plot(volume,'r');
title(ax2,'Volume')

ax3 = subplot(3,1,3);
plot(changeplot,'r');
title(ax3,'Change of volume')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,    %
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,                      %
%  Technical University of Denmark,                                                                     %
%  DK-2800 Lyngby, Denmark.                                                                              %
% Please sent your comments to: sigmund@fam.dtu.dk                                      %
%                                                                                                                          %
% The code is intended for educational purposes and theoretical details             %
% are discussed in the paper                                                                                %
% "Efficient topology optimization in MATLAB using 88 lines of code,                 %
% E. Andreassen, A. Clausen, M. Schevenels,                                                       %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010                            %
% This version is based on earlier 99-line code                                                    %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,             %
% Vol 21, pp. 120--127.                                                                                        %
%                                                                                                                          %
% The code as well as a postscript version of the paper can be                           %
% downloaded from the web-site: http://www.topopt.dtu.dk                                %
%                                                                                                                          %
% Disclaimer:                                                                                                        %
% The authors reserves all rights but do not guaranty that the code is                %
% free from errors. Furthermore, we shall not be liable in any event                     %
% caused by the use of the program.                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%