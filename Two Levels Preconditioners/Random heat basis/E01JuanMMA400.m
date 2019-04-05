%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
nx = 20; ny = 20; Nx = 20; Ny = 20; volfrac = 0.3; penal = 3.0; rmin = 2; ft = 2; % rmin = 0.04*(width of domain)
nelx = nx*Nx; nely = ny*Ny;
max_it_PCG = 300; max_it_opt = 299; tol = 1e-6; lambdamax = 6;
ax = 0; bx = 1; ay = 0; by = 1;
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;    % Minimum Young modulus
nu = 0.3;   % Poisson's ratio
%% CREATE DOMAIN
nvel = (Nx*nx+1)*(Ny*ny+1);
% Create the coeficient matrix
x = repmat(volfrac,nely,nelx);
xPhys = x;
% Mesh construction
dom_ov = overlapmesh(ax,bx,ay,by,Nx,Ny,nx,ny);
dom = finemesh(ax,bx,ay,by,Nx,Ny,nx,ny);
% Define forcing term
F = sparse(2:2:2*(nely+1)*(nelx+1),ones((nely+1)*(nelx+1),1),-0.001*ones((nely+1)*(nelx+1),1),2*(nely+1)*(nelx+1),1);
A11 = [12 3 -6 -3; 3 12 3 0; -6 3 12 -3; -3 0 -3 12];
A12 = [-6 -3 0 3; -3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6];
B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4];
B12 = [ 2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2];
LEaux = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
LE = LEaux([1,3,5,7,2,4,6,8],[1,3,5,7,2,4,6,8]);
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
loop = 0;
change = 1;
obj = NaN(max_it_opt,1);
pcgiterations = NaN(max_it_opt,1);
dimcoarse = NaN(max_it_opt,1);
eigcalc = NaN(max_it_opt,1);
%% ORDERING VECTORS
ord = zeros(nelx+1,1);
ord(1) = 2*(nely+1)-1;
for i = 2:nelx+1
   ord(i) = ord(i-1)+2*(nely+1);
end
ord = repmat(ord,nely+1,1);
inc = repmat(0:2:(2*nely+1),nely+1,1);
inc = inc(:);
ord = ord-inc;
ord = [ord;ord+1];
F = F(ord);
%% INITIALIZE MMA OPTIMIZER
m = 1; % The number of general constraints.
n = nelx*nely;     % The number of design variables x_j.
xmin = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax = ones(n,1); % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);     % xval, one iteration ago (provided that iter>1).
xold2 = x(:);     % xval, two iterations ago (provided that iter>2).
low = ones(n,1); % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp = ones(n,1); % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0 = 1; % The constants a_0 in the term a_0*z.
a = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1e4*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%% START ITERATION
void=[];
save('E01-400var.mat','void')
while loop <= max_it_opt
    loop = loop + 1;
    %% FE-ANALYSIS
    % Local matrices
    [dom,b] = localsettingopt400(dom,Nx,Ny,nvel,nu,penal,E0,Emin,xPhys);
    dom_ov = localsetting_ov2opt400(dom_ov,Nx,Ny,nvel,nu,penal,E0,Emin,xPhys);
    if loop == 1 || loop == 2
        meanpcgiter = 1;
    elseif loop <= 6
        meanpcgiter = 1.2*mean(pcgiterations(1:loop-2));
    else
        meanpcgiter = 1.2*mean(pcgiterations(loop-6:loop-2)); % 20% of the last iterations mean
    end
    if loop == 1 || pcgiterations(loop-1) >= meanpcgiter || mod(loop,10) == 0
        eigcalc(loop) = 1;
        dom_ov = linearonesDB(dom_ov,Nx,Ny,b);
        %fprintf('Computing heat eigenvectors\n');
        %[dom_ov,badmat] = localeigenvectorsDBopt(b,dom_ov,Nx,Ny,lambdamax);
        fprintf('Computing heat random eigenvectors\n');
        maxrandom = 18;
        [dom_ov,badmat] = localeigenvectorsRandomDBopt(b,dom_ov,Nx,Ny,lambdamax,maxrandom);
        disp(badmat)
        %fprintf('Computing elasticity random eigenvectors\n');
        %maxrandom = 38;
        %[dom_ov,badmatE] = localeigenvectorsRandomDBEoptboundary(bE,dom_ov,Nx,Ny,nvel,lambdamax,maxrandom);
        %fprintf('Computing elasticity eigenvectors\n');
        %[dom_ov,badmatE] = localeigenvectorsDBEopt(bE,dom_ov,Nx,Ny,nvel,lambdamax);
        %disp(badmatE)
        %fprintf('Creating energy minimizing basis\n');
        dom_ov = emfpuK22(dom_ov,Nx,Ny,nvel);
        %dom_ov = emfpuK22E(dom_ov,Nx,Ny,nvel);
        %fprintf('Generating the elasticity coarse matrix (elasticity basis)...\n');
        %A0E4 = coarse_matrixBADE4(dom,dom_ov,Nx,Ny,nvel);
        fprintf('Generating the elasticity coarse matrix (heat basis)...\n');
        pause
        A0Erot=coarse_matrixBADE3(dom,dom_ov,Nx,Ny,nvel);
    end
    %dimcoarse(loop) = size(A0E4,1);
    dimcoarse(loop) = size(A0Erot,1);
    bE = F;
    %[U,error2CE,iter2CE,flag2CE,lambdamax2CE,condnumber2CE] = AEpcg_ADD_2LtwoC(bE*0,bE,bE,max_it_PCG,tol,dom,dom_ov,Nx,Ny,A0E4,nvel);
    [U,error2CE,iter2CE,flag2CE,lambdamax2CE,condnumber2CE] = AEpcg_ADD_2LtwoCRot(bE*0,bE,bE,max_it_PCG,tol,dom,dom_ov,Nx,Ny,A0Erot,nvel);
    pcgiterations(loop) = iter2CE;
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dc = zeros(nely,nelx);
    cs = zeros(nely,nelx);
    for i1 = 1:Nx
        for i2 = 1:Ny
            M = dom(i1,i2).M;
            v = dom(i1,i2).v;
            mesh = dom(i1,i2).mesh;
            Ig = dom(i1,i2).Ig;
            for i = 1:mesh.ne
                lcol = M(i,:);
                lcol2 = [Ig(lcol), nvel+Ig(lcol)];
                Ue = U(lcol2);
                DCEi = Ue'*LE*Ue;
                Vi = v(lcol(1:4),:);
                IK = floor(400*mean(Vi))+1;
                dc(IK(1),IK(2)) = -penal*(E0-Emin)*xPhys(IK(1),IK(2)).^(penal-1).*DCEi;
                cs(IK(1),IK(2)) = (Emin+xPhys(IK(1),IK(2)).^penal*(E0-Emin)).*DCEi;
            end
        end
    end
    c = sum(sum(cs));
    dv = ones(nely,nelx);
    obj(loop+1) = c;
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
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
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    elseif ft == 3
        xTilde(:) = (H*xnew(:))./Hs;
        xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    end
    xold2 = xold1(:);
    xold1 = x(:);
    x = xnew;
    %% PRINT RESULTS
    fprintf('It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,full(c),full(mean(xPhys(:))),full(change));
    %% SAVE VARIABLES
    if mod(loop,25) == 0
        S1.(strcat('badmat',num2str(loop))) = badmat;
        S2.(strcat('eigenvalues',num2str(loop))) = [dom_ov(:).lambda];
        %S3.(strcat('psibad',num2str(loop))) = [dom_ov(:).psibadE]; % Not working
        save('E01-400var.mat','S1','S2','-append')
    end
end

%% EXTRA PLOTS
obj(isnan(obj)) = [];
pcgiterations(isnan(pcgiterations)) = [];
dimcoarse(isnan(dimcoarse)) = [];
eigcalc = eigcalc(1:loop);
xaxisplot = 1:1:size(obj,1);

%% SAVE FINAL VARIABLES
save('E01-400varfinal.mat','obj','pcgiterations','eigcalc')

%% SAVE FINAL RESULTS
fid = fopen('E01output400.txt','w');
fprintf(fid,'Obj.:%11.4f Total coarse basis calculations: %3i',obj(end),sum(eigcalc,'omitnan'));
fclose(fid);

%% SAVE FIGURES
plotfig(1) = figure;
ax1 = subplot(3,1,1);
plot(obj,'r');
title(ax1,'Objective function')

ax2 = subplot(3,1,2);
plot(xaxisplot,pcgiterations,'r');
hold on;
plot(xaxisplot,mean(pcgiterations)*eigcalc,'ob','MarkerFaceColor','b');
title(ax2,'PCG iterations')

ax3 = subplot(3,1,3);
plot(xaxisplot,dimcoarse,'r');
title(ax3,'Dimension of the coarse space')

xlim([ax1 ax2 ax3],[1 size(obj,1)])

plotfig(2) = figure;
colormap(gray); imagesc(1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

savefig(plotfig,'E01-400.fig','compact')
close(plotfig)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund, Department of Solid Mechanics,             %
% Technical University of Denmark,                                         %
% DK-2800 Lyngby, Denmark.                                                 %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
