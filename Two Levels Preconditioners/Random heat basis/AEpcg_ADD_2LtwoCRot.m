function [x, error, iter, flag, lambdamax, condnumber] = AEpcg_ADD_2LtwoCRot(x,b,br,max_it,tol,dom,dom_ov,Nx,Ny,A0,nvel)
%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; et templates.ps).
%
%  [x, error, iter, flag] = cg(A, x, b, M, max_it, tol)
%
% cg.m solves the symmetric positive definite linear system Ax=b
% using the Conjugate Gradient method with preconditioning.
%
% input   A        REAL symmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         M        REAL preconditioner matrix
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

flag = 0;                                 % initialization
iter = 0;

bnrm2 = norm(br);
if  (bnrm2 == 0.0), bnrm2 = 1.0; end
Ax=apply_AE(x,dom,Nx,Ny,nvel);
% r=0*b;
r = b-Ax;
r1=zerodir(r(1:nvel),dom,Nx,Ny);
r2=zerodir(r(1+nvel:2*nvel),dom,Nx,Ny);
r=[r1;r2];
%  plot_vector(r,dom,Nx,Ny)

%  error = norm(r) / bnrm2;
bnrm2=norm(r);
error=1;

fprintf('PCG2L2C residual(%d) = %g\n',iter,error)
if (error < tol), return, end

for iter = 1:max_it
    
    zlocal = apply_AE2ladd_inv(r,dom_ov,Nx,Ny,nvel);
    z0 = apply_P0BADE3(r,dom_ov,Nx,Ny,A0);
    
    z=zlocal+z0;    

    rho = (r'*z);
    
    if (iter > 1)
        beta(iter) = rho / rho_1;
        p = z + beta(iter)*p;
    else
        p = z;
    end
    Ap= apply_AE(p,dom,Nx,Ny,nvel);
    %     q=0*p;
    q=Ap;
    alpha(iter) = rho / (p'*q);
    x = x + alpha(iter) * p;                    % update approximation vector
    %       subplot(2,2,3)
    %       plot_vector(alpha(iter)*p, dom,Nx,Ny,mu,[]); title('\alpha p')
    %       hold off
    %
    %      subplot(2,2,4)
    %      plot_vectorE(x, dom,Nx,Ny);
    %      title('x')
    %      pause
    %      hold off
    % pause(0.2)
    
    % % % %  I0x=apply_I0new(x,dom,dom_ov,Nx,Ny,mu,k1perm);
    % % % %  Ax_I0= applay_A(x,dom,Nx,Ny);
    % % % %  AI0x=applay_A(I0x,dom,Nx,Ny);
    % % % %  AxmI0x=applay_A(x-I0x,dom,Nx,Ny);
    % % % %  CAI0A=full(dot(I0x,AI0x)/dot(x,Ax_I0));
    % % % %  CAI0A2=full(dot(x-I0x,AxmI0x)/dot(x,Ax_I0));
    % % % %  fprintf('Cstab=%.2f, Capprox=%.2f \n',CAI0A,CAI0A2);
    
    
    r = r - alpha(iter)*q;                      % compute residual
    error = norm(r) / bnrm2;            % check convergence
    fprintf('[ %f, %f ] PCG2L2C residual(%d) = %g\n',norm(r),bnrm2,iter,error)
    if (error <= tol), break, end
    
    rho_1 = rho;
    %      plot_vector(x,dom,Nx,Ny)
    %      MovieCG(iter)=getframe;
    %plot_vector(x,dom_ov,Nx,Ny);
    %axis([-1,1,-1,1,-1,1])
    %pause(0.2)
    %hold off
end
%
%  compute eigenvalues and codition number
%
d(1)= 1/alpha(1);
for i = 2: iter
    d(i)=beta(i)/alpha(i-1)+1/alpha(i);
end
for i = 1: iter-1
    s(i)=-1*sqrt(beta(i+1))/alpha(i);
end
%
T = sparse(zeros(iter,iter));
T(1,1)= d(1);
for i=2:iter
    T(i,i) = d(i);
end
for i = 1:iter-1
    T(i,i+1)= s(i); T(i+1,i) = T(i,i+1);
end
lambda = eig(full(T));
lambdamax = max(lambda);
lambdamin = min(lambda);
condnumber = lambdamax/lambdamin;

if (error > tol)
    flag = 1; end         % no convergence

% END cg.m
