function dom=localsetting_ov2opt400(dom,Nx,Ny,nvel,nu,penal,E0,Emin,xPhys)

%b=sparse(nvel,1);

fprintf('Creating matrices for the coarse grid\n');
for i1=1:Nx-1
    for i2=1:Ny-1
        M=dom(i1,i2).M;
        v=dom(i1,i2).v;
        mesh=dom(i1,i2).mesh;
        [Asd,Asde]=Nmatrixopt400(M,v,mesh,nu,penal,E0,Emin,xPhys);
        [Massd,Mass1]=NMassmatrixopt400(M,v,dom,mesh,penal,E0,Emin,xPhys);
        dom(i1,i2).Mass = Massd; % mass matrix for heat
        %dom(i1,i2).MassE = MassdE; % mass matrix for elasticity
        dom(i1,i2).Mass1 = Mass1; % mass matrix with coefficient 1
        %dom(i1,i2).Mass1E = Mass1E; % mass matrix with coefficient 1 for elasticity
        dom(i1,i2).A = Asd; % stiffness matrix for heat
        dom(i1,i2).AE = Asde; % stiffness matrix for elasticity
    end
end
