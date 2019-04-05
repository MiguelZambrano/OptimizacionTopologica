function [dom,b]=localsettingopt400(dom,Nx,Ny,nvel,nu,penal,E0,Emin,xPhys)

b=sparse(nvel,1);;

fprintf('Creating fine local matrices\n');
for i1=1:Nx
    for i2=1:Ny
        M=dom(i1,i2).M;
        v=dom(i1,i2).v;
        mesh=dom(i1,i2).mesh;
        %free=dom(i1,i2).free;
        [Asd,Asde]=Nmatrixopt400(M,v,mesh,nu,penal,E0,Emin,xPhys);
        [Massd,~]=NMassmatrixopt400(M,v,dom,mesh,penal,E0,Emin,xPhys);
        dom(i1,i2).Mass=Massd; %mass matrix for heat
        dom(i1,i2).A=Asd; % stiffness matrix for heat
        %dom(i1,i2).b=bsd(free); %right hand side for heat
        %dom(i1,i2).MassE=MassdE; % mass matrix for elasticity
        dom(i1,i2).AE=Asde; % stiffness matrix for elasticity
        %dom(i1,i2).be=bsde(free); %right hand side for elasticity
        %colI=dom(i1,i2).Igfree;
        %b(colI)=b(colI)+bsd(free); % assembling of right hand side
        %bE(colI)=bE(colI)+bsde(free);
    end
end
