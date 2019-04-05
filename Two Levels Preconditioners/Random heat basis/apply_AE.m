function w=apply_AE(u,dom,Nx,Ny,nvel)
w=zeros(size(u,1),size(u,2));
for i1=1:Nx
    for i2=1:Ny
        A=dom(i1,i2).AE;
        free=dom(i1,i2).free;
        mesh=dom(i1,i2).mesh;
        AI=A([free,free+mesh.nv],[free,free+mesh.nv]);
        Igfree=dom(i1,i2).Igfree;   
        lu=u([Igfree,Igfree+nvel]);
        lAu=AI*lu;
        w([Igfree,Igfree+nvel])=w([Igfree,Igfree+nvel])+lAu;
    end
end