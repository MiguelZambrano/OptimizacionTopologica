function w=apply_AE2ladd_inv(u,dom_ov,Nx,Ny,nvel)
w=zeros(size(u,1),size(u,2));
for i1=1:Nx-1
    for i2=1:Ny-1
        A=dom_ov(i1,i2).AE;
        mesh=dom_ov(i1,i2).mesh;
        free=dom_ov(i1,i2).free;
        free=horzcat(free,free+mesh.nv);
        AI=A(free,free);
        Igfree=dom_ov(i1,i2).Igfree;
        Igfree=horzcat(Igfree,Igfree+nvel);
        lu=u(Igfree);
        lAu=AI\lu;
        w(Igfree)=w(Igfree)+lAu;
    end
end