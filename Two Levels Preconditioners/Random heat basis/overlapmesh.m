function dom = overlapmesh(ax,bx,ay,by,Nx,Ny,nx,ny)
fprintf('Setting overlap mesh \n');
for i1 = 1:Nx-1
    for i2 = 1:Ny-1
        intx = [ax+(i1-1)*(bx-ax)/Nx,ax+(i1+1)*(bx-ax)/Nx];
        inty = [ay+(i2-1)*(by-ay)/Ny,ay+(i2+1)*(by-ay)/Ny];
        [M,B,v,mesh,h,Igfree,Ig,free] = createmeshAS2(nx,ny,intx,inty,Nx,Ny,i1,i2);
        dom(i1,i2).M = M;
        dom(i1,i2).B = B;
        dom(i1,i2).v = v;
        dom(i1,i2).mesh = mesh;
        dom(i1,i2).h = h;
        dom(i1,i2).Ig = Ig;
        dom(i1,i2).Igfree = Igfree;
        dom(i1,i2).c = [ax+(i1)*(bx-ax)/Nx,ay+(i2)*(by-ay)/Ny];
        dom(i1,i2).diam.x = (2)*(bx-ax)/Nx;
        dom(i1,i2).diam.y = (2)*(by-ay)/Ny;
        dom(i1,i2).free = free;
    end
end
