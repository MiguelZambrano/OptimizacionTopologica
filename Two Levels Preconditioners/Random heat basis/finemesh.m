function dom = finemesh(ax,bx,ay,by,Nx,Ny,nx,ny)
fprintf('Setting fine mesh \n');
hx = (bx-ax)/(Nx*nx);
hy = (by-ay)/(Ny*ny);
ov = 0;
for i1 = 1:Nx
    for i2 = 1:Ny
        intx = [ax+(i1-1)*(bx-ax)/Nx-ov*hx*(i1>1),ax+(i1)*(bx-ax)/Nx+ov*hx*(i1<Nx)];
        inty = [ay+(i2-1)*(by-ay)/Ny-ov*hy*(i2>1),ay+(i2)*(by-ay)/Ny+ov*hy*(i2<Nx)];
        [M,B,v,mesh,h,Igfree,Ig,free] = createmeshAS(nx,ny,intx,inty,Nx,Ny,i1,i2,ov);
        dom(i1,i2).M = M; %mesh
        dom(i1,i2).B = B; % boundary vertices and elements
        dom(i1,i2).v = v; % vertices
        dom(i1,i2).mesh = mesh; % mesh info
        dom(i1,i2).h = h; % h info
        dom(i1,i2).Ig = Ig; % global nodes
        dom(i1,i2).Igfree = Ig; % free nodes
        dom(i1,i2).free = 1:size(Ig,2); %local enumeration
    end
end
for i1 = 1
    for i2 = 1:Ny
        vleft = dom(i1,i2).B.vleft;
        Ig = dom(i1,i2).Ig;
        free = dom(i1,i2).free;
        free = setdiff(free,vleft);
        Igfree = Ig(free);
        dom(i1,i2).Igfree = Igfree;
        dom(i1,i2).free = free;
    end
end
clear vleft
for i1 = Nx
    for i2 = 1:Ny
        vright = dom(i1,i2).B.vright;
        Ig = dom(i1,i2).Ig;
        free = dom(i1,i2).free;
        free = setdiff(free,vright);
        Igfree = Ig(free);
        dom(i1,i2).Igfree = Igfree;
        dom(i1,i2).free = free;
    end
end
clear vright
for i1 = 1:Nx
    for i2 = 1
        vdown = dom(i1,i2).B.vdown;
        Ig = dom(i1,i2).Ig;
        free = dom(i1,i2).free;
        free = setdiff(free,vdown);
        Igfree = Ig(free);
        dom(i1,i2).Igfree = Igfree;
        dom(i1,i2).free = free;
    end
end
clear vdown
for i1 = 1:Nx
    for i2 = Ny
        vup = dom(i1,i2).B.vup;
        Ig = dom(i1,i2).Ig;
        free = dom(i1,i2).free;
        free = setdiff(free,vup);
        Igfree = Ig(free);
        dom(i1,i2).Igfree = Igfree;
        dom(i1,i2).free = free;
    end
end