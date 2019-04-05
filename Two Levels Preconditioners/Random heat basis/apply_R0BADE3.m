function x0=apply_R0BADE3(x,dom_ov,Nx,Ny)
dimC=0;
for i1=1:Nx-1
    for i2=1:Ny-1
        dimC=dimC+dom_ov(i1,i2).Nbad;
    end
end

x0=zeros(3*dimC,1);
I=0;

for i1=1:Nx-1
    for i2=1:Ny-1
        Nbad=dom_ov(i1,i2).Nbad;
        for l=1:Nbad
            I=I+1;         
            phi1=dom_ov(i1,i2).cb(l).phi;
            rx=dom_ov(i1,i2).rx;
            ry=dom_ov(i1,i2).ry;
            x0(I)=dot(x,[phi1;phi1*0]);
            x0(I+dimC)=dot(x,[0*phi1;phi1]);
            x0(I+2*dimC)=dot(x,[phi1.*rx;phi1.*ry]);
        end
    end
end

