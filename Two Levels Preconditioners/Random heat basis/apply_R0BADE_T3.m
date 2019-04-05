function x=apply_R0BADE_T3(x0,dom_ov,Nx,Ny)

x=0*(dom_ov(1,1).philin);
x=[x;x];
dimC=length(x0)/3;
%x0=sparse((Nx-1)*(Ny-1),1);
I=0;
for i1=1:Nx-1
    for i2=1:Ny-1
        Nbad=dom_ov(i1,i2).Nbad;
        for l=1:Nbad
            I=I+1;
            phi1=dom_ov(i1,i2).cb(l).phi;
            rx=dom_ov(i1,i2).rx;
            ry=dom_ov(i1,i2).ry;
            x=x+x0(I)*[phi1;phi1*0]+x0(I+dimC)*[phi1*0;phi1]+x0(I+2*dimC)*[phi1.*rx;phi1.*ry];
        end
    end
end
