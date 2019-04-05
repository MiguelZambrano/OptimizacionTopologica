function A0=coarse_matrixBADE3(dom,dom_ov,Nx,Ny,nvel)
dimC=0;
for i1=1:Nx-1
    for i2=1:Ny-1
        dimC=dimC+dom_ov(i1,i2).Nbad;
    end
end

A0=sparse(3*dimC,3*dimC);
I=0;

for i1=1:Nx-1
    for i2=1:Ny-1
        NbadI=dom_ov(i1,i2).Nbad;
        for l=1:NbadI
            I=I+1;
            rx=dom_ov(i1,i2).rx;
            ry=dom_ov(i1,i2).ry;
            phi1=dom_ov(i1,i2).cb(l).phi;

            Aphi1x=apply_AE([phi1; phi1*0],dom,Nx,Ny,nvel);
            Aphi1y=apply_AE([0*phi1; phi1],dom,Nx,Ny,nvel);
            Aphi1r=apply_AE([rx.*phi1; ry.*phi1],dom,Nx,Ny,nvel);
            J=0;
            for j1=1:Nx-1
                for j2=1:Ny-1
                    NbadJ=dom_ov(j1,j2).Nbad;
                    for m=1:NbadJ
                        phi2=dom_ov(j1,j2).cb(m).phi;
                        rx=dom_ov(j1,j2).rx;
                        ry=dom_ov(j1,j2).ry;
                        
                        phi2x=[phi2;phi2*0];
                        phi2y=[phi2*0;phi2];
                        phi2r=[rx.*phi2;ry.*phi2];
                        
                        J=J+1;
                        
                        A0(I,J)=sum(Aphi1x.*phi2x);
                        A0(I+dimC,J)=sum(Aphi1y.*phi2x);                         A0(I,J+dimC)=sum(Aphi1x.*phi2y);
                        A0(I+dimC,J+dimC)=sum(Aphi1y.*phi2y);
                        
                        A0(I+2*dimC,J)=sum(Aphi1r.*phi2x);
                        A0(I,J+2*dimC)=sum(Aphi1x.*phi2r);
                        
                        A0(I+2*dimC,J+dimC)=sum(Aphi1r.*phi2y);
                        A0(I+dimC,J+2*dimC)=sum(Aphi1y.*phi2r);
                        
                        A0(I+2*dimC,J+2*dimC)=sum(Aphi1r.*phi2r);
                    end
                end
            end
        end
    end
end
