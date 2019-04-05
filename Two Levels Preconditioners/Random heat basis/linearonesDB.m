function dom_ov=linearonesDB(dom_ov,Nx,Ny,b)

for i1=1:Nx-1
    for i2=1:Ny-1
        xlin=b*0;
        %xlinx=b*0;
        %xliny=b*0;
        c=dom_ov(i1,i2).c;
        diamx=dom_ov(i1,i2).diam.x;
        diamy=dom_ov(i1,i2).diam.y;
        Ig=dom_ov(i1,i2).Ig;
        v=dom_ov(i1,i2).v;
        vx=v(:,1);
        vy=v(:,2);
        xlin(Ig)=phiI(c,diamx,diamy,vx,vy);
        %bdr=0;
        
        if i1==1
            xlin(Ig)=phiIx(c,diamx,diamy,vx,vy);
            %bdr=1;
        end
        if i1==Nx-1
            xlin(Ig)=phiIxlast(c,diamx,diamy,vx,vy);
            %bdr=1;
        end
        if i2==1
            xlin(Ig)=phiIy(c,diamx,diamy,vx,vy);
            %bdr=1;
        end
        if i2==Ny-1
            xlin(Ig)=phiIylast(c,diamx,diamy,vx,vy);
            %bdr=1;
        end
        if i1==1 && i2==1
            xlin(Ig)=phiIxy(c,diamx,diamy,vx,vy);
        end
        if i1==1 && i2==Ny-1
            xlin(Ig)=phiIxylast(c,diamx,diamy,vx,vy);
        end
        if i1==Nx-1 && i2==1
            xlin(Ig)=phiIxlasty(c,diamx,diamy,vx,vy);
        end
        if i1==Nx-1 && i2==Ny-1
            xlin(Ig)=phiIxlastylast(c,diamx,diamy,vx,vy);
        end
        
        %         if %bdr==1
        %         xlin=xlinx;
        %         end
        
        dom_ov(i1,i2).philin=xlin;
        %        plot_vector(xlin,dom,Nx,Ny);
        %        pause
        %        hold off
        rxg=b*0;
        ryg=b*0;
        rxg(Ig)=vy;
        ryg(Ig)=-vx;
        dom_ov(i1,i2).rx=rxg;
        dom_ov(i1,i2).ry=ryg;     
    end
end