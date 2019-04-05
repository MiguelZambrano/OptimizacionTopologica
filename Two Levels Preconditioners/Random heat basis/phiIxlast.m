function z=phiIxlast(c,diamx,diamy,vx,vy)

x=min(vx*0+1,(1+2*(vx-c(1))/diamx));
y=(1- 2*abs(vy-c(2))/diamy);
z=x.*y;
