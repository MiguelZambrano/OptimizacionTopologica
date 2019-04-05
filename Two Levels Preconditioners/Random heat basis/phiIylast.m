function z=phiIylast(c,diamx,diamy,vx,vy)

x=(1- 2*abs(vx-c(1))/diamx);
y= min(1,(1+ 2*(vy-c(2))/diamy));
z=x.*y;
