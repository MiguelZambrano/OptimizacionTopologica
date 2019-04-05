function z=phiIxlastylast(c,diamx,diamy,vx,vy)

x=min(vx*0+1,(1+2*(vx-c(1))/diamx));
y= min(1,(1+ 2*(vy-c(2))/diamy));
z=x.*y;
