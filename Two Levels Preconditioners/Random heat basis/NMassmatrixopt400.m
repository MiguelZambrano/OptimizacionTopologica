function [A,A1]=NMassmatrixopt400(M,v,~,mesh,penal,E0,Emin,xPhys)

nvel=mesh.nv; % velocity degrees of freedom==number of vertices.
A=sparse(nvel,nvel); 
A1=sparse(nvel,nvel);
L=(1/9)*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4];
%hx=dom(1,1).h.hx;
%hy=dom(1,1).h.hy;
%L=(hx/2)*(hy/2)*L;

for i=1:mesh.ne
    lcol=M(i,:); % subdomain indexes of bases.
    Vi=v(lcol,:);
    %IK=mean(Vi);
    IK=floor(400*mean(Vi))+1;    
    A(lcol,lcol)=A(lcol,lcol)+(Emin+xPhys(IK(1),IK(2)).^penal*(E0-Emin))*L;
    A1(lcol,lcol)=A1(lcol,lcol)+ones(size(xPhys(IK(1),IK(2))))*L;
end

%z=zeros(nvel,nvel); 
%AE=[A z;z A];
%AE1=[A1 z;z A1];

%AE=sparse(2*nvel,2*nvel); % grad*grad
% v=[v;v];
% 
% for i=1:mesh.ne
%     lcol=[M(i,:) M(i,:)+nvel]; % subdomain indexes of bases.
%     lAE=localMass(M(i,:),v,mu,k1perm); % compute local part of AE
%     size(lAE)
%     size(AE(lcol,lcol))
%     AE(lcol,lcol)=AE(lcol,lcol)+lAE;
% end
