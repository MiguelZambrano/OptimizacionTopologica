function [A,AE]=Nmatrixopt400(M,v,mesh,nu,penal,E0,Emin,xPhys)
%% Heat stiffness matrix
nvel=mesh.nv; % velocity degrees of freedom=number of vertices
A=sparse(nvel,nvel); % grad*grad   
%b=sparse(nvel,1); % right hand side part
L=[4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4];
for i=1:mesh.ne
    %%%%%%%%%%%%%%%%%%% GLOBAL NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lcol=M(i,:); % subdomain indexes of bases.
    %%%%%%%%%%%%%%%%%%% VELOCITY*VELOCITY STIFNESS %%%%%%%%%%%%%%%%%%%%
    Vi=v(lcol,:);
    %IK=mean(Vi);
    IK=floor(400*mean(Vi))+1;
    A(lcol,lcol)=A(lcol,lcol)+(Emin+xPhys(IK(1),IK(2)).^penal*(E0-Emin))*L;
    %%%%%%%%%%%%%%%%%%%%%%%% LOAD VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %lb=localb(M(i,:),v,mesh);
    %b(lcol,1)=b(lcol,1)+lb';
end

%% Elasticity stiffness matrix
AE=sparse(2*nvel,2*nvel); % grad*grad   
%be=sparse(2*nvel,1); % right hand side part
A11=[12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12=[-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11=[-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12=[ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
LEaux=1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
LE=LEaux([1,3,5,7,2,4,6,8],[1,3,5,7,2,4,6,8]);
for i=1:mesh.ne
    %%%%%%%%%%%%%%%%%%% GLOBAL NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lcol=[M(i,:) M(i,:)+nvel]; % subdomain indexes of bases.
    %%%%%%%%%%%%%%%%%%% VELOCITY*VELOCITY STIFNESS %%%%%%%%%%%%%%%%%%%%
    Vi=v(lcol(1:4),:);
    IK=floor(400*mean(Vi))+1;  
    AE(lcol,lcol)=AE(lcol,lcol)+(Emin+xPhys(IK(1),IK(2)).^penal*(E0-Emin))*LE;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% LOAD VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plb=localb(M(i,:),v,mesh);
    %lbe=[plb plb];
    %be(lcol,1)=be(lcol,1)+lbe';
end
