function [dom_ov,badmat]=localeigenvectorsRandomDBopt(x,dom_ov,Nx,Ny,lambdamax,maxrandom)

for i1=1:Nx-1
    for i2=1:Ny-1
        Af=dom_ov(i1,i2).A;
        Mf=dom_ov(i1,i2).Mass;
        M1f=dom_ov(i1,i2).Mass1; % matriz de masa con coefficiente=1
        
        Ig=dom_ov(i1,i2).Ig;
        B=dom_ov(i1,i2).B;
        bry=[B.vup;B.vleft;B.vdown;B.vright];
        It=1:size(Ig,2);
        G=unique(bry);
        G2=G;
        if (i1==1)
            G2=setdiff(G2,B.vleft);
        end
        if (i1==(Nx-1))
            G2=setdiff(G2,B.vright);
        end
        if (i2==1)
            G2=setdiff(G2,B.vdown);
        end
        if (i2==(Ny-1))
            G2=setdiff(G2,B.vup);
        end
        G2=setdiff(G,G2);
       
        Ifree=setdiff(It,G2);
        Afree=Af(Ifree,Ifree);
        Mfree=Mf(Ifree,Ifree);
        M1free=M1f(Ifree,Ifree);
        
        mf=sum(M1free,2); % en cada nodo vale la integral de la funcion base
        naux=length(Afree);
        U=sparse(naux,300);
        U(:,1)=1;
        for r=2:maxrandom % 6/11/16
            f=rand(naux,1);
            f=f-dot(mf,f);
            u=[Afree,mf;mf',0]\[f;0];
            U(:,r)=u(1:naux);
        end
        A=U'*Afree*U;
        M=U'*Mfree*U;
        
        A=0.5*A+0.5*A';
        M=0.5*M+0.5*M';
        
        [Qaux,D]=eig(full(A),full(M));
        [lambda,I]=sort(diag(D),'ascend');
        Qaux=Qaux(I,I);
        Q=U(:,I)*Qaux;
        
        lambda=lambda(1:lambdamax);
        nlambda=max(size(lambda));
        incre=diff(lambda);
        Iincre=incre<0.001;
        I=1:nlambda-1;
        Ig90=I(Iincre);
        
        if size(Ig90,2)>3
            Nbad=size(Ig90,2);
            if Nbad >= lambdamax
                Nbad = lambdamax;
            end
        else
            Nbad=3;
        end
        
        if i1==1||i1==Nx-1||i2==1||i2==Ny-1
            Nbad=1;
        end

        n=size(x,1);
        psi1=sparse(n,Nbad);
        psi1(Ig(Ifree),1:(Nbad+1))=Q(:,1:(Nbad+1));
        badmat(i1,i2)=Nbad;
        dom_ov(i1,i2).psibad=psi1;
        dom_ov(i1,i2).Nbad=Nbad;
        dom_ov(i1,i2).lambda=lambda;
    end
end
badmat=rot90(badmat,-1);
