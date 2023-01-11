clear all;
[X,P]=readObj('bozbezbozzel100K_cut');
%% [X,P]=readObj('camelhead');
B = findBoundary(X, P);
nP=size(P,1);
nX=size(X,1);
%% Tutte_参数化
L = laplacian(X, P, 'uniform');
I = setdiff(1:nX, B);
z = zeros(nX,1);
z(B) = exp(2i*pi*(1:numel(B))'/numel(B));
z(I) = -L(I,I)\(L(I,B)*z(B));
Y=[real(z),imag(z)];
%% plot mesh

%% 初始化
n=1;
My_Mesh=mesh_Matrix(P,X);
Order=[P,P+nX]';
mexoption = struct('energy_type',0,  'hessian_projection', 1, 'energy_param',1, 'verbose', 0);
K=0;
%% Progressive Parameterizations
t2=tic;
qwe=1;
end_qwe=10000;
iter=10000;
while(qwe<=end_qwe)
    [svd_U,svd_D,svd_V,Z]=Initialization2(P,X,Y);
    svd_VT=pagetranspose(svd_V);
    ZU=pagemtimes(Z,svd_U);
    delta1=squeeze(svd_D(1,1,:));
    delta2=squeeze(svd_D(2,2,:));
    K=10000;
    delta_t1=log10((K+sqrt(K^2-4))/2)./log10(delta1(delta1>1))/2;
    delta_t2=log10((K-sqrt(K^2-4))/2)./log10(delta1(delta1<1))/2;
    delta_t3=log10((K+sqrt(K^2-4))/2)./log10(delta2(delta2>1))/2;
    delta_t4=log10((K-sqrt(K^2-4))/2)./log10(delta2(delta2<1))/2;
    delta_t4(delta_t4<=0)=[];
    delta_t2(delta_t2<=0)=[];
    tempt=min([delta_t1;delta_t2;delta_t3;delta_t4]);
    K=max(delta1.^(2*tempt)+delta1.^(-2*tempt)+delta2.^(2*tempt)+delta2.^(-2*tempt));
    if(tempt>=0.8)
        K=2;
        tempt=1;
        end_qwe=qwe;
    end
    tempt
    tempS=reshape([delta1';delta2'].^(2*tempt),1,2,[]);
    qwe=qwe+1;
    Tri=pagemtimes(ZU.*tempS,svd_VT);
    %%
    z1=zeros(1,1,nP);
    z2=complex(Tri(1,1,:),Tri(1,2,:));
    z3=complex(Tri(2,1,:),Tri(2,2,:));
    cv_e=conj([z2-z3,z3,-z2]);
    v_z=[z1,z2,z3];
    d=pagemtimes(cv_e,pagetranspose(v_z));
    D=reshape(cv_e./d,3,[]);
    %% iter
    D2t=conj(D);
    D=D.';
    ComplexY=complex(Y(:,1),Y(:,2));
    CYP= ComplexY(P);
    fz=reshape(sum(D.*CYP,2),[],1);
    gz=reshape(sum(conj(D).*CYP,2),[],1);
    NewSDEngery=SD_engery(fz,gz);
    for iter=1:1000
        option_triangle=[];
        for count=1:2
            %% Hessian and Gradient
            if(count==2)
                oldhs=hs;
            end
            tt=tic;
            OldSDEngery=NewSDEngery;
            [~, g, hs] = meshIsometricEnergyC(fz(:), gz(:), D2t, ones(nP,1), mexoption);
            new_vecH=reshape(hs,6,6,[]);
            %         [new_vecH,g]=Get_Step(My_Mesh,nX);
            LaplaceE=accumarray(Order(:),g(:));
            %% option update
            option_triangle=[];
            if(count~=1)
                F=judgement(hs,oldhs,zeros(nP,1));
                elip=0.2;
                option_triangle=find( F>elip);
                if(numel(option_triangle)<100)
                    option_triangle=find( F>0.5*max(F,[],'all'));
                end
                if(numel(option_triangle)>0.25*nP)
                    option_triangle=[];
                end
            end
            Hupdate(My_Mesh,new_vecH, option_triangle);
            %% H\E
            temp_Step=solveE(My_Mesh,LaplaceE, option_triangle);
            Step=reshape(temp_Step,[],2);

            %% Gradient
            ComplexE=complex(Step(:,1),Step(:,2));
            CEP= ComplexE(P);
            dfz=sum(D.*CEP,2);
            dgz=sum(conj(D).*CEP,2);
            %% line search energy decreasing
            alpha= min( [maxtForPositiveArea( fz, gz, dfz, dgz )*0.95, 1] );
            c=0.1;
            E=sum(OldSDEngery);
            NewSDEngery=SD_engery(fz-dfz*alpha,gz-dgz*alpha);
            max_norm=norm(Step);
            while(sum(NewSDEngery)>E-c*alpha*(LaplaceE.')*Step(:) && max_norm*alpha>1e-12)
                alpha=alpha*0.5;
                NewSDEngery=SD_engery(fz-dfz*alpha,gz-dgz*alpha);
            end
            fz=fz-dfz*alpha;
            gz=gz-dgz*alpha;
            CYP=CYP-CEP*alpha;
            Y=Y-alpha*Step;
            n=n+1;
            %% loop end judgement
            fprintf("count=%d,slovetime=%f",count,toc(tt));
            fprintf(",N=%d,SDEngery=%f.\n",n,sum(NewSDEngery));
            if( sum(NewSDEngery)<=nP*10 &&tempt~=1)
                break;
            end
        end
        if( sum(NewSDEngery)<=nP*10 &&tempt~=1)
            break;
        end
    end
end

    function E=SD_engery(fz,gz)
    Nf=real(fz(:)).^2+imag(fz(:)).^2;
    Ng=real(gz(:)).^2+imag(gz(:)).^2;
    E=(Nf+Ng).*(1+1./((Nf-Ng).^2));
    end



