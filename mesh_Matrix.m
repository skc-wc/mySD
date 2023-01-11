classdef mesh_Matrix < handle
    properties(SetAccess = protected)
        f1;  % matrix, nonzero pattern
        f2;
        P;
         hes_vec;
        hess=[];
        r;
        invr;
        oldhess=[];
        ldl_L=[];
        ldl_D;
        X;
        ETREE;
        ldl_LMODEL=[];
        ldl_LMODEL_T=[];
        Hess_Model=[];
        level;
        levelr;
        leveli;
    end
    methods
        function obj=mesh_Matrix(P,X)
             %% obj.matrix
            %% obj.matrix
            nX=size(X,1);
            obj.P=P;
            obj.X=X;
            %% obj.hess
            Hvec=ones(6,6,size(obj.P,1));
            obj.hes_vec= Hvec;
            new_hess(obj,Hvec);
            obj.r=dissect(obj.hess);
            obj.invr=zeros(size(obj.r));
            for i=1:numel(obj.r)
                obj.invr(obj.r(i))=i;
            end 
            obj.hess=obj.hess(obj.r,obj.r);
            [~,~,obj.ETREE,~,obj.ldl_LMODEL] = symbfact(obj.hess,'lo','lower');
            obj.ldl_D=ones(1,size(obj.hess,1));
            obj.ldl_LMODEL=my_modify(obj.ldl_LMODEL);
           
            
            H_Model=my_modify(obj.hess);
            H_Model=H_Model(obj.invr,obj.invr);
            Order=[obj.P,obj.P+nX]';
            H_i=repelem(Order,6,1) ;
            H_j=repmat(Order,6,1);
            my_order=sub2ind(size(H_Model),H_i(:),H_j(:));
            obj.Hess_Model=int32(full(H_Model(my_order)));
            obj.ldl_L=obj.ldl_LMODEL;
            obj.ldl_LMODEL_T=obj.ldl_LMODEL';
             getplevel(obj,obj.ETREE);
              [obj.levelr,obj.leveli]=option_lv(obj,(1:size(obj.hess)));
        end

        function Hupdate(obj,newHvec,transtri)
%                         if (isempty(transtri))
%                             obj.oldhess=[];
%             %                 obj.ldl_D=ones(1,size(obj.hess,1));
%                             obj.hes_vec=newHvec;
%                         else
%                             obj.oldhess=obj.hess;
%                             obj.hes_vec(:,:,transtri)=newHvec(:,:,transtri);
%                         end
%                         obj.hess=H_update(obj.hess,obj.Hess_Model,obj.hes_vec);
                        if (isempty(transtri))
                            obj.oldhess=[];
            %                 obj.ldl_D=ones(1,size(obj.hess,1));
                            obj.hes_vec=newHvec;
                            obj.hess=H_update2(obj.hess,obj.Hess_Model,obj.hes_vec,newHvec);
                             obj.hess(1,1)=obj.hess(1,1)+1e-5;
                        else
                            obj.oldhess=obj.hess;
                            obj.hess=H_update2(obj.hess,obj.Hess_Model,obj.hes_vec,newHvec,int32(transtri));
                            obj.hes_vec(:,:,transtri)=newHvec(:,:,transtri); 
%                             obj.hess=H_update(obj.hess,obj.Hess_Model,obj.hes_vec);
                        end
                        
        end

        function new_hess(obj,hes_vec)
            nX=size(obj.X,1);
            Order=[obj.P,obj.P+nX]';
            H_i=repelem(Order,6,1) ;
            H_j=repmat(Order,6,1); 
            obj.hess=sparse(H_i(:),H_j(:),hes_vec(:),nX*2,nX*2)+spdiags(0.1*ones(2*nX,1),0,nX*2,nX*2);
        end

        function getplevel(obj,p)
            n=numel(p);
            obj.level=ones(n,1);
            for i=1:n
                if(p(i)~=0)
                    obj.level(p(i))=max([obj.level(p(i)),obj.level(i)+1]);
                end
            end
        end

        function step=solveE(obj,E, option_tri)  
            tt=tic;
            if(isempty(obj.oldhess))
                lr=obj.levelr;
                li=obj.leveli;
                obj.ldl_L=obj.ldl_LMODEL;
                obj.ldl_D=ones(1,size(obj.hess,1));
            else

%                 [order,~,~] = find(obj.oldhess-obj.hess);
%                 order=unique(order);
                temp=unique([obj.P(option_tri,:)]);
                temp=[temp;temp+size(obj.X,1)];
                temp2=unique(obj.invr(temp));
                order=temp2(:);
                p=obj.ETREE;
                test_I=zeros(size(p,2),1);
                for k=1:numel(order)
                    i=order(k);
                    while(test_I(i)==0)
                        test_I(i)=1;
                        if(p(i)~=0)
                            i=p(i);
                        else
                            break;
                        end
                    end
                end
                I2=(find(test_I==1))';
                [lr,li]=option_lv(obj,I2);
            end
 %            [obj.ldl_L,obj.ldl_D]=test2(obj.ldl_LMODEL,ones(1,size(obj.hess,1)),tril(obj.hess),obj.ldl_LMODEL_T,obj.level_r,obj.level_i);
%     obj.hess=tril(obj.hess);  
%     obj.oldhess=tril(obj.oldhess);
             
            [obj.ldl_L,obj.ldl_D]=SparseUpdate(obj.ldl_L,obj.ldl_D,tril(obj.hess),obj.ldl_LMODEL_T,lr,li);
            
%              fprintf("row=%f,time1=%f",numel(li),toc(t2));
            %% LTDLx=y
            step=(obj.ldl_L')\((obj.ldl_L\E(obj.r))./(obj.ldl_D'));
            step=step(obj.invr);
            fprintf("Node2=%d,ldltime=%f",numel(li),toc(tt));
        end

        function [level_r,level_i]=option_lv(obj,order)
            N=max(obj.level);
            n=numel(obj.level);
            c=zeros(N,1);
            for i=1:numel(order)
                c(obj.level(order(i)))=c(obj.level(order(i)))+1;
            end
            level_i=zeros(numel(order),1);
            level_r=zeros(N+1,1);
            level_r(1)=1;
            for i=2:N+1
                level_r(i)=level_r(i-1)+c(i-1);
                c(i-1)=level_r(i-1);
            end
            for j=1:numel(order)
                i=order(j);
                level_i(c(obj.level(i)))=i-1;
                c(obj.level(i))=c(obj.level(i))+1;
            end
            level_r = int32(level_r);
            level_i = int32(level_i);
        end

%         function option_tri=option_tri_num(obj,Y)
%             old_J=obj.Jaco;
%             updateY(obj,Y);
%             new_J=obj.Jaco;
% 
%             delta=squeeze(sum((old_J-new_J).^2,[1,2])); 
%             normF=squeeze(sum((old_J).^2,[1,2]));
%             my_temp=delta./normF;
%             max(my_temp);
%             option_tri=find( my_temp>0.1);
% %             size(obj.T,1)-numel(option_tri);
%         end
        
%         function hess=gethess(obj)
%             hess=obj.hess(obj.invr,obj.invr);
%         end

%         function updateY(obj,new)
%             Complex_Y=complex(new(:,1),new(:,2));
%             vector_Y=[zeros(size(Complex_Y(obj.T(:,1)))),Complex_Y(obj.T(:,2))-Complex_Y(obj.T(:,1)),Complex_Y(obj.T(:,3))-Complex_Y(obj.T(:,1))].';
%              Compute_Y=reshape(vector_Y,3,size(obj.T,1));
%             obj.f1=reshape(sum(obj.D.*Compute_Y),[],1);
%             obj.f2=reshape(sum(conj(obj.D).*Compute_Y),[],1);
%              obj.Jaco=[real(obj.f1+obj.f2);imag(-obj.f1+obj.f2); imag(obj.f1+obj.f2);real(obj.f1-obj.f2)];
%         end

    end
end