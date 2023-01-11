function [svd_U,svd_D,svd_V,Z]=Initialization2(P,X,Y)
Triangel=Para_Tri2(P(:,1),P(:,2),P(:,3),X);
Z=[Y(P(:,2),:)-Y(P(:,1),:),Y(P(:,3),:)-Y(P(:,1),:)]';
Z([2,3],:)=Z([3,2],:);
Invtemp=1./(Z(4,:).*Z(1,:)-Z(2,:).*Z(3,:)).*Z;
Invtemp([2,3],:)=-Invtemp([2,3],:);
Invtemp([1,4],:)=Invtemp([4,1],:);
tempinv=reshape(Invtemp,2,2,[]);
Z=reshape(Z,2,2,[]);
% for i=1:nP
%    temp(:,:,i)=Z(:,:,i)\Triangel(:,:,i);
% end
% for i=1:nP
%    [svd_U(:,:,i),svd_D(:,:,i),svd_V(:,:,i)]=svd(temp(:,:,i));
% end
[svd_U,svd_D,svd_V]=projectJacobian2Rotation(Triangel,tempinv);
end

% function [x]=Get_Singular(A)
% a=A(1,1);
% b=A(1,2);
% c=A(2,1);
% d=A(2,2);
% 
% x(1)=0.5*(sqrt((b+c)^2+(a-d)^2)+sqrt((b-c)^2+(a+d)^2));
% x(2)=0.5*(sqrt((b+c)^2+(a-d)^2)-sqrt((b-c)^2+(a+d)^2));
% end
% function A=Para_Tri1(x1,x2,x3,Position)
%     OA=Position(x2,:)-Position(x1,:);
%     OB=Position(x3,:)-Position(x1,:);
%     theta = acos(dot(OA,OB)/(norm(OA)*norm(OB)));
%     A=[norm(OA)          ,0;
%       norm(OB)*cos(theta),norm(OB)*sin(theta)];
% end
function A=Para_Tri2(p1,p2,p3,X)
    OA=X(p2,:)-X(p1,:);
    OB=X(p3,:)-X(p1,:);
    normA=vecnorm(OA,2,2);
    normB=vecnorm(OB,2,2);
    theta = acos(dot(OA,OB,2)./(normA.*normB));
    A=reshape([normA,normB.*cos(theta),zeros(size(p1,1),1),normB.*sin(theta)]',2,2,[]);
end
