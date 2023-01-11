% clear 'all  '
%
clear 'all  ';
%  mex -V -R2018a CFLAGS="$CFLAGS -std=c99 -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" judgement.cpp

k=100000;
A=rand(36,k);
B=rand(36,k);
C=zeros(k,1);
tic
for i=1:10
f1=judgement(A,B,C);
end
toc
tic
for i=1:10
f2=sum(abs(A-B))./sum(abs(B));
end
toc
find(abs(f1-f2')~=0)









% 
% %  mex -R2018a pre_treat.cpp
% %   mex  -R2018a my_modify.cpp -I.
% %   mex  -R2018a update.cpp
% %  mex -R2018a fulltosparseIC.c              
%   mex -R2018a  test2.cpp -I.
% n = 10;
% e =ones(n,1);
% A = spdiags([0.2*e,e*0,1*e,e*0,0.2*e] ,-2:2,n,n);
% A(1,2)=0.1;
% A(2,1)=0.1;
% r=dissect(A);
% A=A(r,r);
% k=(1:100);
% [count,h,parent,post,M] = symbfact(A,'lo','lower');
% M=my_modify(M);
% [level_r,level_i]=paralevel(parent);
% L=M;
% D=ones(1,n);
% [L,D]=test2(L,D,tril(A),M',level_r,level_i);
% full(L);
% 
% 
% function [level_r,level_i]=paralevel(p)
% n=numel(p);
% level=ones(n,1);
% 
% for i=1:n
%     if(p(i)~=0)
%         level(p(i))=max([level(p(i)),level(i)+1]);
%     end
% end
% N=max(level);
% c=zeros(N,1);
% level_i=zeros(n,1);
% for i=1:n
%     c(level(i))=c(level(i))+1;
% end
% level_r=zeros(N+1,1);
% level_r(1)=1;
% for i=2:N+1
%     level_r(i)=level_r(i-1)+c(i-1);
%     c(i-1)=level_r(i-1);
% end
% for i=1:n
%     level_i(c(level(i)))=i;
%     c(level(i))=c(level(i))+1;
% end
% end







