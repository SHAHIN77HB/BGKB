function [ c_b , c , q_l , p_lk , p_l1k ] = BGKBfun( A , B , l , k )
%% Block Golub-Kahan Bidiagonalization For AX=B
%% Github.com/shahin77hb
%Inputs : Matrix A
%         Right handside matrix B
%         l : stepsize of the process

%Output : Bidiagonal matrix C
%         Orthonormal matrices p and q
% In this function k is number of iterations in the main process of algorithms of
% reconstruction of an image
%% 

[ss s]=size(B);

[p(:,1:s) beta(:,1:s)]=qr(B,0);
w=A'*p(:,1:s);
[q(:,1:s) alpha(:,1:s)]=qr(w,0);
alpha(:,1:s)=alpha(:,1:s)';

for j=2:l+1
    z=A*q(:,((j-2)*s+1):(j-1)*s)-p(:,((j-2)*s+1):(j-1)*s)*alpha(:,(j-2)*s+1:(j-1)*s);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p(:,(j-1)*s+1:j*s) beta(:,(j-1)*s+1:j*s)]=qr(z,0);
    w=A'*p(:,(j-1)*s+1:j*s)-q(:,(j-2)*s+1:(j-1)*s)*beta(:,(j-1)*s+1:j*s);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Uncomment For Full Reorthogonalization
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [q(:,(j-1)*s+1:j*s) alpha(:,(j-1)*s+1:j*s)]=qr(w,0);
    alpha(:,(j-1)*s+1:j*s)=alpha(:,(j-1)*s+1:j*s)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p=GSO(p);
    q=GSO(q);
end


c1=alpha(:,1:s);
for i=2:l+1
    c1= blkdiag(c1,alpha(:,(i-1)*s+1:i*s));
end


c2=beta(:,1:s);
for i=2:l+1
    c2=blkdiag(c2,beta(:,(i-1)*s+1:i*s));
end

c1=c1(:,1:s*l);
c2=c2(:,s+1:(l+1)*s);
c=c1+c2;

q_l=q(:,1:l*k);
p_l1k=p(:,1:(l+1)*k);
p_lk=p(:,1:l*k);
c_b=c(1:(l+1)*k,1:l*k);
c=c(1:l*k,1:l*k);
