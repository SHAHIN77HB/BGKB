%% Block Golub-Kahan Bidiagonalization For AX=B
%% Github.com/shahin77hb
%Inputs : Matrix A
%         Right handside matrix B
%         l : stepsize of the process

%Output : Bidiagonal matrix C
%         Orthonormal matrices p and q



%%
clear
clc
a=input('A=');
b=input('b=');
l=input('stepsize=');
[ss s]=size(b);

%%

[p(:,1:s) beta(:,1:s)]=qr(b,0);
w=a'*p(:,1:s);
[q(:,1:s) alpha(:,1:s)]=qr(w,0);
alpha(:,1:s)=alpha(:,1:s)';

for j=2:l+1
    z=a*q(:,((j-2)*s+1):(j-1)*s)-p(:,((j-2)*s+1):(j-1)*s)*alpha(:,(j-2)*s+1:(j-1)*s);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Uncomment For Full Reorthogonalization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for xx=1:s+1
%     z=z-p(:,((j-2)*s+1):(j-1)*s)*(p(:,((j-2)*s+1):(j-1)*s)'*z);
%     %z=z-p(:,((j-2)*s+1):(j-1)*s)*(p(:,((j-2)*s+1):(j-1)*s)'*z);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p(:,(j-1)*s+1:j*s) beta(:,(j-1)*s+1:j*s)]=qr(z,0);
    w=a'*p(:,(j-1)*s+1:j*s)-q(:,(j-2)*s+1:(j-1)*s)*beta(:,(j-1)*s+1:j*s);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Uncomment For Full Reorthogonalization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for zz=1:s+1
%     w=w-q(:,((j-2)*s+1):(j-1)*s)*(q(:,((j-2)*s+1):(j-1)*s)'*w);
%     %w=w-q(:,((j-2)*s+1):(j-1)*s)*(q(:,((j-2)*s+1):(j-1)*s)'*w);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%
disp('c=')
disp(c)
disp('p=')
disp(p)
disp('q=')
disp(q)
