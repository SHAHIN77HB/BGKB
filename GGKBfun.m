function [ c_b , c , v_lk , u_lk , u_lk1 ] = GGKBfun( a , b , l  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,k]=size(b);

sigma(1)=norm(b,'fro');
u(:,1:k)=b/sigma(1);
v(:,1:k)=zeros(size(u));

for j=1:l
    %4
    vv=a'*u(:,(j-1)*k+1:j*k)-sigma(j)*v(:,(j-1)*k+1:j*k);
    rho(j)=norm(vv,'fro');
    v(:,j*k+1:(j+1)*k)=vv/rho(j);
    %5
    uu=a*v(:,j*k+1:(j+1)*k)-rho(j)*u(:,(j-1)*k+1:j*k);
    sigma(j+1)=norm(uu,'fro');
    u(:,j*k+1:(j+1)*k)=uu/sigma(j+1);
end
    v=v(:,k+1:end);
%     u=GSO(u);
%     v=GSO(v);



c1=[diag(rho(1:l));zeros(1,l)];
c2=[zeros(1,l);diag(sigma(2:l+1))];
c=c1+c2;

u_lk1=u(:,1:(l+1)*k);
u_lk=u(:,1:l*k);
v_lk=v(:,1:l*k);
c_b=c;
c=c(1:l,1:l);