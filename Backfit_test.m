clear all
close all

ker=@(x,y) exp(-pdist2(x',y'));
%ker=@(x,y) (1+pdist2(x',y')*5).*exp(-pdist2(x',y')*5);
%ker=@(x,y) (1+pdist2(x',y')*20+400*pdist2(x',y').^2/3).*exp(-pdist2(x',y')*20);


n=500;
m=ceil(n.^.5);
D=50;
lambda=1;
I=eye(n);
I_m=eye(m);
gamma=1;

%U=rand(D,m);
%X=[U rand(D,n-m)];
X=lhsdesign(n,D)';

S=[];
S_m=[];
P_U=[];
G=zeros(n,1);
for d=1:D
    P{d}=inv(ker(X(d,:),X(d,:)));
    G=G+mvnrnd(zeros(n,1),ker(X(d,:),X(d,:)))'/D;

    %P_U=blkdiag(P_U,inv(ker(U(d,:),U(d,:))));
    %S=[S;I];
    %I_mn{d}=(ker(U(d,:),U(d,:)))\ker(U(d,:),X(d,:));
    %I_nm{d}=ker(X(d,:),U(d,:))/ker(U(d,:),U(d,:));
    %S_m=[S_m;I_m];
end
%Sigma_m=gamma*S_m*S_m';

for d=1:D
    %V(:,d)=ones(n,1);
    V(:,d)=G+randn(n,1);
end
%V(:,1)=zeros(n,1);

for t=1:100

    for d=1:D
        V_temp=V;
        V_temp(:,d)=[];
        V(:,d)=-lambda*(P{d}+lambda*I)\(sum(V_temp,2));
        %R(:,d)=P{d}*V(:,d);
    end
    %{
    V_sum=sum(V,2);
    R_m=[];
    for i=1:D
        R(:,d)=-(R(:,d)+lambda*V_sum)/lambda;
        R_m=[R_m;I_mn{d}*R(:,d)];
    end
    delta=lambda*(P_U+Sigma_m)\R_m;
    for d=1:D
        V(:,d)=V(:,d)+I_nm{d}*delta((1:m)+(d-1)*m);
    end

    %}

    Y(t)=(norm(V));
    plot(1:t,Y)
    xlabel('iterations')
    ylabel('error')
    drawnow
end



