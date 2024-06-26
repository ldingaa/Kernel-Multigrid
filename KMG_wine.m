
clear all
close all

%ker=@(x,y) exp(-pdist2(x',y'));
ker=@(x,y) (1+pdist2(x',y')*30).*exp(-pdist2(x',y')*30);
%ker=@(x,y) (1+pdist2(x',y')*10+100*pdist2(x',y').^2/3).*exp(-pdist2(x',y')*10);

load wine.mat
n=1599;
m=10;
T=100;
D=11;
lambda=1;
I=speye(n);

ind=randperm(1599,n);
X=X(:,ind)+randn(D,n)/1000;
Y=Y(ind);
ind_m=randperm(n,m);
X_H=X(:,ind_m);
%X=[X_H rand(D,n-m)];

S=[];
A_diag=[];
A_H_diag=[];
for d=1:D
%K{d}=ker(X(d,:),X(d,:));
A{d}=inv(ker(X(d,:),X(d,:)));
%A{d}(abs(A{d})<.001)=0;
%A{d}=sparse(A{d});
%K_H{d}=ker(X_H(d,:),X_H(d,:));
A_H{d}=inv(ker(X_H(d,:),X_H(d,:)));
%A_H{d}(abs(A_H{d})<.001)=0;
%A_H{d}=sparse(A_H{d});
I_hH{d}=ker(X(d,:),X_H(d,:))*A_H{d};
%I_hH{d}(abs(I_hH{d})<.001)=0;
%I_hH{d}=sparse(I_hH{d});

A_diag=blkdiag(A_diag,A{d});
S=[S;I];
U{d}=zeros(n,1);
A_H_diag=blkdiag(A_H_diag,A_H{d});

end



M=A_diag+lambda*S*S';
%Y=randn(n,1)*D;
U_true=M\(S*Y);

U=randn(n,D);

U_BB=randn(n,D);

for d1=1:D
    ind=(d1-1)*n+1;
    ind_end=d1*n;
    U_true_d(:,d1)=U_true(ind:ind_end);
    for d2=1:D
        ind_1=(d1-1)*m+1;
        ind_1_end=d1*m;
        ind_2=(d2-1)*m+1;
        ind_2_end=d2*m;
        Sigma_mm(ind_1:ind_1_end,ind_2:ind_2_end)=I_hH{d1}'*I_hH{d2};
    end
end
%M_HH=[A_H{1}+lambda*I_hH{1}'*I_hH{1} lambda*I_hH{1}'*I_hH{2};lambda*I_hH{2}'*I_hH{1} A_H{2}+lambda*I_hH{2}'*I_hH{2}];
%M_HH=A_H_diag+[lambda*I_hH{1}'*I_hH{1} lambda*I_hH{1}'*I_hH{2};lambda*I_hH{2}'*I_hH{1} lambda*I_hH{2}'*I_hH{2}];
M_HH=A_H_diag+lambda*Sigma_mm;
for t=1:T

    for d=1:D
        U_temp=U;
        U_temp(:,d)=[];
        U(:,d)=lambda*(A{d}+lambda*I)\(Y-sum(U_temp,2));

        U_temp=U_BB;
        U_temp(:,d)=[];
        U_BB(:,d)=lambda*(A{d}+lambda*I)\(Y-sum(U_temp,2));
    end
    
    %r=[Y;Y]-M*[U(:,1);U(:,2)];
    r=S*Y-M*reshape(U,[D*n,1]);
    r_m=zeros(D*m,1);
    for d=1:D
       ind_n=(d-1)*n+1;
       ind_n_end=d*n;
       ind_m=(d-1)*m+1;
       ind_m_end=d*m;
       r_m(ind_m:ind_m_end)=I_hH{d}'*r(ind_n:ind_n_end);
    end
    %r=[I_hH{1}'*r(1:n);I_hH{2}'*r(n+1:end)];
    r=r_m;
    r=(M_HH)\r;
    for d=1:D
       ind_m=(d-1)*m+1;
       ind_m_end=d*m;
       delta(:,d)=I_hH{d}*r(ind_m:ind_m_end);
       U(:,d)=U(:,d)+delta(:,d);

    end
    
    %U(:,1)=U(:,1)+delta(:,1);
    %U(:,2)=U(:,2)+delta(:,2);
 %}
    Error(t)=norm(U_true-reshape(U,[D*n,1]))/(D*n)^.5;
    Error_BB(t)=norm(U_true-reshape(U_BB,[D*n,1]))/(D*n)^.5;
    
    hold on
    plot(1:t,log(Error),'Color','b')
    plot(1:t,log(Error_BB),'Color','r')
    legend('KMG','Backfit')
    drawnow
end

