
clear all
close all

w=1;
W=10000;
perturb=0.0001;
%ker=@(x,y) exp(-pdist2(x',y')*w);
ker=@(x,y) (1+pdist2(x',y')*w).*exp(-pdist2(x',y')*w);
%ker=@(x,y) (1+pdist2(x',y')*10+100*pdist2(x',y').^2/3).*exp(-pdist2(x',y')*10);

load susy.mat
n=3000000;
m=20;
T=100;
D=8;
lambda=100;
I=speye(n);
I_Dn=speye(D*n);

ind=randperm(5000000,n);
X=X(:,ind);
%X=W*X;
for d=1:D
    X_temp=sort(rand(1,n));
    [~,ind]=sort(X(d,:));
    X(d,:)=X_temp(ind)*W;
end

Y=Y(ind);
ind_m=randperm(n,m);
X_H=X(:,ind_m);
%X=[X_H rand(D,n-m)];

S=[];
A_diag=[];
invK_1_2_diag=[];
invK_H_diag=[];
for d=1:D
%K{d}=ker(X(d,:),X(d,:));

A{d}=sparse_A(X(d,:))*n^.5;
invK_1_2{d}=sparse_invK(X(d,:),w)/n^.25;
invK_1_2_square{d}=invK_1_2{d}*invK_1_2{d};
%K_H{d}=ker(X_H(d,:),X_H(d,:));
invK_H{d}=inv(ker(X_H(d,:),X_H(d,:)));
%invK_H{d}=sparse_invK(X_H(d,:),w);
%invK_H{d}(abs(invK_H{d})<.0000000001)=0;
%invK_H{d}=sparse(invK_H{d});
I_hH{d}=ker(X(d,:),X_H(d,:))*invK_H{d};


A_diag=blkdiag(A_diag,A{d});
invK_1_2_diag=blkdiag(invK_1_2_diag,invK_1_2{d});
S=[S;I];
U{d}=zeros(n,1);
invK_H_diag=blkdiag(invK_H_diag,invK_H{d});

end



%M=A_diag/lambda+(S*S');
%Y=randn(n,1)*D;
%U_true=M\(S*Y);
Y_S=S*Y;

U=rand(n,D);

U_BB=randn(n,D);

for d1=1:D
    ind=(d1-1)*n+1;
    ind_end=d1*n;
    %U_true_d(:,d1)=U_true(ind:ind_end);
    for d2=1:D
        ind_1=(d1-1)*m+1;
        ind_1_end=d1*m;
        ind_2=(d2-1)*m+1;
        ind_2_end=d2*m;
        Sigma_mm(ind_1:ind_1_end,ind_2:ind_2_end)=I_hH{d1}'*I_hH{d2};
    end
end
%M_HH=[invK_H{1}+lambda*I_hH{1}'*I_hH{1} lambda*I_hH{1}'*I_hH{2};lambda*I_hH{2}'*I_hH{1} invK_H{2}+lambda*I_hH{2}'*I_hH{2}];
%M_HH=invK_H_diag+[lambda*I_hH{1}'*I_hH{1} lambda*I_hH{1}'*I_hH{2};lambda*I_hH{2}'*I_hH{1} lambda*I_hH{2}'*I_hH{2}];
M_HH=invK_H_diag/lambda+Sigma_mm;
for t=1:T

    for d=1:D
        U_temp=U;
        U_temp(:,d)=[];
        Y_temp=Y-sum(U_temp,2);
        U(:,d)=invK_1_2{d}*Y_temp;
        U(:,d)=(A{d}*lambda+invK_1_2_square{d}+perturb*I)\U(:,d);
        U(:,d)=invK_1_2{d}*U(:,d);
        U(:,d)=Y_temp-U(:,d);
        %prod(isnan(U(:,d)))

        U_temp=U_BB;
        U_temp(:,d)=[];
        Y_temp=Y-sum(U_temp,2);
        U_BB(:,d)=invK_1_2{d}*Y_temp;
        U_BB(:,d)=(A{d}*lambda+invK_1_2_square{d}+perturb*I)\U_BB(:,d);
        U_BB(:,d)=invK_1_2{d}*U_BB(:,d);
        U_BB(:,d)=Y_temp-U_BB(:,d);
    end
    
    %r=[Y;Y]-M*[U(:,1);U(:,2)];
    r1_temp=S*S'*reshape(U,[D*n,1]);
    r2_temp=invK_1_2_diag*reshape(U,[D*n,1]);
    r2_temp=(A_diag*lambda+perturb*I_Dn)\r2_temp;
    r2_temp=invK_1_2_diag*r2_temp;
    r=Y_S-r1_temp-r2_temp;
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
    Y_BB_temp1=S*S'*reshape(U_BB,[D*n,1]);
    Y_BB_temp2=invK_1_2_diag*reshape(U_BB,[D*n,1]);
    Y_BB_temp2=(A_diag*lambda)\Y_BB_temp2;
    Y_BB_temp2=invK_1_2_diag*Y_BB_temp2;

    Y_temp1=S*S'*reshape(U,[D*n,1]);
    Y_temp2=invK_1_2_diag*reshape(U,[D*n,1]);
    Y_temp2=(A_diag*lambda)\Y_temp2;
    Y_temp2=invK_1_2_diag*Y_temp2;
    Error(t)=norm(Y_S-Y_temp2-Y_temp1)/(D*n)^.5;
    Error_BB(t)=norm(Y_S-Y_BB_temp2-Y_BB_temp1)/(D*n)^.5;
    Error_BB
    Error
    
hold on
    plot(1:t,log(Error),'Color','b')
    plot(1:t,log(Error_BB),'Color','r')
    legend('KMG','Backfit')
    drawnow
 %}
 t
end

