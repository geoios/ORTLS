function [X,cf]=Fang(x,y,py,px)
% function [X,cf,sigma2_hat, vA_final]=Fang(x,y,py,px)

    %计算Q矩阵
    [n,m]=size(x);
    q=zeros(n,n);
    Py=diag(py);
    PA=diag(px);
    QA=inv(PA);
    Qy=inv(Py);
    Qxy=q;
%     for i=1:n
%         Qxy(i,i)=rho(i)*sqrt(1/pA(i))*sqrt(1/py(i));
%     end
    Q=[QA,zeros(n),Qxy;zeros(n,3*n);Qxy,zeros(n),Qy];
    
    %最小二乘初值
    A=[x ones(n,1)];
    X0=pinv(A' *Py* A )*A'*Py*y;
    X0=X0-zeros(2,1);
    
    %求初值
    I=eye(n);
    
    cf=0;
    cita=1;
    
     while cita>10^(-10)
    
         cf=cf+1;
         B=[kron(X0',I),-I];
         r=pinv(B*Q*B')*(y-A*X0)-zeros(n,1);
         vA=[QA zeros(n) Qxy;
             zeros(n,3*n) ]*B'*r-zeros(2*n,1);
         X=pinv(A'*pinv(B*Q*B')*A)*(kron(eye(2),r')*vA+A'*pinv(B*Q*B')*y);
         cita=norm(X-X0);
         X=X-zeros(2,1);
         X0=X;

     end
%     %  计算单位权方差估计（公式2.21的WTLS形式）
%     B_final = [kron(X', I), -I];               % 收敛后的 B 矩阵
%     r_final=pinv(B_final*Q*B_final')*(y-A*X)-zeros(n,1);
%     vA_final=[QA zeros(n) Qxy;
%         zeros(n,3*n) ]*B_final'*r_final-zeros(2*n,1);
%     res = y - A * X;                           % 观测残差
%     omega = res' / (B_final * Q * B_final') * res; % 等价于 v^T * P * v
%     sigma2_hat = omega / (n - 2);               % 单位权方差估计

end
function [X,ci]=Wang(n,A1,L,py,pA)

%计算Q矩阵
q=zeros(n,n);
Py=diag(py);
PA=diag(pA);
QA=inv(PA);
Qy=inv(Py);
Qxy=q;
% for i=1:n
%     Qxy(i,i)=rho(i)*sqrt(1/pA(i))*sqrt(1/py(i));
% end
Q=[QA,Qxy;Qxy,Qy];


%最小二乘初值
A=[A1 ones(n,1)];
X0=pinv(A' *Py* A )*A'*Py*L;

%求初值
I=eye(n);
% Q=[I,zeros(n);
%     zeros(n),I];
C1=eye(2*n);
h=[zeros(n,1);ones(n,1)];
b=[1;0];
B=kron(b,I);

ci=0;
cita=1;
 while cita>10^(-10)
%      tic
     ci=ci+1;
     C2=[-kron(X0',I)*B,I]*C1;
     lambda=(C2*Q*C2')\(L-kron(X0',I)*(h+B*A1));
     e=-Q*C2'*lambda;
     eA=e(1:n);
     EA=reshape(B*eA,n,2);
     A0=A+EA;
     y1=L+EA*X0;
     X=pinv(A0'*pinv(C2*Q*C2')*A0)*A0'*pinv(C2*Q*C2')*y1;
     cita=norm(X-X0);
     X0=X;
%      t=toc;
 end

end