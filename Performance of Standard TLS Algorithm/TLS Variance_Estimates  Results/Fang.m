
 function [X,cf,sigma2_hat]=Fang(A,y,py,px)

    %计算Q矩阵
    [n,m]=size(A);
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
         X0=X;

     end

    %  计算单位权方差估计（公式2.21的WTLS形式）
    B_final = [kron(X', I), -I];               % 收敛后的 B 矩阵
    r_final=pinv(B_final*Q*B_final')*(y-A*X)-zeros(n,1);
    vA_final=[QA zeros(n) Qxy;
        zeros(n,3*n) ]*B_final'*r_final-zeros(2*n,1);
    res = y - A * X;                           % 观测残差
    omega = res' / (B_final * Q * B_final') * res; % 等价于 v^T * P * v
    sigma2_hat = omega / (n - 2);               % 单位权方差估计

end
