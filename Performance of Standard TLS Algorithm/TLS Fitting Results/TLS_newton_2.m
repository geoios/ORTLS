 function [X,ci] = TLS_newton_2(A,L,P)

    % 输入:
    %   A: 矩阵 [ai1, ai2, ..., ain]
    %   L: 列向量 [y1, y2, ..., yn]'
    %   P: [L|A]对应的权阵 [pi0, ai1, ..., aim]
    % 输出:
    %   X 的值

    [m,n]=size(P);
    PP=P(2:m,:);
    for i=1:n
        Pi{i}=diag(PP(:,i));
        P_i{i}=pinv(Pi{i});
    end
    
    %最小二乘求初值
    p0=P(1,:);
    P0=diag(p0);
     
    X0=pinv(A' *P0* A )*A'*P0*L;

    cita=1;
    ci=0;
   
    while cita>10^(-10)
        tic;
        ci=ci+1;
        v=L-A*X0;
        C1=0;
        for i=1:n
            w(i)=p0(i)/(1+p0(i)*X0'*P_i{i}*X0);
            e{i}=w(i)*v(i)*P_i{i}*X0;
            E(i,:)=e{i}';
            C1=C1+w(i)^2*v(i)^2*P_i{i};
        end
        W=diag(w);
        C2=-4*A'*W*E;
        C3=-A'*W*A;
        C4=-4*E'*W*E;
        F=(A+E)'*W*v;
        H=C1+C2+C3+C4;
        X=X0-pinv(H)*F;
        cita=norm(X-X0);
        X0=X;
        t=toc;
    end
    
end