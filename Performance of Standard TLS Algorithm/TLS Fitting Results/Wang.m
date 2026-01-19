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

% function [X, sigma2, Qxx, ci, e] = Wang(n, A, L, py, pA, is_random_col2)
% % 修改后的Partial EIV模型求解函数，支持第二列为随机列的情况
% % 输入:
% %   n: 观测数量
% %   A: 系数矩阵（包含所有列）
% %   L: 观测向量
% %   py: y观测值的权
% %   pA: x观测值的权（若第二列为随机，pA应为2n维向量）
% %   is_random_col2: 逻辑值，1表示第二列为随机列，0表示第二列为常数列
% % 输出:
% %   X: 参数估计值
% %   sigma2: 单位权方差估计值
% %   Qxx: 参数协因数阵
% %   ci: 迭代次数
% %   e: 残差向量
% 
% % 构建对角协因数阵（不考虑相关性）
% if is_random_col2
%     % 第二列为随机列，随机元素总数为2n
%     t = 2 * n;
%     Q = zeros(t + n, t + n);  % 维度为(3n, 3n)
%     
%     % 设置对角线元素
%     for i = 1:t
%         Q(i, i) = 1 / pA(i);
%     end
%     for i = 1:n
%         Q(t+i, t+i) = 1 / py(i);
%     end
% else
%     % 第二列为常数列，随机元素总数为n
%     t = n;
%     Q = zeros(2*n, 2*n);  % 维度为(2n, 2n)
%     
%     % 设置对角线元素
%     for i = 1:n
%         Q(i, i) = 1 / pA(i);
%         Q(n+i, n+i) = 1 / py(i);
%     end
% end
% 
% % 最小二乘初值
% Py = diag(py);
% X0 = pinv(A' * Py * A) * A' * Py * L;
% 
% % 初始化
% I = eye(n);
% C1 = eye(n + t);  % 转换矩阵
% 
% if is_random_col2
%     % 第二列为随机列：没有常数部分，h=0
%     h = zeros(n*2, 1);  % vec(A)的长度为2n
%     % B为单位矩阵，因为所有元素都是随机的
%     B = eye(2*n);
% else
%     % 第二列为常数列：h包含常数部分
%     h = [zeros(n, 1); ones(n, 1)];  % 常数部分：第一列为0，第二列为1
%     b = [1; 0];  % 随机元素只出现在第一列
%     B = kron(b, I);
% end
% 
% ci = 0;
% cita = 1;
% e = zeros(n + t, 1);
% 
% % 从A中提取随机元素向量a
% if is_random_col2
%     % 所有列都是随机的，a为vec(A)
%     a = reshape(A, 2*n, 1);
% else
%     % 只有第一列是随机的
%     a = A(:, 1);
% end
% 
% % 迭代求解
% while cita > 1e-10
%     ci = ci + 1;
%     
%     % 构造C2矩阵
%     C2 = [-kron(X0', I) * B, I] * C1;
%     
%     % 计算拉格朗日乘子
%     Q1 = C2 * Q * C2';
%     lambda = Q1 \ (L - kron(X0', I) * (h + B * a));
%     
%     % 计算残差
%     e = -Q * C2' * lambda;
%     
%     % 提取随机元素的残差
%     eA = e(1:t);
%     
%     % 构建误差矩阵EA
%     EA = reshape(B * eA, n, 2);
%     % 更新系数矩阵和观测向量
%     A0 = A + EA;
%     y1 = L + EA * X0;
%     
%     % 计算新参数估计
%     X = pinv(A0' / Q1 * A0) * A0' / Q1 * y1;
%     
%     % 检查收敛
%     cita = norm(X - X0);
%     X0 = X;
% end
% 
% 
% % 构造C2矩阵
% C2 = [-kron(X', I) * B, I] * C1;
% 
% % 计算拉格朗日乘子
% Q1 = C2 * Q * C2';
% lambda = Q1 \ (L - kron(X', I) * (h + B * a));
% 
% % 计算残差
% e = -Q * C2' * lambda;
% % 提取随机元素的残差
% eA = e(1:t);
% 
% % 构建误差矩阵EA
% EA = reshape(B * eA, n, 2);
% % 更新系数矩阵和观测向量
% A0 = A + EA;
% % 计算单位权方差
% % 自由度：观测值个数(n+t) - 参数个数(t+m)
% % 对于直线拟合：n=观测点数量，t=随机元素个数，m=参数个数=2
% if is_random_col2
%     r = (n + 2*n) - (2*n + 2);  % 当第二列为随机时，t=2n
% else
%     r = (n + n) - (n + 2);      % 当第二列为常数时，t=n
% end
% 
% W = inv(Q);  % 权阵
% sigma2 = (e' * W * e) / r;
% 
% % 计算参数协因数阵
% Qxx = inv(A0' / Q1 * A0);
% 
% end