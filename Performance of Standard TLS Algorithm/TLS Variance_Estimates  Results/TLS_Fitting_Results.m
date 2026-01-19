
clc;clear;

dd = importdata('PTLS-data.txt'); %读取数据
n = size(dd,1);
m = 2;
k = m + 1;
A_obs = [dd(:,1),ones(n,1)];
L_obs = dd(:,3);



%% ============== 牛顿法 ==============
tic;
%计算权阵
P = [dd(:,4),dd(:,2),10^15*ones(n,1)]';
PP = diag(P(:));
[X_newton,P_v,cn] = TLS_XG_newton3(A_obs, L_obs, PP);
v_newton = L_obs - A_obs * X_newton;
% 单位权方差
sit_n = (v_newton' * P_v * v_newton) / (n - m);
sit_newton_sqrt = sqrt(sit_n);
tn = toc;



%% ============== Schaffrin法 ==============
%% Schaffrin (2015)算法实现 - 包含残差计算
% 构建观测向量的协因数矩阵 Q_y
tic;
p_y = dd(:, 4);  % L的权值
p_A = dd(:, 2);  % A第一列的权值
Q_y = diag(1 ./ p_y);  % n×n 对角矩阵

% 构建设计矩阵A的协因数矩阵 Q_A
p_A_vec = [1 ./ p_A; 1e-15 * ones(n, 1)];  % 前n个是x的协因数，后n个是常数项的协因数
Q_A = diag(p_A_vec);  % 2n×2n 对角矩阵

% 交叉协方差矩阵（假设为零）
Q_yA = zeros(n, 2*n);


max_iter = 100;
tol = 1e-10;
[X_s, sigma_S, e_y_Schaffrin, e_A_Schaffrin, iterations, converged, cj] = schaffrin_line_fit_with_residuals(A_obs, L_obs, Q_y, Q_A, Q_yA, max_iter, tol);
sit_schaffrin_sqrt = sqrt(sigma_S);
t_schaffrin = toc;
%% ============== Wang法 ==============
tic
P = [dd(:,4),dd(:,2),10^15*ones(n,1)]';
py = P(1,:);

pA = P(2,:);
[X_W, sigma_W, Qxx, ci] = Wang_sigma(n, A_obs, L_obs, py, pA,0);
sit_Wang_sqrt = sqrt(sigma_W);
t_wang = toc;

%% ============== Xu法 ==============
tic;
P = [dd(:,4),dd(:,2),10^15*ones(n,1)]';
py = P(1,:);
Py = diag(py);
[X_xu, a_hat, c_Xu, conv_flag, S0, A0, PA] = Xu_m(A_obs, L_obs, P);


a_hat0 = A_obs(:,1);
r_a = a_hat0(:) - a_hat;            % 随机元素的残差（使用观测值a）
r_y = L_obs - A0 * X_xu;       % 观测向量的残差

% 计算单位权方差
sigma_Xu = (r_a' * PA * r_a + r_y' * Py * r_y) / (n - m);
sit_Xu_sqrt = sqrt(sigma_Xu);
t_xu = toc;



%% ============== Fang法 ==============
tic
py = dd(:,4);
pA = dd(:,2);
[X,cf,sigma_Fang] = Fang(A_obs,L_obs,py,pA);
sit_Fang_sqrt = sqrt(sigma_Fang);
t_Fang = toc;



%% ============== jazaeri法 ==============
% 初始化权矩阵
tic
py = dd(:,4);
pA = dd(:,2);
Q_y = diag(1./py);
Q_A = zeros(m*n, m*n);
Q_A(1:n, 1:n) = diag(1./pA);
Q_A(n+1:m*n, n+1:m*n) = zeros(n);

% 从加权最小二乘开始
P0 = diag(py);
x_jaz = inv(A_obs' * P0 * A_obs) * (A_obs' * P0 * L_obs);
epsilon = 1e-10;

% 开始计时

for iter = 1:20
    % 估计残差
    e_hat = L_obs - A_obs * x_jaz;
    
    % 计算 Q_y_tilde = Q_y + (x^T ⊗ I_m) Q_A (x ⊗ I_m)
    x_kron_T = kron(x_jaz', eye(n));  % m × mn
    x_kron = kron(x_jaz, eye(n));      % mn × m
    Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
    Q_y_tilde_inv = inv(Q_y_tilde);
    
    % 计算 E_A_hat = -vec^{-1}(Q_A (x ⊗ I_m) Q_y_tilde^{-1} e_hat)
    vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * e_hat;
    E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)]; 
    
    % 预测值
    A_tilde = A_obs - E_A_hat;
    y_tilde = L_obs - E_A_hat * x_jaz;
    
    % 更新参数
    x_hat_new = (A_tilde' * Q_y_tilde_inv * A_tilde) \ (A_tilde' * Q_y_tilde_inv * y_tilde);
    
    % 收敛检查
    delta = norm(x_hat_new - x_jaz);
    
    x_jaz = x_hat_new;
    
    if delta < epsilon
        break;
    end
end

%% 计算单位权方差和参数协方差矩阵
cj=iter;%jaz法迭代次数
v_jaz = L_obs - A_obs * x_jaz;

% 重新计算Q_y_tilde用于最终参数
x_kron_T = kron(x_jaz', eye(n));
x_kron = kron(x_jaz, eye(n));
Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
Q_y_tilde_inv = inv(Q_y_tilde);

sigma_0_sq = (v_jaz' * Q_y_tilde_inv * v_jaz) / (n - m);
sigma_jaz_sqrt = sqrt(sigma_0_sq);
t_jaz = toc;
