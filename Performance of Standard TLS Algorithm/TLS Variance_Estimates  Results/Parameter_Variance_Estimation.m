dd = importdata('PTLS-data.txt'); %读取数据
n = size(dd,1);
m = 2;
k = m + 1;
A_obs = [dd(:,1),ones(n,1)];
L_obs = dd(:,3);

%计算权阵
P = [dd(:,4),dd(:,2),10^15*ones(n,1)]';
py = P(1,:);
Py = diag(py);
pA = P(2,:);
A1 = dd(:,1);
PP = diag(P(:));

X = TLS_XG_newton3(A_obs, L_obs, PP);
v = L_obs - A_obs * X;
% 计算误差传播
Q_e = inv(PP);
Q_total = [P(1,:),P(2,:),P(3,:)];
Q_total = inv(diag(Q_total));
[H, e_A, B, e] = Hessian(A_obs, L_obs, PP, X);
e_hat_reshaped = reshape(e, k, n)';
e_L = e_hat_reshaped(:, 1);


% 计算J矩阵[J0,J1,J2]2*30
J_total = zeros(m, k*n);

% 1.计算 P_v = (B * Q_e * B')^{-1}
Pv_inv = B * Q_e * B';
if rcond(Pv_inv) < eps
    P_v = pinv(Pv_inv);
else
    P_v = inv(Pv_inv);
end

% 2.计算 dx_dL = -H \ (A + e_A)' * P_v;
Gamma = (A_obs + e_A)' * P_v;
if rcond(H) < eps
    dx_dL = -pinv(H) * Gamma;  % m×n
else
    dx_dL = -H \ Gamma;        % 更高效，等价于 -inv(H) * Gamma
end
% 3.计算J0

J0 =dx_dL;
% 4.计算Ji
for param_idx = 1:m  % 对每个参数列（A1, A2, ...）
    %         Sigma_a_i = sit0_1 * inv(Pa{param_idx});
    for obs_idx = 1:n  % 对每个观测值
        % 计算R1和R2矩阵（针对特定参数和观测值）
        R1 = zeros(n, m);
        R1(obs_idx, param_idx) = 1;  % 第obs_idx行，第param_idx列设为1

        R2 = zeros(n, 1);
        R2(obs_idx) = -X(param_idx);  % -X(param_idx)

        % 计算∂e/∂a_{param_idx,obs_idx}
        de_da = Q_e * B' * P_v * R2;
        dET_da = zeros(m, n);
        for j = 1:n
            dET_da(:, j) = de_da((m+1)*(j-1)+2:(m+1)*(j-1)+(m+1));
        end
        dF_da_single = (R1' + dET_da) * P_v * v + (A_obs + e_A)' * P_v * R2;
        dF_da_i(:, obs_idx) = dF_da_single;
        if rcond(H) < eps
            dx_da_i = -pinv(H) * dF_da_i;  % m×n
        else
            dx_da_i = -H \ dF_da_i;        % 更高效，等价于 -inv(H) * de_da_i
        end
        da=zeros(n,m);
        da(obs_idx,param_idx)=1;
        term1(:,obs_idx) = da * X ;
    end
    Ja{param_idx} = dx_da_i;

end
J_total = [J0, Ja{1}, Ja{2}];

% 单位权方差
sit0_1 = (v' * P_v * v) / (n - m);
Dx_newton = sit0_1 *  J_total * Q_total * J_total';

%% ============== Wang法 ==============

[X, sigma2, Qxx, ci] = Wang_sigma(n, A_obs, L_obs, py, pA,0);
mm_Wang=sqrt(sigma2);
Dx_Wang = sigma2*Qxx;



%% ============== Xu法 ==============
[beta_hat, a_hat, iter, conv_flag, S0, A0, PA] = Xu_m(A_obs, L_obs, P);

% 修正：使用观测值（而非真实值）计算单位权方差
% 注意：实际应用中无法获得真实值A1_true，只能使用观测值a
a_hat0 = A_obs(:,1);
r_a = a_hat0(:) - a_hat;            % 随机元素的残差（使用观测值a）
r_y = L_obs - A0 * beta_hat;       % 观测向量的残差

% 计算单位权方差
sigma2_hat = (r_a' * PA * r_a + r_y' * Py * r_y) / (n - m);
% 计算法方程矩阵
N11 = A0' * Py * A0;
N12 = A0' * Py * S0;
N21 = S0' * Py * A0;
N22 = S0' * Py * S0 + PA;

% 计算参数β的方差-协方差矩阵
Dx_xu = inv(N11 - N12 * inv(N22) * N21) * sigma2_hat;

%% ============== jazaeri法 ==============
% 初始化权矩阵
Q_y = diag(1./py);
Q_A = zeros(m*n, m*n);
% Q_A结构：vec(A)按列堆叠，所以前m个元素是x列，后m个是常数列
Q_A(1:n, 1:n) = diag(1./pA);           % x列有误差
Q_A(n+1:m*n, n+1:m*n) = zeros(n); % 常数列无误差

% 从加权最小二乘开始
P0 = diag(py);
x_jaz = inv(A_obs' * P0 * A_obs) * (A_obs' * P0 * L_obs);
epsilon = 1e-10;

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

% 重新计算Q_y_tilde用于最终参数
v_jaz = L_obs - A_obs * x_jaz;

% 重新计算Q_y_tilde用于最终参数
x_kron_T = kron(x_jaz', eye(n));
x_kron = kron(x_jaz, eye(n));
Q_y_tilde = Q_y + x_kron_T * Q_A * x_kron;
Q_y_tilde_inv = inv(Q_y_tilde);

sigma_0_sq = (v_jaz' * Q_y_tilde_inv * v_jaz) / (n - m);
sigma_0 = sqrt(sigma_0_sq);

% 参数协方差矩阵计算
vec_E_A = -Q_A * x_kron * Q_y_tilde_inv * v_jaz;
E_A_hat = [vec_E_A(1:n), vec_E_A(n+1:end)];
A_tilde = A_obs - E_A_hat;
Q_x = inv(A_tilde' * Q_y_tilde_inv * A_tilde);

Sigma_x = sigma_0_sq * Q_x;