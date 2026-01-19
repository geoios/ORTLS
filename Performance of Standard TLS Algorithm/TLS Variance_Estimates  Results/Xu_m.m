function [X, a_hat, iter, conv_flag, S0, A0, PA] = Xu_m(A, L, P)
    % 通用部分EIV模型TLS算法
    % 输入:
    %   A: 设计矩阵 (N×M), 第M列为全1 (常数项)
    %   L: 观测向量 (N×1)
    %   P: 权矩阵 (M×N), 第1行: y的权重, 第2到M行: 各随机自变量的权重
    % 输出:
    %   beta_hat: 参数估计 (M×1)
    %   a_hat: 随机变量真值估计 (N*(M-1)×1)
    %   iter: 迭代次数


    % 初始化收敛标志为0（未收敛）
    conv_flag = 0;
    % 获取问题维度
    [N, M] = size(A);        % N=观测点数, M=参数个数
    rand_count = M - 1;       % 随机自变量个数
    t = N * rand_count;       % 随机元素总个数
    
    % 提取权和构建对角权矩阵
    py = P(1, :)';           % y的权值 (N×1)
    P_rand = P(2:M, :);      % 随机自变量的权值 ((M-1)×N)
    
    Py = diag(py);           % y的权矩阵 (N×N)
    
    % 构建随机自变量的权重向量
    w_rand = zeros(t, 1);
    for i = 1:rand_count
        start_idx = (i-1)*N + 1;
        end_idx = i*N;
        w_rand(start_idx:end_idx) = P_rand(i, :)';
    end
    PA = diag(w_rand);       % 随机变量的权矩阵 (t×t)
    
    % 构建部分EIV模型参数
    a_obs = A(:, 1:rand_count);  % 随机自变量观测值
    a_obs = a_obs(:);            % 转换为列向量 (t×1)
    
    % 构建h向量 (N*M×1)
    h = zeros(N*M, 1);
    h((M-1)*N+1:end) = ones(N, 1);  % 常数项位置设为1
    
    % 构建B矩阵 (N*M×t)
    B = zeros(N*M, t);
    for i = 1:rand_count
        start_row = (i-1)*N + 1;
        end_row = i*N;
        start_col = (i-1)*N + 1;
        end_col = i*N;
        B(start_row:end_row, start_col:end_col) = eye(N);
    end
    
    % 将h和B分块存储
    h_blocks = mat2cell(h, N*ones(M, 1), 1);    % M个N×1向量
    B_blocks = mat2cell(B, N*ones(M, 1), t);    % M个N×t矩阵
    
    % 初始化
    beta_hat0 = (A'*Py*A)\(A'*Py*L);  % LS初始值
    a_hat = a_obs;                    % 初始随机变量估计
    tol = 1e-10;
    max_iter = 1000;
    iter = 0;
    converged = false;  % 新增收敛状态变量
    X_first = [];
    tic;
    for iter = 1:max_iter
        % ========== 计算beta_hat (公式31) ==========
        N_h = zeros(M, M);
        N_B = zeros(M, M);
        N_Bh = zeros(M, M);
        N_hB = zeros(M, M);
        u_h = zeros(M, 1);
        u_B = zeros(M, 1);
        
        for i = 1:M
            for j = 1:M
                N_h(i, j) = h_blocks{i}' * Py * h_blocks{j};
                N_B(i, j) = a_hat' * B_blocks{i}' * Py * B_blocks{j} * a_hat;
                N_Bh(i, j) = a_hat' * B_blocks{i}' * Py * h_blocks{j};
                N_hB(i, j) = h_blocks{i}' * Py * B_blocks{j} * a_hat;
            end
            u_h(i) = h_blocks{i}' * Py * L;
            u_B(i) = a_hat' * B_blocks{i}' * Py * L;
        end
        
        % 计算当前beta_hat
        N_total = N_h + N_B + N_Bh + N_hB;
        u_total = u_h + u_B;
        X = N_total \ u_total;
        % ========== 保存第一次迭代的X值 ==========
        if iter == 1
            X_first = X;
        end
        
        % ========== 更新a_hat (公式28) ==========
        % 计算S_beta = Σβ_i*B_i
        S_beta = zeros(N, t);
        sum_beta_h = zeros(N, 1);
        for i = 1:M
            S_beta = S_beta + X(i) * B_blocks{i};
            sum_beta_h = sum_beta_h + X(i) * h_blocks{i};
        end
        % 应用公式(28)
        a_hat_new = (PA + S_beta' * Py * S_beta) \ ...
                   (PA * a_obs - S_beta' * Py * sum_beta_h + S_beta' * Py * L);
        
        % ========== 检查收敛 ==========
        delta_beta = norm(X - beta_hat0);
        delta_a = norm(a_hat_new - a_hat);
        
        if max(delta_beta, delta_a) < tol
            a_hat = a_hat_new;
            converged = true;  % 标记为收敛
            break;
        end
        
        % 更新迭代变量
        beta_hat0 = X;
        a_hat = a_hat_new;
    end
    time = toc;
    % 设置收敛标志
    if converged
        conv_flag = 1;  % 收敛
    else
        conv_flag = 0;  % 未收敛
        warning('未在最大迭代次数内收敛');
    end
        % ========== 计算输出S0和A0 (用于公式57) ==========
    % 1. 计算S0 (公式29): S_β = Σβ_i*B_i
    S0 = zeros(N, t);
    for i = 1:M
        S0 = S0 + X_first(i) * B_blocks{i};
    end
   

    % 2. 计算A0: 设计矩阵真值估计
    % 将a_hat重塑为N×(M-1)矩阵
    A_rand_hat = reshape(a_hat, N, rand_count);
    % 构建完整设计矩阵 (添加常数项列)
    A0 = [A_rand_hat, ones(N, 1)];  % 最后一项为常数项
end
