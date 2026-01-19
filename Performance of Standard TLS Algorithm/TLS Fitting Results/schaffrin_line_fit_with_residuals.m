%% Schaffrin (2015)算法实现 - 包含残差计算
function [X, sigma2, ey, eA, iterations, converged,cc] = schaffrin_line_fit_with_residuals(A, L, Q_y, Q_A, Q_yA, max_iter, tol)
    % 修正后的Schaffrin (2015)算法D实现
    % 严格按照公式(9c)实现，并添加残差计算
    
    n = size(A, 1);
    m = size(A, 2);  % 对于直线拟合 m=2
    
    % 定义Q_Ay
    Q_Ay = Q_yA';
    
    % 初始估计：使用加权最小二乘
    if rank(Q_y) == size(Q_y, 1)
        X = (A' / Q_y * A) \ (A' / Q_y * L);
    else
        X = A \ L;
    end
    
    % 预计算
    I_n = speye(n);
    I_m = speye(m);
    
    % 迭代求解
    converged = false;
    cc=0;
    for iterations = 1:max_iter
        cc=cc+1;
        X_old = X;
        
        % 1. 构建B矩阵: B = [I_n, -(\xi^T ⊗ I_n)]
        kron_term = kron(X', I_n);  % n × (n*m)
        B = [I_n, -kron_term];  % n × (n + n*m)
        
        % 2. 构建完整的协因数矩阵Q
        Q = [Q_y, Q_yA; 
             Q_Ay, Q_A];
        
        % 3. 计算Q1 = B * Q * B'
        Q1 = B * Q * B';
        
        % 4. 检查并处理病态Q1
        if cond(Q1) > 1e12
            Q1 = Q1 + 1e-8 * speye(size(Q1));
        end
        
        % 5. 计算拉格朗日乘子λ = Q1^{-1} * (L - A*X)
        residual = L - A * X;
        
        % 使用Cholesky分解高效求解
        try
            R = chol(Q1);
            lambda = R \ (R' \ residual);
            Q1_inv_A = R \ (R' \ A);
            Q1_inv_L = R \ (R' \ L);
        catch
            % 如果Cholesky失败，使用直接求解
            lambda = Q1 \ residual;
            Q1_inv_A = Q1 \ A;
            Q1_inv_L = Q1 \ L;
        end
        
        % 6. 构建 (I_m ⊗ λ)
        I_lambda = kron(I_m, lambda);  % n*m × m
        
        % 7. 计算M和N矩阵（根据公式9c）
        % M = A^T * Q1^{-1} * A - (I_lambda)^T * Q_A * I_lambda
        M_part1 = A' * Q1_inv_A;
        M_part2 = I_lambda' * Q_A * I_lambda;
        M = M_part1 - M_part2;
        
        % N = A^T * Q1^{-1} * L - (I_lambda)^T * Q_Ay * λ
        N_part1 = A' * Q1_inv_L;
        N_part2 = I_lambda' * Q_Ay * lambda;
        N = N_part1 - N_part2;
        
        % 8. 求解新的参数估计
        if cond(M) > 1e12
            M = M + 1e-8 * eye(m);
        end
        
        X_new = M \ N;
        
        % 9. 检查收敛
        delta_X = norm(X_new - X_old);
        if delta_X < tol
            X = X_new;
            converged = true;
            break;
        end
        
        X = X_new;
    end
    
    if iterations == max_iter && ~converged
        warning('达到最大迭代次数 %d 仍未收敛', max_iter);
    end
    
    % ==============================================
    % 计算最终残差ey和eA（根据公式8c）
    % [ey; eA] = Q * B' * Q1^{-1} * (y - A*ξ)
    % ==============================================
    
    % 使用最终参数重新计算相关矩阵
    kron_term_final = kron(X', I_n);
    B_final = [I_n, -kron_term_final];  % n × (n + n*m)
    
    % 构建完整协因数矩阵Q
    Q_full = [Q_y, Q_yA; 
              Q_Ay, Q_A];  % (n + n*m) × (n + n*m)
    
    % 计算Q1_final = B * Q * B'
    Q1_final = B_final * Q_full * B_final';
    
    % 计算最终残差向量
    final_residual = L - A * X;
    
    % 计算λ_final = Q1_final^{-1} * (L - A*X)
    try
        R_final = chol(Q1_final);
        lambda_final = R_final \ (R_final' \ final_residual);
    catch
        lambda_final = Q1_final \ final_residual;
    end
    
    % 计算残差向量 [ey; eA] = Q * B' * λ_final
    res_vec = Q_full * B_final' * lambda_final;
    
    % 分解残差向量
    ey = res_vec(1:n);  % y的残差，n×1向量
    eA_vec = res_vec(n+1:end);  % A的残差向量化形式，n*m×1向量
    
    % 将eA向量重构为n×m矩阵
    eA = reshape(eA_vec, n, m);
    
    % 计算方差分量估计（公式8d）
    % σ² = (y - Aξ)^T * Q1^{-1} * (y - Aξ) / (n - m)
    sigma2 = (final_residual' / Q1_final * final_residual) / (n - m);
   
end