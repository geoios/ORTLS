function  e_hat = E_hat(A, L, P, x)



    % 获取问题维度
    [n, ~] = size(A);

    
    % 只计算一次Q_e和I_n，不预提取块
    Q_e = pinv(P);
    I_n = eye(n);

    % TLS牛顿迭代法

   % 1. 计算残差向量
    v = L - A * x;
        
    % 2. 构建b向量和B矩阵
    b = [1, x'];
    B = kron(I_n, b);
        
    % 3. 计算 P_v = (B * Q_e * B')^{-1}
    Pv_inv = B * Q_e * B';
    if rcond(Pv_inv) < eps
        P_v = pinv(Pv_inv);
    else
        P_v = inv(Pv_inv);
    end
        
    % 预计算常用项
    vp = P_v * v;
    B_T_vp = B' * vp;
        
    % 4. 计算误差估计 ê
    e_hat = Q_e * B_T_vp;
end