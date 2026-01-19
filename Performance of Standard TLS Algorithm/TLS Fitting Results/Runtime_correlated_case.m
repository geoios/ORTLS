% clear all; clc; close all;

% 实验参数设置
n_values = [10, 100:100:1000]; % 观测值数量数组
m = 2;                          % 参数个数
k = m + 1;                      % 每个观测点的变量数(L + A1 + A2)
num_n = length(n_values);       % 观测维度数量

% 设置噪声方差
sigma_x = 0.001;                % x的噪声方差
sigma_L = 0.001;                % L的噪声方差

% 初始化存储运行时间的数组
time_newton_full = zeros(num_n, 1);      % 牛顿法（全相关）
time_newton_indep = zeros(num_n, 1);     % 牛顿法（独立）

% 循环处理不同观测维度
for n_idx = 1:num_n
    n = n_values(n_idx);
    fprintf('=== 实验: 观测维度 n = %d ===\n', n);
    
    % 生成真实数据
    A1_true = 10 * rand(n, 1);      % [0,10]区间均匀分布
    A2_true = ones(n, 1);          % 常数项
    L_true = 2 * A1_true - A2_true;
    
    % 生成带噪声的观测数据
    A1_obs = A1_true + sigma_x * randn(n, 1);
    L_obs = L_true + sigma_L * randn(n, 1);
    A2_obs = A2_true;  % A2没有误差
    A_obs = [A1_obs, A2_obs];
    
    % 设置权阵
    py = (1/sigma_L) * ones(n, 1);  % L的权
    px1 = (1/sigma_x) * ones(n, 1); % A1的权
    px2 = 10^12 * ones(n, 1);       % A2的权（近似无噪声）
    pA = [px1, px2];
    P = [py, px1, px2]';            % 3 x n
    PP = diag(P(:));
    
    %% ============== 牛顿法（全相关）计时 ==============
    fprintf('  计算牛顿法（全相关）...\n');
    tic;
    X_newton_full = TLS_XG_newton3(A_obs, L_obs, PP);
    time_newton_full(n_idx) = toc;
    
    %% ============== 牛顿法（独立）计时 ==============
    fprintf('  计算牛顿法（独立）...\n');
    tic;
    X_newton_indep = TLS_newton_2(A_obs, L_obs, P);
    time_newton_indep(n_idx) = toc;
    
    fprintf('  观测维度 n = %d 完成\n\n', n);
end
time_newton_indep=time_newton_indep*1000;
%% ============== 绘制运行时间对比图 ==============
fprintf('绘制运行时间对比图...\n');

% 创建图形窗口
figure('Position', [100, 100, 1200, 800]);

% 设置图形参数
markerSize = 10;
lineWidth = 2;
fontSize = 20;
tickFontSize = 18;  % 坐标刻度字号设置为8

% 创建双纵轴
yyaxis left;
plot(n_values, time_newton_full, 'b-o', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'DisplayName', 'Newton (Full Correlation)');
ylabel('Runtime/s (correlated case)', 'FontSize', fontSize, 'FontWeight', 'bold');
ylim_left = [min(time_newton_full)*0.9, max(time_newton_full)*1.1];
ylim(ylim_left);
grid on;
hold on;

% 右侧纵轴
yyaxis right;
plot(n_values, time_newton_indep, 'r-s', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'DisplayName', 'Newton (Independent)');
ylabel('Runtime/ms (uncorrelated case)', 'FontSize', fontSize, 'FontWeight', 'bold');
ylim_right = [0, max(time_newton_indep) * 1.1];
ylim(ylim_right);

% 设置图形属性
xlabel('Observation Dimension (n)', 'FontSize', fontSize, 'FontWeight', 'bold');

% 设置坐标刻度字号为8
ax = gca;
ax.FontSize = tickFontSize;

% 添加图例（左上角）
legend('Location', 'northwest', 'FontSize', fontSize-2, 'NumColumns', 2);

% 设置坐标轴颜色以匹配线条颜色
ax.YAxis(1).Color = 'b';  % 左侧纵轴为蓝色
ax.YAxis(2).Color = 'r';  % 右侧纵轴为红色

% 添加网格
grid on;
box on;

% ========== 保存为600 DPI PNG图像 ==========
filename = 'runtime_comparison_two_newton_methods.png';
resolution = 600;
print('-dpng', ['-r' num2str(resolution)], filename);
disp(['图像已保存为: ' filename ' (600 DPI PNG格式)']);

%% ============== 显示运行时间统计信息 ==============
fprintf('\n运行时间统计信息:\n');
fprintf('观测维度: %s\n', mat2str(n_values));
fprintf('牛顿法（全相关）运行时间: %s 秒\n', mat2str(time_newton_full'));
fprintf('牛顿法（独立）运行时间: %s 秒\n', mat2str(time_newton_indep'));

fprintf('\n所有实验完成！\n');


%% ============== 显示运行时间统计信息 ==============
fprintf('\n运行时间统计信息:\n');
fprintf('观测维度: %s\n', mat2str(n_values));
fprintf('牛顿法（全相关）运行时间: %s 秒\n', mat2str(time_newton_full'));
fprintf('牛顿法（独立）运行时间: %s 秒\n', mat2str(time_newton_indep'));

fprintf('\n所有实验完成！\n');