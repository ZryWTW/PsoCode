%% 清空环境变量
clc;
clear;

%% 网络参数
L = 60;                  % 区域边长
n = 70;                  % 节点个数
R = 10;                  % 通信半径
data = 1;                % 离散粒度（若执行速度太慢，可尝试调大，例如 data = 2 或 3）

%% 粒子群参数
maxgen = 500;            % 迭代次数
sizepop = 30;            % 粒子规模
Wmax = 0.9;
Wmin = 0.4;

%% 参数初始化
c1 = 2;                  % 自我认知参数
c2 = 2;                  % 社会认知参数
Vmax = 2;                % 最大速度
Vmin = -2;               % 最小速度
popmax = L;              % 位置最大值
popmin = 0;              % 位置最小值

%% 预先计算网格点（只计算一次）
[m_grid, n_grid, p_grid] = meshgrid(0:data:L);
grid_points = [m_grid(:), n_grid(:), p_grid(:)];  % 将网格点转换为 M x 3 的矩阵

%% 初始化粒子群及适应度
empty_pop.Position = [];
empty_pop.V = [];
pop = repmat(empty_pop, sizepop, 1);
fitness = zeros(sizepop, 1);  % 保存每个粒子的适应度

for i = 1:sizepop
    pop(i).Position = rand(n, 3) * L;         % 随机初始位置 (n x 3)
    pop(i).V = randi([-2, 2], n, 3);            % 随机初始速度 (n x 3)
    % 使用向量化适应度函数计算覆盖率
    fitness(i) = fun_vec(pop(i).Position(:,1), pop(i).Position(:,2), pop(i).Position(:,3), grid_points, R);
end

[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex).Position;    % 群体最优位置
gbest_initial = gbest;              % 保存初始群体最优位置，用于对比显示
pbest = pop;                        % 个体最优位置初始化
fitnessgbest = bestfitness;         % 种群最优适应度
fitnesspbest = fitness;             % 个体最优适应度

%% 显示初始结果
disp('初始位置：');
disp(gbest_initial);
disp(['初始覆盖率：', num2str(fitnessgbest)]);

%% 迭代寻优（PSO算法）
zz = zeros(maxgen,1);   % 保存每一代最优适应度
for gen = 1:maxgen
    % 更新惯性权重
    W = Wmax - ((Wmax - Wmin) / maxgen) * gen;
    
    for j = 1:sizepop
        %% 速度更新
        pop(j).V = W * pop(j).V + ...
                   c1 * rand(n, 3) .* (pbest(j).Position - pop(j).Position) + ...
                   c2 * rand(n, 3) .* (gbest - pop(j).Position);
        % 速度边界处理
        pop(j).V = max(pop(j).V, Vmin);
        pop(j).V = min(pop(j).V, Vmax);
        
        %% 位置更新
        pop(j).Position = pop(j).Position + pop(j).V;
        % 位置边界处理
        pop(j).Position = max(pop(j).Position, popmin);
        pop(j).Position = min(pop(j).Position, popmax);
        
        %% 适应度更新（调用向量化适应度函数）
        fitness(j) = fun_vec(pop(j).Position(:,1), pop(j).Position(:,2), pop(j).Position(:,3), grid_points, R);
    end
    
    %% 个体和群体最优更新
    for j = 1:sizepop
        % 更新个体最优
        if fitness(j) > fitnesspbest(j)
            pbest(j).Position = pop(j).Position;
            fitnesspbest(j) = fitness(j);
        end
        % 更新全局最优
        if fitness(j) > fitnessgbest
            gbest = pop(j).Position;
            fitnessgbest = fitness(j);
        end
    end
    
    zz(gen) = fitnessgbest;
end

%% 显示最终结果
disp('最优位置：');
disp(gbest);
disp(['最优覆盖率：', num2str(zz(end))]);

%% 绘图部分

% 图1：初始节点位置及覆盖区域
figure(1);
hold on;
for i = 1:n
    [X, Y, Z] = sphere;  % 生成单位球体
    surf(R*X + gbest_initial(i,1), R*Y + gbest_initial(i,2), R*Z + gbest_initial(i,3), ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
plot3(gbest_initial(:,1), gbest_initial(:,2), gbest_initial(:,3), 'r*', 'MarkerSize', 8);
axis([0 L 0 L 0 L]);
axis equal;
view(3);
rotate3d on;
title('初始节点位置及覆盖区域');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
hold off;

% 图2：算法迭代过程（覆盖率曲线）
figure(2);
plot(1:maxgen, zz, 'r', 'LineWidth', 2);
title('算法过程', 'FontSize', 12);
xlabel('迭代次数', 'FontSize', 12);
ylabel('粒子覆盖率', 'FontSize', 12);
grid on;

% 图3：PSO优化后节点位置及覆盖区域
figure(3);
hold on;
for i = 1:n
    [X, Y, Z] = sphere;  % 生成单位球体
    surf(R*X + gbest(i,1), R*Y + gbest(i,2), R*Z + gbest(i,3), ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');
end
plot3(gbest(:,1), gbest(:,2), gbest(:,3), 'r*', 'MarkerSize', 8);
axis([0 L 0 L 0 L]);
axis equal;
view(3);
rotate3d on;
title('PSO优化后节点位置及覆盖区域');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
hold off;

%% 向量化适应度函数：计算WSNs覆盖率
function scale = fun_vec(x, y, z, grid_points, R)
    % 输入：
    %   x, y, z      - 节点坐标的列向量 (n x 1)
    %   grid_points  - 网格中所有点的坐标矩阵 (M x 3)
    %   R            - 通信半径
    % 输出：
    %   scale        - 覆盖率（覆盖的网格点比例）
    
    sensor_pos = [x, y, z];  % n x 3
    % 计算每个网格点到所有节点的欧氏距离，返回 M x n 矩阵
    distances = pdist2(grid_points, sensor_pos);
    % 判断每个网格点是否至少被一个节点覆盖
    covered = any(distances <= R, 2);
    scale = sum(covered) / size(grid_points, 1);
end
