function IDBO_UWSN_3D_scenario2
    % 场景二的仿真参数
    dim = 3; % 三维空间
    pop_size = 30; % 种群数量
    max_iter = 1000; % 最大迭代次数
    n_sensors = 100; % 节点数量
    sensing_radius = 10; % 感知半径
    communication_radius = 10; % 通信半径
    search_space = [0 60; 0 60; 0 60]; % 仿真区域 (x, y, z 的上下界)
    R = 1 - (1:max_iter) / max_iter; % 自适应权重因子
    coverage_history = zeros(max_iter, 1); % 记录每次迭代的覆盖率
    move_distances = zeros(max_iter, 1); % 记录节点的移动距离
    run_times = zeros(max_iter, 1); % 记录每次迭代的运行时间
    
    % 随机初始化种群
    pop = chebyshev_map(pop_size, n_sensors * dim, search_space);
    fitness = zeros(pop_size, 1);
    
    % 评估初始种群并绘制初始节点分布图
    for i = 1:pop_size
        sensor_positions = reshape(pop(i, :), [n_sensors, dim]); % 将平铺的一维数组转换为n_sensors x 3矩阵
        fitness(i) = coverage_metric(sensor_positions, n_sensors, sensing_radius, search_space);
    end
    
    % 绘制初始节点分布图，并表示每个节点的感知范围
    sensor_positions = reshape(pop(1, :), [n_sensors, dim]);
    figure;
    plot_nodes_with_range(sensor_positions, sensing_radius);
    title('初始节点分布');
    xlabel('X (米)');
    ylabel('Y (米)');
    zlabel('Z (米)');
    grid on;
    
    % 初始化全局最优解
    [best_fitness, best_idx] = min(fitness);
    gbest = pop(best_idx, :);
    
    % IDBO主循环
    for iter = 1:max_iter
        tic; % 开始计时
        prev_gbest = gbest; % 保存全局最优解的前一位置
        
        for i = 1:pop_size
            % 更新蜣螂位置，使用Levy飞行策略
            pop(i, :) = levy_flight(pop(i, :), gbest, R(iter), search_space, n_sensors, dim);
            
            % 重新评估适应度
            sensor_positions = reshape(pop(i, :), [n_sensors, dim]); % 重新整理为n_sensors x 3矩阵
            fitness(i) = coverage_metric(sensor_positions, n_sensors, sensing_radius, search_space);
        end
        
        % 更新全局最优解
        [current_best_fitness, best_idx] = min(fitness);
        if current_best_fitness < best_fitness
            best_fitness = current_best_fitness;
            gbest = pop(best_idx, :);
        end
        
        % 记录每次迭代的覆盖率
        coverage_history(iter) = -best_fitness;
        
        % 计算并记录节点的移动距离
        move_distances(iter) = sum(vecnorm(gbest - prev_gbest));
        
        % 记录每次迭代的运行时间
        run_times(iter) = toc;
        
        % 显示当前迭代信息
        disp(['迭代次数: ' num2str(iter) ', 最优覆盖率: ' sprintf('%.2f', -best_fitness*100) '%']);
    end

    % 绘制IDBO算法优化后的节点分布图，并表示每个节点的感知范围
    sensor_positions = reshape(gbest, [n_sensors, dim]); % 获取最终全局最优的节点位置
    figure;
    plot_nodes_with_range(sensor_positions, sensing_radius);
    title('IDBO算法优化后的节点分布');
    xlabel('X (米)');
    ylabel('Y (米)');
    zlabel('Z (米)');
    grid on;
    
    % 绘制覆盖率曲线图
    figure;
    plot(1:max_iter, coverage_history, 'LineWidth', 2);
    title('覆盖率曲线');
    xlabel('迭代次数');
    ylabel('覆盖率');
    grid on;
    
    % 绘制算法的运行时间图
    figure;
    plot(1:max_iter, run_times, 'LineWidth', 2);
    title('算法运行时间');
    xlabel('迭代次数');
    ylabel('时间 (秒)');
    grid on;
    
    % 绘制算法节点的移动距离图
    figure;
    plot(1:max_iter, move_distances, 'LineWidth', 2);
    title('节点移动距离');
    xlabel('迭代次数');
    ylabel('移动距离');
    grid on;
end

% 绘制节点及其感知范围的函数
function plot_nodes_with_range(sensor_positions, sensing_radius)

%输入：sensor_positions 100x3的矩阵，表示100个节点的坐标
%输入：sensing_radius 半径
hold on;
[num_sensors, ~] = size(sensor_positions);
for i = 1:num_sensors
    sensor_pos = sensor_positions(i, :);
    [X, Y, Z] = sphere(30); % 30 是球体的分辨率，可以调节哦，30够了
    X = X * sensing_radius + sensor_pos(1);
    Y = Y * sensing_radius + sensor_pos(2);
    Z = Z * sensing_radius + sensor_pos(3);
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', '[0.8, 0.8, 0.8]');
end

% 绘制传感器节点（用不同的颜色和大小表示）
scatter3(sensor_positions(:,1), sensor_positions(:,2), sensor_positions(:,3), 50, 'r', 'filled');

% 设置图形属性
xlabel('X');
ylabel('Y');
zlabel('Z');
title('三维传感器网络节点分布及感知范围');
grid on;
axis equal;
view(45, 30); % 调整视角
colormap(jet); % 设置颜色图

% 添加颜色条
colorbar;
caxis([0, num_sensors]);

hold off;

end

% 这里就用的随机，混不混沌影响不大，论文里可以水
function pop = chebyshev_map(pop_size, n_dims, bounds)
    ub = bounds(1,2);lb = bounds(1,1);
    pop = [];
    for i = 1:pop_size
        temp = initialization(n_dims/3,3,ub,lb);
        temp = reshape(temp, 1, []);
        pop(i,:) = temp;
    end
end

% Levy飞行位置更新函数
function new_pos = levy_flight(pos, gbest, inertia, bounds, n_sensors, dim)
    beta = 1.5; % Levy飞行参数
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(size(pos)) * sigma;
    v = randn(size(pos));
    step = u ./ abs(v).^(1 / beta);
    new_pos = pos + inertia .* step .* (pos - gbest);
    
    % 确保新位置在搜索空间内
    for i = 1:dim
        idx = (i-1)*n_sensors + (1:n_sensors); % 获取对应维度的索引
        new_pos(idx) = max(new_pos(idx), bounds(i, 1));
        new_pos(idx) = min(new_pos(idx), bounds(i, 2));
    end
end

% 目标函数：布尔感知模型逼近覆盖率
function coverage = coverage_metric(sensor_positions, n_sensors, sensing_radius, bounds)
    %输入：sensor_positions 100x3的矩阵，表示100个点的坐标
    %输入：n_sensors 没用到
    %输入：sensing_radius 半径
    %输入：bounds 边界条件（0，60）
    
    select_num = 1000;%选取1000个点作为参考点，越大精度越高，如果感觉速度慢可以取500
    matrix = init_jiadianji(select_num,3,bounds(1,2),bounds(1,1));%用佳点集均匀取点
    counters = 0;
    for i = 1:select_num
        label = any((sqrt((sensor_positions(:,1) - matrix(i,1)).^2 + (sensor_positions(:,2) - matrix(i,2)).^2 + ...
            (sensor_positions(:,3) - matrix(i,3)).^2)<sensing_radius));
        if label == 1
            counters = counters+1;
        end
    end
    coverage = round(-counters/select_num,8);
end
