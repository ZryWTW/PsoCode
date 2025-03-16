%% 清空环境变量
clc
clear

%% 网络参数
L = 60;                  % 区域边长
n = 70;                  % 节点个数
R = 10;                  % 通信半径
data = 1;                % 离散粒度

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
gbest = zeros(n, 3);    % 最优解

%% 随机生成群体位置、速度和对应的适应度值
empty_pop.Position = [];
empty_pop.V = [];
pop = repmat(empty_pop, sizepop, 1);

for i = 1:sizepop
    pop(i).Position = rand(n, 3) .* L;  % 初始种群位置 (三维)
    pop(i).V = randi([-2, 2], n, 3);     % 初始化速度 (范围为-2到2)
    fitness(i) = fun(pop(i).Position(:, 1), pop(i).Position(:, 2), pop(i).Position(:, 3), L, R, data);  % 粒子群的适应度
end

[bestfitness, bestindex] = max(fitness);
gbest = pop(bestindex).Position;    % 群体最优极值
pbest = pop;                        % 个体最优极值
fitnessgbest = bestfitness;         % 种群最优适应度值
fitnesspbest = fitness;             % 个体最优适应度值

%% 初始结果显示
disp('初始位置：' );
disp([num2str(gbest)]);
disp(['初始覆盖率：', num2str(fitnessgbest)]);

% 初始覆盖图
figure
for i = 1:n
    axis([0 L 0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    z = gbest(:, 3);
    theta = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    % 画出每个节点的覆盖区域
    [X, Y, Z] = sphere;  % 生成单位球体
    surf(R*X + x(i), R*Y + y(i), R*Z + z(i), 'FaceAlpha', 0.3); % 绘制透明覆盖区域
end
p1 = plot3(gbest(:, 1), gbest(:, 2), gbest(:, 3), 'r*');
legend([p1], {'WSNs节点'});
title '初始节点位置及覆盖区域';
xlabel('x');
ylabel('y');
zlabel('z');
grid on;

%% 迭代寻优
for i = 1:maxgen
    W = Wmax - ((Wmax - Wmin) / maxgen) * i; % 更新惯性权重
    for j = 1:sizepop
        %% 速度更新
        pop(j).V = W * pop(j).V + c1 * rand(n, 3) .* (pbest(j).Position - pop(j).Position) + c2 * rand(n, 3) .* (gbest - pop(j).Position);
        
        % 边界处理
        pop(j).V = max(pop(j).V, Vmin);
        pop(j).V = min(pop(j).V, Vmax);
        
        %% 位置更新
        pop(j).Position = pop(j).Position + pop(j).V;
        
        % 边界处理
        pop(j).Position = max(pop(j).Position, popmin);
        pop(j).Position = min(pop(j).Position, popmax);
        
        %% 适应度值更新
        fitness(j) = fun(pop(j).Position(:, 1), pop(j).Position(:, 2), pop(j).Position(:, 3), L, R, data);    
    end
    
    %% 个体和群体极值更新
    for j = 1:sizepop
        % 个体极值更新
        if fitness(j) > fitnesspbest(j)
            pbest(j).Position = pop(j).Position;
            fitnesspbest(j) = fitness(j);
        end
        
        % 群体极值更新
        if fitness(j) > fitnessgbest
            gbest = pop(j).Position;
            fitnessgbest = fitness(j);
        end
    end
    
    %% 每一代群体最优值存入zz数组
    zz(i) = fitnessgbest;
end

%% 结果显示
disp('最优位置：' );
disp([num2str(gbest)]);
disp(['最优覆盖率：', num2str(zz(end))]);

%% 绘图
figure;
plot(zz, 'r', 'lineWidth', 2);          % 画出迭代图
title('算法过程', 'fontsize', 12);
xlabel('迭代次数', 'fontsize', 12);
ylabel('粒子覆盖率', 'fontsize', 12);

figure
for i = 1:n
    axis([0 L 0 L 0 L]);        % 限制坐标范围
    x = gbest(:, 1);
    y = gbest(:, 2);
    z = gbest(:, 3);
    theta = 0:pi/100:2*pi;   % 角度[0, 2*pi]
    hold on;
    % 画出每个节点的覆盖区域
    [X, Y, Z] = sphere;  % 生成单位球体
    surf(R*X + x(i), R*Y + y(i), R*Z + z(i), 'FaceAlpha', 0.3); % 绘制透明覆盖区域
end
p1 = plot3(gbest(:, 1), gbest(:, 2), gbest(:, 3), 'r*');
legend([p1], {'WSNs节点'});
title 'PSO优化后节点位置及覆盖区域';
xlabel('x');
ylabel('y');
zlabel('z');
grid on;

%% 适应度函数
function z = fun(x, y, z, L, R, data)
%% 适应度函数：WSNs的覆盖率
% input：
% x        圆心横坐标
% y        圆心纵坐标
% z        圆心纵深坐标
% L        区域边长
% R        通信半径
% data     离散粒度
% output:
% z        覆盖率

N = length(x);                      % 节点总个数
[m, n, p] = meshgrid(0:data:L);    % 离散化区域内的点
[row, col, dep] = size(m);
M = zeros(row, col, dep);          % 初始化覆盖状态矩阵

for i = 1:N
    D = sqrt((m-x(i)).^2 + (n-y(i)).^2 + (p-z(i)).^2);   % 计算坐标点到圆心的距离
    M(D <= R) = 1;                  % 改变覆盖状态
end
scale = sum(M(:)) / (row * col * dep); % 计算覆盖比例
z = scale;
end 