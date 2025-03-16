function IDBO_UWSN_3D_scenario2
    % �������ķ������
    dim = 3; % ��ά�ռ�
    pop_size = 30; % ��Ⱥ����
    max_iter = 1000; % ����������
    n_sensors = 100; % �ڵ�����
    sensing_radius = 10; % ��֪�뾶
    communication_radius = 10; % ͨ�Ű뾶
    search_space = [0 60; 0 60; 0 60]; % �������� (x, y, z �����½�)
    R = 1 - (1:max_iter) / max_iter; % ����ӦȨ������
    coverage_history = zeros(max_iter, 1); % ��¼ÿ�ε����ĸ�����
    move_distances = zeros(max_iter, 1); % ��¼�ڵ���ƶ�����
    run_times = zeros(max_iter, 1); % ��¼ÿ�ε���������ʱ��
    
    % �����ʼ����Ⱥ
    pop = chebyshev_map(pop_size, n_sensors * dim, search_space);
    fitness = zeros(pop_size, 1);
    
    % ������ʼ��Ⱥ�����Ƴ�ʼ�ڵ�ֲ�ͼ
    for i = 1:pop_size
        sensor_positions = reshape(pop(i, :), [n_sensors, dim]); % ��ƽ�̵�һά����ת��Ϊn_sensors x 3����
        fitness(i) = coverage_metric(sensor_positions, n_sensors, sensing_radius, search_space);
    end
    
    % ���Ƴ�ʼ�ڵ�ֲ�ͼ������ʾÿ���ڵ�ĸ�֪��Χ
    sensor_positions = reshape(pop(1, :), [n_sensors, dim]);
    figure;
    plot_nodes_with_range(sensor_positions, sensing_radius);
    title('��ʼ�ڵ�ֲ�');
    xlabel('X (��)');
    ylabel('Y (��)');
    zlabel('Z (��)');
    grid on;
    
    % ��ʼ��ȫ�����Ž�
    [best_fitness, best_idx] = min(fitness);
    gbest = pop(best_idx, :);
    
    % IDBO��ѭ��
    for iter = 1:max_iter
        tic; % ��ʼ��ʱ
        prev_gbest = gbest; % ����ȫ�����Ž��ǰһλ��
        
        for i = 1:pop_size
            % ��������λ�ã�ʹ��Levy���в���
            pop(i, :) = levy_flight(pop(i, :), gbest, R(iter), search_space, n_sensors, dim);
            
            % ����������Ӧ��
            sensor_positions = reshape(pop(i, :), [n_sensors, dim]); % ��������Ϊn_sensors x 3����
            fitness(i) = coverage_metric(sensor_positions, n_sensors, sensing_radius, search_space);
        end
        
        % ����ȫ�����Ž�
        [current_best_fitness, best_idx] = min(fitness);
        if current_best_fitness < best_fitness
            best_fitness = current_best_fitness;
            gbest = pop(best_idx, :);
        end
        
        % ��¼ÿ�ε����ĸ�����
        coverage_history(iter) = -best_fitness;
        
        % ���㲢��¼�ڵ���ƶ�����
        move_distances(iter) = sum(vecnorm(gbest - prev_gbest));
        
        % ��¼ÿ�ε���������ʱ��
        run_times(iter) = toc;
        
        % ��ʾ��ǰ������Ϣ
        disp(['��������: ' num2str(iter) ', ���Ÿ�����: ' sprintf('%.2f', -best_fitness*100) '%']);
    end

    % ����IDBO�㷨�Ż���Ľڵ�ֲ�ͼ������ʾÿ���ڵ�ĸ�֪��Χ
    sensor_positions = reshape(gbest, [n_sensors, dim]); % ��ȡ����ȫ�����ŵĽڵ�λ��
    figure;
    plot_nodes_with_range(sensor_positions, sensing_radius);
    title('IDBO�㷨�Ż���Ľڵ�ֲ�');
    xlabel('X (��)');
    ylabel('Y (��)');
    zlabel('Z (��)');
    grid on;
    
    % ���Ƹ���������ͼ
    figure;
    plot(1:max_iter, coverage_history, 'LineWidth', 2);
    title('����������');
    xlabel('��������');
    ylabel('������');
    grid on;
    
    % �����㷨������ʱ��ͼ
    figure;
    plot(1:max_iter, run_times, 'LineWidth', 2);
    title('�㷨����ʱ��');
    xlabel('��������');
    ylabel('ʱ�� (��)');
    grid on;
    
    % �����㷨�ڵ���ƶ�����ͼ
    figure;
    plot(1:max_iter, move_distances, 'LineWidth', 2);
    title('�ڵ��ƶ�����');
    xlabel('��������');
    ylabel('�ƶ�����');
    grid on;
end

% ���ƽڵ㼰���֪��Χ�ĺ���
function plot_nodes_with_range(sensor_positions, sensing_radius)

%���룺sensor_positions 100x3�ľ��󣬱�ʾ100���ڵ������
%���룺sensing_radius �뾶
hold on;
[num_sensors, ~] = size(sensor_positions);
for i = 1:num_sensors
    sensor_pos = sensor_positions(i, :);
    [X, Y, Z] = sphere(30); % 30 ������ķֱ��ʣ����Ե���Ŷ��30����
    X = X * sensing_radius + sensor_pos(1);
    Y = Y * sensing_radius + sensor_pos(2);
    Z = Z * sensing_radius + sensor_pos(3);
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', '[0.8, 0.8, 0.8]');
end

% ���ƴ������ڵ㣨�ò�ͬ����ɫ�ʹ�С��ʾ��
scatter3(sensor_positions(:,1), sensor_positions(:,2), sensor_positions(:,3), 50, 'r', 'filled');

% ����ͼ������
xlabel('X');
ylabel('Y');
zlabel('Z');
title('��ά����������ڵ�ֲ�����֪��Χ');
grid on;
axis equal;
view(45, 30); % �����ӽ�
colormap(jet); % ������ɫͼ

% �����ɫ��
colorbar;
caxis([0, num_sensors]);

hold off;

end

% ������õ�������첻����Ӱ�첻�����������ˮ
function pop = chebyshev_map(pop_size, n_dims, bounds)
    ub = bounds(1,2);lb = bounds(1,1);
    pop = [];
    for i = 1:pop_size
        temp = initialization(n_dims/3,3,ub,lb);
        temp = reshape(temp, 1, []);
        pop(i,:) = temp;
    end
end

% Levy����λ�ø��º���
function new_pos = levy_flight(pos, gbest, inertia, bounds, n_sensors, dim)
    beta = 1.5; % Levy���в���
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)))^(1 / beta);
    u = randn(size(pos)) * sigma;
    v = randn(size(pos));
    step = u ./ abs(v).^(1 / beta);
    new_pos = pos + inertia .* step .* (pos - gbest);
    
    % ȷ����λ���������ռ���
    for i = 1:dim
        idx = (i-1)*n_sensors + (1:n_sensors); % ��ȡ��Ӧά�ȵ�����
        new_pos(idx) = max(new_pos(idx), bounds(i, 1));
        new_pos(idx) = min(new_pos(idx), bounds(i, 2));
    end
end

% Ŀ�꺯����������֪ģ�ͱƽ�������
function coverage = coverage_metric(sensor_positions, n_sensors, sensing_radius, bounds)
    %���룺sensor_positions 100x3�ľ��󣬱�ʾ100���������
    %���룺n_sensors û�õ�
    %���룺sensing_radius �뾶
    %���룺bounds �߽�������0��60��
    
    select_num = 1000;%ѡȡ1000������Ϊ�ο��㣬Խ�󾫶�Խ�ߣ�����о��ٶ�������ȡ500
    matrix = init_jiadianji(select_num,3,bounds(1,2),bounds(1,1));%�üѵ㼯����ȡ��
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
