% -------------------------------------------------------------------------
% 线性变换的体素史诗
% -------------------------------------------------------------------------

% --- 1. 设定 ---
clc; clear; close all;

D = diag([2, 0.5, -1.5]);
[V, ~] = qr(randn(3)); 
V_inv = V';
A = V * D * V_inv;
I = eye(3); 

x = [1.5; 1; 2]; % 待变换的向量

fprintf('--- 变换设定 ---\n');
fprintf('复杂变换矩阵 A:\n'); disp(A);
fprintf('自然坐标轴 (特征向量) V:\n'); disp(V);
fprintf('动画准备就绪。这是一个无限循环动画。\n');
fprintf('请关闭图形窗口或在命令窗口按 Ctrl+C 来终止程序。\n');
pause(1);

% --- 2. 颜色与几何的关键帧 ---
color_I = [0.7, 0.7, 0.7]; 
color_A = [0.4, 0.9, 0.4]; 
color_V = [0.9, 0.9, 0.4]; 

T1_Identity = I;
T2_A_Warped = A;
T3_Eigen_Aligned = V;

% 向量x的路径关键帧
x1_start = x;
x2_midpoint = A * x;
x3_endpoint = x; 

% --- 3. 搭建舞台 ---
fig = figure('Name', '向量与空间的同步演化史诗', 'Position', [50, 50, 1200, 900], 'Color', 'w');
ax = axes('Parent', fig, 'Position', [0.05 0.05 0.9 0.85]);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
start_azimuth = 220;
start_elevation = 25;
view(ax, [start_azimuth, start_elevation]);
axis_limit = 5;
axis(ax, [-1 1 -1 1 -1 1] * axis_limit);
set(ax, 'FontSize', 12, 'CameraViewAngle', 7);
xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z');

% 静态标题
h_textbox = annotation('textbox', [0.05, 0.85, 0.4, 0.1], ...
    'String', '准备阶段: 初始标准体素空间 (灰色)', ...
    'FontSize', 16, 'EdgeColor', 'none', 'BackgroundColor', 'w', 'FaceAlpha', 0.5);

% --- 4. 定义并绘制初始体素网格 ---
grid_range = -4:2:4;
[X_base, Y_base, Z_base] = meshgrid(grid_range);
num_voxels = numel(X_base);
unit_cube_corners = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]'-0.5;
faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 4 3 2; 5 8 7 6];
initial_voxel_vertices = zeros(3, 8, num_voxels);
for i = 1:num_voxels
    offset = [X_base(i); Y_base(i); Z_base(i)];
    initial_voxel_vertices(:,:,i) = unit_cube_corners*2 + offset;
end

% --- 5. 绘制会动的“主角” ---
% 初始化动态体素网格
h_patches_active = create_grid(ax, initial_voxel_vertices, color_I, 0.2);

% 初始化动态基向量
h_basis_i = quiver3(ax, 0,0,0, I(1,1), I(2,1), I(3,1), 0, 'r', 'LineWidth', 4, 'MaxHeadSize', 0.5);
h_basis_j = quiver3(ax, 0,0,0, I(1,2), I(2,2), I(3,2), 0, 'g', 'LineWidth', 4, 'MaxHeadSize', 0.5);
h_basis_k = quiver3(ax, 0,0,0, I(1,3), I(2,3), I(3,3), 0, 'b', 'LineWidth', 4, 'MaxHeadSize', 0.5);

% 初始化动态向量x及其轨迹
h_vector_active = quiver3(ax, 0,0,0, x(1), x(2), x(3), 0, 'm', 'LineWidth', 4, 'MaxHeadSize', 0.3);
h_trajectory = plot3(ax, x(1), x(2), x(3), 'm--', 'LineWidth', 2.5);
path_data = []; % 轨迹点将在循环中生成

drawnow;
pause(2);

% --- 6. 动画开始 (无限循环) ---
num_steps_per_act = 150;
total_step_counter = 0;

while ishandle(fig)
    
    path_data = x; % 每个循环重置轨迹
    set(h_trajectory, 'XData', path_data(1,:), 'YData', path_data(2,:), 'ZData', path_data(3,:));

    % --- 创建并冻结标准基(I)的幽灵 ---
    h_patches_ghost_I = create_grid(ax, initial_voxel_vertices, color_I*0.95, 0.05, '--');
    
    % --- 第一幕: I (灰色) -> A (绿色) ---
    for step = 1:num_steps_per_act
        if ~ishandle(fig), return; end
        s = step / num_steps_per_act;
        T_current = (1-s)*T1_Identity + s*T2_A_Warped;
        color_current = (1-s)*color_I + s*color_A;
        set(h_textbox, 'String', sprintf('第一幕: 标准空间 -> A扭曲空间 (%.0f%%)', s*100));
        
        x_current = (1-s)*x1_start + s*x2_midpoint;
        
        update_all_visuals(h_patches_active, {h_basis_i, h_basis_j, h_basis_k}, h_vector_active, h_trajectory, ...
                           T_current, color_current, initial_voxel_vertices, x_current, path_data);
        path_data = [path_data, x_current];
        
        total_step_counter = total_step_counter + 1;
        view(ax, start_azimuth + total_step_counter*0.1, start_elevation);
        drawnow limitrate;
    end
    delete(h_patches_ghost_I); % 删除灰色幽灵
    pause(1);

    % --- 创建并冻结A基的幽灵 (体素和向量) ---
    final_vertices_A = T2_A_Warped * reshape(initial_voxel_vertices, 3, []);
    h_patches_ghost_A = create_grid(ax, final_vertices_A, color_A, 0.05, '--');
    h_vector_ghost_A = quiver3(ax, 0,0,0, x2_midpoint(1), x2_midpoint(2), x2_midpoint(3), 0, 'Color', color_A, 'LineWidth', 2, 'LineStyle', '--');
    pause(1.5);

    % --- 第二幕: A (绿色) -> V (黄色) ---
    for step = 1:num_steps_per_act
        if ~ishandle(fig), return; end
        s = step / num_steps_per_act;
        T_current = (1-s)*T2_A_Warped + s*T3_Eigen_Aligned;
        color_current = (1-s)*color_A + s*color_V;
        set(h_textbox, 'String', sprintf('第二幕: 扭曲空间 -> 对齐至特征基 (%.0f%%)', s*100));
        
        x_current = (1-s)*x2_midpoint + s*x3_endpoint;

        update_all_visuals(h_patches_active, {h_basis_i, h_basis_j, h_basis_k}, h_vector_active, h_trajectory, ...
                           T_current, color_current, initial_voxel_vertices, x_current, path_data);
        path_data = [path_data, x_current];
        
        total_step_counter = total_step_counter + 1;
        view(ax, start_azimuth + total_step_counter*0.1, start_elevation);
        drawnow limitrate;
    end
    delete(h_patches_ghost_A); delete(h_vector_ghost_A); % 删除绿色幽灵
    pause(1.5);
    
    % --- 第三幕: 重置 V (黄色) -> I (灰色) ---
    set(h_textbox, 'String', '重置循环...');
    for step = 1:num_steps_per_act
        if ~ishandle(fig), return; end
        s = step / num_steps_per_act;
        T_current = (1-s)*T3_Eigen_Aligned + s*T1_Identity;
        color_current = (1-s)*color_V + s*color_I;

        update_all_visuals(h_patches_active, {h_basis_i, h_basis_j, h_basis_k}, h_vector_active, h_trajectory, ...
                           T_current, color_current, initial_voxel_vertices, x, path_data); % 向量x保持不动
        
        total_step_counter = total_step_counter + 1;
        view(ax, start_azimuth + total_step_counter*0.1, start_elevation);
        drawnow limitrate;
    end
    pause(1);
end

% -------------------------------------------------------------------------
% --- 本地函数定义 ---
% -------------------------------------------------------------------------

% 创建一个网格并返回句柄
function h_patches = create_grid(ax, vertices, color, face_alpha, edge_style)
    if nargin < 5, edge_style = '-'; end
    if nargin < 4, face_alpha = 0.1; end
    
    num_v = size(vertices, 3);
    h_patches = gobjects(num_v, 1);
    faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 4 3 2; 5 8 7 6];
    
    for i = 1:num_v
        % *** 此处是修正点: EdgeLineStyle -> LineStyle ***
        h_patches(i) = patch(ax, 'Vertices', vertices(:,:,i)', 'Faces', faces, ...
                                 'FaceColor', color, 'EdgeColor', 'k', 'LineStyle', edge_style, ...
                                 'FaceAlpha', face_alpha, 'LineWidth', 0.5);
    end
end

% 统一更新所有动态视觉元素
function update_all_visuals(h_patches, h_basis, h_vector, h_traj, T, color, initial_verts, x_vec, path)
    % 1. 更新体素网格
    current_all_verts = T * reshape(initial_verts, 3, []);
    verts_reshaped = permute(reshape(current_all_verts, 3, 8, []), [2, 1, 3]);
    for idx = 1:length(h_patches)
        set(h_patches(idx), 'Vertices', verts_reshaped(:,:,idx), 'FaceColor', color);
    end
    % 2. 更新基向量
    set(h_basis{1}, 'UData', T(1,1), 'VData', T(2,1), 'WData', T(3,1));
    set(h_basis{2}, 'UData', T(1,2), 'VData', T(2,2), 'WData', T(3,2));
    set(h_basis{3}, 'UData', T(1,3), 'VData', T(2,3), 'WData', T(3,3));
    % 3. 更新向量x和轨迹
    set(h_vector, 'UData', x_vec(1), 'VData', x_vec(2), 'WData', x_vec(3));
    set(h_traj, 'XData', path(1,:), 'YData', path(2,:), 'ZData', path(3,:));
end