function plot_flow_component(xvec_speed, yvec_speed, xvec_pressure, yvec_pressure, u1_vec, u2_vec, p_vec, style)
    % 支持速度和压力使用不同网格的流场绘图函数（已移除等高线）
    % 输入参数：
    %   xvec_speed, yvec_speed - 速度场的x、y方向一维坐标向量
    %   xvec_pressure, yvec_pressure - 压力场的x、y方向一维坐标向量
    %   u1_vec, u2_vec - 速度x、y分量的一维向量（对应速度网格）
    %   p_vec - 压力的一维向量（对应压力网格）
    %   style - 绘图风格：'general'（矩形域）或 'channel'（管道流）
    
    % 中文显示设置
    set(0, 'DefaultTextFontName', 'SimHei');
    set(0, 'DefaultAxesFontName', 'SimHei');
    
    % 默认风格处理
    if nargin < 8 || isempty(style)
        style = 'general';
    end
    
    %% 1. 处理速度场数据（独立网格）
    % 速度网格维度计算
    nx_speed = size(xvec_speed, 2);
    ny_speed = size(yvec_speed, 1);
    total_nodes_speed = nx_speed * ny_speed;
    % 生成速度场二维网格
    [x_speed, y_speed] = meshgrid(xvec_speed(1,:), yvec_speed(:,1));
    % 检查速度数据长度并转换为二维
    if length(u1_vec) ~= total_nodes_speed || length(u2_vec) ~= total_nodes_speed
        error('速度向量长度与速度网格点数不匹配（应为%d）', total_nodes_speed);
    end
    u1 = reshape(u1_vec, ny_speed, nx_speed);
    u2 = reshape(u2_vec, ny_speed, nx_speed);
    
    %% 2. 处理压力场数据（独立网格）
    % 压力网格维度计算
    nx_pressure = size(xvec_pressure, 2);
    ny_pressure = size(yvec_pressure, 1);
    total_nodes_pressure = nx_pressure * ny_pressure;
    % 生成压力场二维网格
    [x_pressure, y_pressure] = meshgrid(xvec_pressure(1,:), yvec_pressure(:,1));
    % 检查压力数据长度并转换为二维
    if length(p_vec) ~= total_nodes_pressure
        error('压力向量长度与压力网格点数不匹配（应为%d）', total_nodes_pressure);
    end
    p = reshape(p_vec, ny_pressure, nx_pressure);
    
    %% 3. 根据风格绘制合并窗口
    switch lower(style)
        case 'general'
            % 普通矩形域：2行2列布局
            figure('Position', [150, 150, 1000, 800]);
            sgtitle('速度分量（独立网格）与压力分布（独立网格）', 'FontSize', 14, 'FontWeight', 'bold');
            
            % 1. 速度x分量（使用速度网格）
            subplot(2, 2, 1);
            contourf(x_speed, y_speed, u1, 30, 'LineStyle', 'none');
            cb1 = colorbar;  % 先创建颜色条对象
            cb1.Label.String = '速度值';  % 再设置标签（兼容所有版本）
            cb1.FontSize = 9;
            title('速度x方向分量 (u₁)', 'FontSize', 12);
            xlabel('x坐标'); ylabel('y坐标');
            axis equal tight; box on;
            
            % 2. 速度y分量（使用速度网格）
            subplot(2, 2, 2);
            contourf(x_speed, y_speed, u2, 30, 'LineStyle', 'none');
            cb2 = colorbar;
            cb2.Label.String = '速度值';
            cb2.FontSize = 9;
            title('速度y方向分量 (u₂)', 'FontSize', 12);
            xlabel('x坐标'); ylabel('y坐标');
            axis equal tight; box on;
            
            % 3. 压力分布（使用压力网格）
            subplot(2, 2, [3, 4]);
            contourf(x_pressure, y_pressure, p, 30, 'LineStyle', 'none');
            cb3 = colorbar;
            cb3.Label.String = '压力值';
            cb3.FontSize = 9;
            title('压力分布 (p)', 'FontSize', 12);
            xlabel('x坐标'); ylabel('y坐标');
            axis equal tight; box on;
            
            tight_layout(pad=1.2);
            
        case 'channel'
            % 管道流：3行1列布局
            figure('Position', [100, 200, 1200, 350]);
            sgtitle('管道流速度（独立网格）与压力（独立网格）分布', 'FontSize', 14, 'FontWeight', 'bold');
            cmap = jet(64);  % 统一配色
            
            % 1. 速度x分量（速度网格）
            subplot(3,1, 1);
            contourf(x_speed, y_speed, u1, 30, 'LineStyle', 'none');
            colormap(cmap);
            cb1 = colorbar;
            cb1.Label.String = 'u1值';
            cb1.FontSize = 8;
            title('速度x分量', 'FontSize', 11);
            axis equal tight; box on;
            set(gca, 'XTickLabel', '', 'YTickLabel', '', 'TickDir', 'out');
            
            % 2. 速度y分量（速度网格）
            subplot(3,1, 2);
            contourf(x_speed, y_speed, u2, 30, 'LineStyle', 'none');
            colormap(cmap);
            cb2 = colorbar;
            cb2.Label.String = 'u2值';
            cb2.FontSize = 8;
            title('速度y分量', 'FontSize', 11);
            axis equal tight; box on;
            set(gca, 'XTickLabel', '', 'YTickLabel', '', 'TickDir', 'out');
            
            % 3. 压力分布（压力网格）
            subplot(3,1, 3);
            contourf(x_pressure, y_pressure, p, 30, 'LineStyle', 'none');
            colormap(cmap);
            cb3 = colorbar;
            cb3.Label.String = 'p值';
            cb3.FontSize = 8;
            title('压力分布', 'FontSize', 11);
            axis equal tight; box on;
            set(gca, 'XTickLabel', '', 'YTickLabel', '', 'TickDir', 'out');
    end
end
