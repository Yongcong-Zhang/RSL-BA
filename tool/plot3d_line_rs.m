function[] = plot3d_line_rs(camera_cell, param)

camera_cell_size = size(camera_cell, 2);
Scene3DFigure = figure;

for cell_id = 1:camera_cell_size
    for i = 1:param.num_camera
        id = (i - 1) * 12 + 1;
        
        gt_r_state_vector = camera_cell{1, cell_id}.camera_state_vector(id : id + 2);
        gt_rotation = axang2rotm([gt_r_state_vector / norm(gt_r_state_vector), norm(gt_r_state_vector)]);
        gt_translation = camera_cell{1, cell_id}.camera_state_vector(id + 3: id + 5)';
    
        cam = plotCamera('Location',  -inv(gt_rotation) * gt_translation,'Orientation',(gt_rotation),'Opacity',0, 'Size', 1, 'color',camera_cell{1, cell_id}.color );
        hold on
    
    end
    grid on
    axis equal
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    % psd_gs = gt_camera_param(id : end);
    
    for j = 1: 4
        line_orth =  camera_cell{1, cell_id}.line_state_vector(j,:);
        line_plucker = orth_to_plucker(line_orth);
        p1 = line_plucker*[1; 0; 0; -2.5];
        p2 = line_plucker*[1; 0; 0; 2.5];
        P1 = p1./p1(4); P1 = P1(1:3,1);
        P2 = p2./p2(4); P2 = P2(1:3,1);
        plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], camera_cell{1, cell_id}.color, 'LineWidth', 1);
        hold on;
    end
    for j = 5: 8
        line_orth =  camera_cell{1, cell_id}.line_state_vector(j,:);
        line_plucker = orth_to_plucker(line_orth);
        p1 = line_plucker*[0; 1; 0; -2.5];
        p2 = line_plucker*[0; 1; 0; 2.5];
        P1 = p1./p1(4); P1 = P1(1:3,1);
        P2 = p2./p2(4); P2 = P2(1:3,1);
        plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], camera_cell{1, cell_id}.color, 'LineWidth', 1);
        hold on;
    end
    for j = 9: 12
        line_orth =  camera_cell{1, cell_id}.line_state_vector(j,:);
        line_plucker = orth_to_plucker(line_orth);
        p1 = line_plucker*[0; 0; 1; -2.5];
        p2 = line_plucker*[0; 0; 1; 2.5];
        P1 = p1./p1(4); P1 = P1(1:3,1);
        P2 = p2./p2(4); P2 = P2(1:3,1);
        plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], camera_cell{1, cell_id}.color, 'LineWidth', 1);
        hold on;
    end

    hold on
end
hold off
end