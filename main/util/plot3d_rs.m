function[] = plot3d_rs(camera_cell, param)

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
    
    plot3(camera_cell{1, cell_id}.point_state_vector(:,1), camera_cell{1, cell_id}.point_state_vector(:,2), camera_cell{1, cell_id}.point_state_vector(:,3), camera_cell{1, cell_id}.color + "*");
    hold on
end
% [Rot_Rows, trans_Rows] = linearEgoMotion(camera_set(image_id).gt_oritation, camera_set(image_id).gt_translation,pram.w_gt,pram.d_gt,cameraParameters('IntrinsicMatrix',pram.K_matrix'));
% 
% %% Simulated RS porjection  
% p2d = RSWorld2Image(p3d, Rot_Rows, trans_Rows, cameraParameters('IntrinsicMatrix',pram.K_matrix'));
% [Rot_Rows, trans_Rows] = linearEgoMotion(camera_set(image_id).gt_oritation, camera_set(image_id).gt_translation,[0;0;0],[0;0;0],cameraParameters('IntrinsicMatrix',pram.K_matrix'));
% 
% %% Simulated RS porjection  
% p2d_gs = RSWorld2Image(p3d, Rot_Rows, trans_Rows, cameraParameters('IntrinsicMatrix',pram.K_matrix'));
% 
% %% Create 2D RS image 
% RSImage = figure;
% I = zeros(pram.H, pram.W);
% imshow(I, [1 1 1]);
% plot(p2d(:,1),p2d(:,2),'g*');
% hold on;
% plot(p2d_gs(:,1),p2d_gs(:,2),'r*');
hold off
end