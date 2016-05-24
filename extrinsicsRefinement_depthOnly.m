%% extrinsicsRefinement_depthOnly.m

%Note this is not right

load HERB_projectionData.mat

camera_size = [640, 480]
plane_coeffs = [-0.00080255,  -0.863825,  -0.503791,  0.550167]

tag_cam = [423   354     1];

tag_world = [ -0.0753773175125 0.137904712863 0.916102764935]
tag_world_ori = [ -0.000611770725047 0.85704194419 -0.515219139985 0.00528861317285]

%% 
addpath('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/matlab_code/matpcl')

tbl_cam_xyz = loadpcd('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/dataFrom20160520_173121bag/pcd_files/planecloud_1.pcd');
tbl_cam_xyz = tbl_cam_xyz';

 a =plane_coeffs(1) ;
 b =plane_coeffs(2) ;
 c =plane_coeffs(3) ;
 d =plane_coeffs(4) ;
 
 t = [-a*d, -b*d, -c*d]';
 
  dprod = dot(t,[0,0,1]')

  if ( dprod < 0)
    t = -1.0 * t;
  end

  z = t;
  %try to align the x axis with the x axis of the original frame
  %or the y axis if z and x are too close too each other
  x = [1, 0, 0];
  if ( abs(dot(z, x)) > 1.0 - 1.0e-4) 
      x = [0, 1, 0];
  end
  
  y = cross(z, x);
  x = cross(y, z);

  R =zeros(3,3);
  R(:, 1) = x; 	
  R(:, 2) = y; 	
  R(:, 3) = z; 	
    
tbl_cam_xyz = [tbl_cam_xyz ones(size(tbl_cam_xyz,1), 1)];
    
  %% convert to camera points
  
  tbl_cam_uv = (A * tbl_cam_xyz')';
  tbl_cam_uv = [tbl_cam_uv(:,1)./tbl_cam_uv(:,3) tbl_cam_uv(:,2)./tbl_cam_uv(:,3) ones(size(tbl_cam_uv,1), 1)];
  
    %% convert to world points 
  
  T = [R t; 0 0 0 1];
  tabl_xyz_world = (T * tbl_cam_xyz')';
  
  
  %%
  
 rgb_img = imread('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/dataFrom20160520_173121bag/color/frame0000.jpg');
imshow(rgb_img)

 dpt_img = imread('/home/apirespi/Documents/Thesis/ExtrinsicsCalibration/dataFrom20160520_173121bag/depth/frame0000.jpg');
figure;imshow(dpt_img, [])