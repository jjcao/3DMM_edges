% jjcao @ 2018
%%
clc;clearvars;close all;
addpath utils;
addpath ../jjcao_code/toolbox/jjcao_interact;
addpath ../jjcao_code/toolbox/jjcao_mesh;
addpath ../jjcao_code/toolbox/jjcao_mesh/feature;
addpath ../jjcao_code/toolbox/jjcao_io;

badTextThre = 0.2;
inputFile = 'image_0018.mat';%test_LFW1,image_0018,fface1,sface1
%%
load(['output/' inputFile]);

%% find vertices with poor texture, by analysing inner product between vertex normal and z axis
FVr = FV;
FVr.vertices = (R*FV.vertices')';
vnormal = compute_normal(FVr.vertices,FVr.faces, 1)';
z = repmat([0,0,1], size(FVr.vertices,1),1);
angle = dot(vnormal, z, 2);
pos = abs(angle)<badTextThre;

idx = false(size(FV.vertices,1),1);
idx(pos) = true;
idx_via_normal = idx;

pts = FVr.vertices(idx_via_normal,:);
figure; patch(FVr); axis equal; view3d rot; hold on;
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
title('via normal')
%% display original texture
if(~isa(im,'double'))
    im=double(im)./255;
end

FVr.facevertexcdata = faceTexture(FV,R,t,s,im);
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('original texture')

%% repair texture via normal
FVr.facevertexcdata = faceTexture(FV,R,t,s,im);
FVr.facevertexcdata(idx_via_normal,:) = nan;
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none','FaceLighting', 'phong'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('repair via normal');
%light;


