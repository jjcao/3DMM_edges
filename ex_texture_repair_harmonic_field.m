% ex_texture_repair_harmonic_field
%
% 
% Copyright (c) 2018 Junjie Cao
%%
clear;clc;close all;
MYTOOLBOXROOT='../jjcao_code/toolbox/';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])

DEBUG=1;
badTextThre = 0.3;
inputFile = 'output/lin_texture';%test_LFW1,image_0018,fface1,sface1
%% input
load([inputFile '.mat']);
%[verts,faces] = read_mesh([MYTOOLBOXROOT 'fface1.obj']);
%save('fface1_texture.mat','FV','im', 'R', 't', 's');

verts = FV.vertices;
faces = FV.faces;
nverts = size(verts,1);

figure;
p = patch(FV, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('original texture')

%% find vertices with poor texture by analysing inner product between vertex normal and z axis
FVr = FV;
FVr.vertices = (R*FV.vertices')';
vnormal = compute_normal(FVr.vertices,FVr.faces, 1)';
z = repmat([0,0,1], size(FVr.vertices,1),1);
angle = dot(vnormal, z, 2);
pos = abs(angle)<badTextThre;

idx = false(size(FV.vertices,1),1);
idx(pos) = true;
idx_via_normal = idx;

if DEBUG
    pts = FVr.vertices(idx_via_normal,:);
    figure; patch(FVr); axis equal; view3d rot; hold on;
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
    title('via normal')
end

%% find vertices with none-non texture as handle / boundary conditions
idx = isnan(FV.facevertexcdata(:,1));
idx(idx_via_normal) = true;
idx_handle = ~idx;

constraint_id=find(idx_handle);
constraint_value = FV.facevertexcdata(idx_handle,:);
    
if DEBUG
    figure; 
    %patch(FV);
    axis equal; view3d rot; hold on;
    pts = FV.vertices(idx,:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
    pts = FV.vertices(idx_handle,:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
    title('handles')
end

%% computen Laplacian
options = [];
type = 'conformal';%combinatorial;conformal;%spring
options.use_c_implementation = 1;
tic; L = compute_mesh_laplacian(verts,faces,type);toc

%% repair
options.solver = 1;
options.method = 'hard';% 'hard', 'soft'
b = zeros(size(L,1),3);

tic;fid = compute_least_square_system(L, b, constraint_id, constraint_value,options);toc


%% output & plot
figure;title('repaired texture')
p = patch(FV, 'FaceVertexCData', fid, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
