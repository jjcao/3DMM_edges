% jjcao @ 2018

clc;clearvars;close all;
%MYTOOLBOXROOT='../../jjcao_code/toolbox/';
MYTOOLBOXROOT='E:/jjcao_code/toolbox/';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh/feature'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'kdtree'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])
addpath utils;
%addpath(genpath('../'));
%%
badTextThre = 0.2;
negativeThre = 5000;%25000%15000;%0
stitchBelt = 25000;
inputFile = 'output/image_0018';%test_LFW1,image_0018,fface1,sface1;
DEBUG=1;
%%
load([inputFile '.mat']);
edges = compute_edges(FV.faces);

%% find symmetric vertices of vertices with x < 0
tree = kdtree_build( FV.vertices );

idx_postive = FV.vertices(:,1)>=-negativeThre;
idx_negative = FV.vertices(:,1)<=-negativeThre;

pts = FV.vertices(idx_negative,:);
pts(:,1) = - pts(:,1);
[idx_negative_sym,dst] = kdtree_nearest_neighbor(tree,pts);

idx = zeros(size(FV.vertices,1),1);
idx(idx_negative) = idx_negative_sym;
idx_negative_sym = idx;
pts = FV.vertices(idx_negative_sym(idx_negative),:);

if DEBUG
    figure;patch(FV); axis equal; view3d rot; hold on;
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
    pts = FV.vertices(idx_negative,:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
end
    
kdtree_delete(tree);

%% display original texture
if(~isa(im,'double'))
    im=double(im)./255;
end

FV.facevertexcdata = faceTexture(FV,R,t,s,im);
figure;title('original texture')
p = patch(FV, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;

save([inputFile '_texture.mat'],'FV','im', 'R', 't', 's');

%% repair texture via x coordinate
FVr = FV;
%FVr.facevertexcdata = faceTexture(FV,R,t,s,im);
FVr.facevertexcdata(idx_negative,:) = FVr.facevertexcdata(idx_negative_sym(idx_negative),:)*1.0;
figure;title('repair via coordinate'); 
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;


%% blend
% todo: to hsv, blend hue, then back to rgb.
% todo: non linear blend, gamma correction?
FV_left= FVr; FV_right = FV;
FV_left.facevertexcdata(idx_postive,:) = NaN;
idx_nb = FV.vertices(:,1)<=-negativeThre - stitchBelt;
FV_right.facevertexcdata(idx_nb,:) = NaN;

result = FV_left;
result.facevertexcdata(idx_nb,:) = FV_left.facevertexcdata(idx_nb,:);
result.facevertexcdata(idx_postive,:) = FV_right.facevertexcdata(idx_postive,:);
idx_blend = (result.vertices(:,1) >=-negativeThre - stitchBelt) & (result.vertices(:,1) <=-negativeThre);

xmin = -negativeThre - stitchBelt; xmax = -negativeThre;
weight_right = (result.vertices(idx_blend,1) - xmin)./stitchBelt;
weight_right = repmat(weight_right, 1, 3);
weight_left = 1 - weight_right;

tmp.facevertexcdata = FV_left.facevertexcdata;
tmp.facevertexcdata(idx_blend,:) = FV_left.facevertexcdata(idx_blend,:);%.*weight_left
result.facevertexcdata(idx_blend,:) = weight_left.*FV_left.facevertexcdata(idx_blend,:) + weight_right.*FV_right.facevertexcdata(idx_blend,:);
idx_nan = isnan(result.facevertexcdata(:,1));
idx_blend_nan = idx_nan & idx_blend;
%result.facevertexcdata(idx_blend_nan,:) = tmp.facevertexcdata(idx_blend_nan,:);

figure;title('blend'); 
p = patch(result, 'FaceVertexCData', result.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
if DEBUG % show the blend region
    r1 = result;
    r1.facevertexcdata(idx_nb,:) = NaN;
    r1.facevertexcdata(idx_postive,:) = NaN;
    figure;title('blend region'); scatter3(r1.vertices(:,1),r1.vertices(:,2),r1.vertices(:,3),10,'.','MarkerEdgeColor','g');  hold on; 
    p = patch(r1, 'FaceVertexCData', r1.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp'; view3d rot;
    clear r1;
end
%% repair
%idx1 = FV.vertices(:,3)>=10000;
idx_handle = ~idx_nan;
constraint_id=find(idx_handle);
constraint_value = result.facevertexcdata(idx_handle,:);
idx = ~idx_handle;

if DEBUG
    figure; 
    %patch(FV);
    axis equal; view3d rot; hold on;
    pts = result.vertices(idx,:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
    pts = result.vertices(idx_handle,:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
    title('handles')
end

%% computen Laplacian
options = [];
type = 'conformal';%combinatorial;conformal;%spring
options.use_c_implementation = 1;
tic; 
L = compute_mesh_laplacian(result.vertices,result.faces,type);
toc

% repair
options.solver = 1;
options.method = 'hard';% 'hard', 'soft'
b = zeros(size(L,1),3);

tic;fid = compute_least_square_system(L, b, constraint_id, constraint_value,options);toc

% output & plot
figure;title('repaired texture')
p = patch(result, 'FaceVertexCData', fid, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
