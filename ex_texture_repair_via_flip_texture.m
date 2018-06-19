% jjcao @ 2018

clc;clearvars;close all;
MYTOOLBOXROOT='../../jjcao_code/toolbox/';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh/feature'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'kdtree'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath utils;
%%
badTextThre = 0.2;
negativeThre = 0;%25000%15000;%0
inputFile = 'output/lin';%test_LFW1,image_0018,fface1,sface1;
DEBUG=1;
%%
load([inputFile '.mat']);
edges = compute_edges(FV.faces);

%% find symmetric vertices of vertices with x < 0
% idx_negative_sym(i) is the symmetric vertex index if FV.vertices(i,1) < 0;
% idx_negative_sym(i) is zero if FV.vertices(i,1) > 0;
tree = kdtree_build( FV.vertices );

idx_postive = FV.vertices(:,1)>=-negativeThre;
idx_negative = FV.vertices(:,1)<-negativeThre;
% idx_postive = FV.vertices(:,1)<-negativeThre;
% idx_negative = FV.vertices(:,1)>=-negativeThre;

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

%% find vertices with poor texture & negative x, by analysing ratio between 3D edges and the projected 2d edges
% idx_via_ratio(i) == 1 means that texture of ith vertex is poor & netative x
% note: idx_via_ratio的顶点在脑门和嘴角比原始顶点�?��，缺了很多vertex，是edges算得不对�?还没找到原因�?
% 

% rotpts = R*FV.vertices';
% xdif = rotpts(1,edges(1,:)) - rotpts(1,edges(2,:));
% ydif = rotpts(2,edges(1,:)) - rotpts(2,edges(2,:));
% zdif = rotpts(3,edges(1,:)) - rotpts(3,edges(2,:));
% 
% dist2 = xdif.^2 + ydif.^2;
% dist3 = dist2 + zdif.^2;
% ratio = dist2./dist3;
% pos = ratio<badTextThre;
% 
% idx = false(size(FV.vertices,1),1);
% idx(edges(1,pos)) = true;
% idx(edges(2,pos)) = true;
% idx(idx_postive) = false;
% idx_via_ratio = idx;
% 

% vertices = rotpts';
% pts = vertices(idx_via_ratio,:);
% figure; patch('Faces',FV.faces,'Vertices',vertices); axis equal; view3d rot; hold on;
% scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
% pts = vertices(idx_negative_sym(idx_via_ratio),:);
% scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
% title('ratio between edges')

%% find vertices with poor texture & negative x, by analysing inner product between vertex normal and z axis
% 这样找到的idx_via_normal不稀疏了，但是修复后的texture光照问题明显，不能�?过整体修正改善，必须要做�?��blending
FVr = FV;
FVr.vertices = (R*FV.vertices')';
vnormal = compute_normal(FVr.vertices,FVr.faces, 1)';
z = repmat([0,0,1], size(FVr.vertices,1),1);
angle = dot(vnormal, z, 2);
pos = abs(angle)<0.5;

idx = false(size(FV.vertices,1),1);
idx(pos) = true;
idx(idx_postive) = false;
idx_via_normal = idx;

if DEBUG
    pts = FVr.vertices(idx_via_normal,:);
    figure; patch(FVr); axis equal; view3d rot; hold on;
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
    pts = FVr.vertices(idx_negative_sym(idx_via_normal),:);
    scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
    title('via normal')
end
%% display original texture
if(~isa(im,'double'))
    im=double(im)./255;
end

FVr.facevertexcdata = faceTexture(FV,R,t,s,im);
%FVr.facevertexcdata(idx,:) = FVr.facevertexcdata(idx_negative_sym(idx),:);
% FVr.facevertexcdata(idx,:) = 0.0;
% FVr.facevertexcdata(idx_negative_sym(idx),:) = 0.0;
% FVr.facevertexcdata(idx_negative,:) = 0.0;
% FVr.facevertexcdata(idx_negative_sym(idx_negative),:) = 0.0;
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('original texture')

FV.facevertexcdata = FVr.facevertexcdata;
save([inputFile '_texture.mat'],'FV','im', 'R', 't', 's');

%% repair texture via normal
FVr.facevertexcdata = faceTexture(FV,R,t,s,im);
FVr.facevertexcdata(idx_via_normal,:) = FVr.facevertexcdata(idx_negative_sym(idx_via_normal),:);
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none','FaceLighting', 'phong'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('repair via normal');
%light;

%figure;
% pts = FVr.vertices(idx_negative_sym(idx),:);
% scatter3(pts(:,1),pts(:,2),pts(:,3),100,'.','MarkerEdgeColor','g'); 
% pts = FVr.vertices(idx,:);
% scatter3(pts(:,1),pts(:,2),pts(:,3),100,'.','MarkerEdgeColor','r'); 

%% repair texture via x coordinate
FVr.facevertexcdata(idx_negative,:) = FVr.facevertexcdata(idx_negative_sym(idx_negative),:);
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
title('repair via coordinate'); 
%light;

%%
% Rr = R;
% Rr(4,4)=1;
% Sr = eye(4).*s;
% Tr = eye(4);
% Tr(1:2,4)=t;
% T = Tr*Sr*Rr;
% 
% % Get the extrinsic transformation matrix
% M = T(1: 3, :); % the output does not need to be in homogeneous coordinates
% 
% % Get the vertices
% V           = FV.vertices;
% Nvertices   = size(FV.vertices, 1);
% 
% % Compute the transformed vertices
% V(:, 4)	= 1;        % use homogeneous coordinates for input
% V2   	= V * M.';	% the vertices are transposed
% 
% write_obj('fface1.obj', V2, FV.faces);

%%
figure; 
subplot(1,3,1); 
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); light; axis equal; axis off; p.FaceColor = 'interp';
%patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong');light; axis equal; axis off;
colorbar
subplot(1,3,2); imshow(renderFace(FVr,im,R,t,s,true));
subplot(1,3,3); imshow(renderFace(FV,im,R,t,s,false));