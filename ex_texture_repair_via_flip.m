% jjcao @ 2018

clc;clearvars;close all;
addpath utils;
addpath ../jjcao_code/toolbox/jjcao_interact;
addpath ../jjcao_code/toolbox/jjcao_mesh;
addpath ../jjcao_code/toolbox/jjcao_mesh/feature;
addpath ../jjcao_code/toolbox/kdtree;
%addpath(genpath('../'));
%%
%save test.mat
%save('t1.mat','FV', 'im', 'R', 't', 's')

badTextThre = 0.2;
negativeThre = 25000;%1000

%%
load image_0018.mat;
%load yangmi.mat;
edges = compute_edges(FV.faces);

%% find symmetric vertices of vertices with x < 0
% idx_negative_sym(i) is the symmetric vertex index if FV.vertices(i,1) < 0;
% idx_negative_sym(i) is zero if FV.vertices(i,1) > 0;
tree = kdtree_build( FV.vertices );

idx_postive = FV.vertices(:,1)>=-negativeThre;
idx_negative = FV.vertices(:,1)<-negativeThre;
pts = FV.vertices(idx_negative,:);
pts(:,1) = - pts(:,1);
[idx_negative_sym,dst] = kdtree_nearest_neighbor(tree,pts);

idx = zeros(size(FV.vertices,1),1);
idx(idx_negative) = idx_negative_sym;
idx_negative_sym = idx;
pts = FV.vertices(idx_negative_sym(idx_negative),:);

figure;patch(FV); axis equal; view3d rot; hold on;
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 
pts = FV.vertices(idx_negative,:);
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 

kdtree_delete(tree);

%% find vertices with poor texture & negative x, by analysing ratio between 3D edges and the projected 2d edges
% idx(i) == 1 means that texture of ith vertex is poor & netative x

rotpts = R*FV.vertices';
xdif = rotpts(1,edges(1,:)) - rotpts(1,edges(2,:));
ydif = rotpts(2,edges(1,:)) - rotpts(2,edges(2,:));
zdif = rotpts(3,edges(1,:)) - rotpts(3,edges(2,:));

dist2 = xdif.^2 + ydif.^2;
dist3 = dist2 + zdif.^2;
ratio = dist2./dist3;
pos = ratio<badTextThre;

idx = false(size(FV.vertices,1),1);
idx(edges(1,pos)) = true;
idx(edges(2,pos)) = true;
idx(idx_postive) = false;

FVr = FV;
FVr.vertices = rotpts';
pts = FVr.vertices(idx,:);
figure; patch(FVr); axis equal; view3d rot; hold on;
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
pts = FVr.vertices(idx_negative_sym(idx),:);
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 

%% find vertices with poor texture & negative x, by analysing inner product between vertex normal and z axis
% 这样找到的idx不稀疏了，但是为啥修复后的texture有这么明显的光照问题？
normal = compute_normal(FVr.vertices,FVr.faces, 1)';
z = repmat([0,0,1], size(FVr.vertices,1),1);
angle = dot(normal, z, 2);
pos = abs(angle)<0.5;

idx = false(size(FV.vertices,1),1);
idx(pos) = true;
idx(idx_postive) = false;

FVr = FV;
FVr.vertices = rotpts';
pts = FVr.vertices(idx,:);
figure; patch(FVr); axis equal; view3d rot; hold on;
%scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','r'); 
pts = FVr.vertices(idx_negative_sym(idx),:);
scatter3(pts(:,1),pts(:,2),pts(:,3),10,'.','MarkerEdgeColor','g'); 

%% 对称一半是对的，而对称一部分点，是错的；说明是idx计算有问题；
% idx的顶点在脑门和嘴角比原始顶点稀疏，是edges算得不对？ 
% 改用旋转后的顶点的法向和z轴的夹角，有提高，但还是不自然。
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

%%
FVr.facevertexcdata(idx_negative,:) = FVr.facevertexcdata(idx_negative_sym(idx_negative),:)*1.2;
figure;
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
% pts = FVr.vertices(idx_negative_sym(idx),:);
% scatter3(pts(:,1),pts(:,2),pts(:,3),100,'.','MarkerEdgeColor','g'); 
% pts = FVr.vertices(idx,:);
% scatter3(pts(:,1),pts(:,2),pts(:,3),100,'.','MarkerEdgeColor','r'); 


%%
figure; 
subplot(1,2,1); 
p = patch(FVr, 'FaceVertexCData', FVr.facevertexcdata, 'EdgeColor', 'none'); light; axis equal; axis off; p.FaceColor = 'interp';
%patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong');light; axis equal; axis off;
colorbar
subplot(1,2,2); imshow(renderFace(FVr,im,R,t,s,true));
% subplot(1,3,3); imshow(renderFace(FV,im,R,t,s,true));