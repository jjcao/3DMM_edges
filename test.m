clc;clearvars;close all;
%addpath ../jjcao_code/toolbox/jjcao_plot;
addpath ../jjcao_code/toolbox/jjcao_interact;
addpath ../jjcao_code/toolbox/jjcao_mesh;
%addpath(genpath('../'));
%%
%save test.mat
%save('t1.mat','FV', 'im', 'R', 't', 's')

load t1.mat;
badTextThre = 0.1;
edges = compute_edges(FV.faces);

%%
rotpts = R*FV.vertices';
xdif = rotpts(1,edges(1,:)) - rotpts(1,edges(2,:));
ydif = rotpts(2,edges(1,:)) - rotpts(2,edges(2,:));
zdif = rotpts(3,edges(1,:)) - rotpts(3,edges(2,:));

dist2 = xdif.^2 + ydif.^2;
dist3 = dist2 + zdif.^2;
ratio = dist2./dist3;
pos = ratio<badTextThre;

idx = zeros(size(FV.vertices,1),1);
idx(edges(1,pos)) = 1;
idx(edges(2,pos)) = 1;

%%
if(~isa(im,'double'))
    im=double(im)./255;
end

FV.facevertexcdata = faceTexture(FV,R,t,s,im);

rotpts = R*FV.vertices';
FV.vertices = rotpts';
pts = FV.vertices(idx>0,:);
% pos1 = isnan(FV.facevertexcdata(:,1));
% %[0.594704890760164,0.438963497906832,0.323837899519381]
% FV.facevertexcdata(pos1,1)=0.594704890760164;
% FV.facevertexcdata(pos1,2)=0.438963497906832;
% FV.facevertexcdata(pos1,3)=0.323837899519381;

%%
figure;
p = patch(FV, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); axis equal; axis off; p.FaceColor = 'interp';
view3d rot; hold on;
scatter3(pts(:,1),pts(:,2),pts(:,3),100,'.','MarkerEdgeColor','r'); 


%%
figure; 
subplot(1,2,1); 
p = patch(FV, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); light; axis equal; axis off; p.FaceColor = 'interp';
%patch(FV, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong');light; axis equal; axis off;
colorbar
subplot(1,2,2); imshow(renderFace(FV,im,R,t,s,false));
% subplot(1,3,3); imshow(renderFace(FV,im,R,t,s,true));