% jjcao @ 2018

clc;clearvars;close all;
%MYTOOLBOXROOT='../jjcao_code/toolbox/';
MYTOOLBOXROOT='E:/jjcao_code/toolbox/';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh/feature'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'kdtree'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])
DEBUG=1;
SHOW_IMG = 0;

%% input
inputFile = 'output/fface1_texture';%test_LFW1,image_0018,fface1,sface1
load([inputFile '.mat']);

if SHOW_IMG
    inputImg = 'inputImages/fface1.png';
    rgb = imread(inputImg);
    rgb = flip(rgb,1);
    [I,map] = rgb2ind(rgb, 65536);

    [X,Y] = meshgrid(0:size(I,2),0:size(I,1));
    Z = -(X.^2 + Y.^2); 
    Z(:) = 1;

    figure; warp(X,Y,Z,I,map);
    axis equal; view3d rot; hold on;
end 
%%
Rr = R;
Rr(4,4)=1;
Sr = eye(4).*s;
Tr = eye(4);
Tr(1:2,4)=t;
T = Tr*Sr*Rr;

% Get the extrinsic transformation matrix
M = T(1: 3, :); % the output does not need to be in homogeneous coordinates

% Get the vertices
V           = FV.vertices;
Nvertices   = size(FV.vertices, 1);

% Compute the transformed vertices
V(:, 4)	= 1;        % use homogeneous coordinates for input
V2   	= V * M.';	% the vertices are transposed

%% normalization 
cmean = mean(V2);
V2 = V2 - repmat(cmean, size(V2,1),1);
if SHOW_IMG
    X = X - cmean(1); Y = Y - cmean(2); Z = Z - cmean(3);
end

%% rotation
k = 20;
filenames = 'rota.gif';
if SHOW_IMG
    figure; warp(X,Y,Z,I,map);
end
set(gcf,'color','w'); hold on;axis equal; axis off; view3d rot;
for i = 1:4:k
    theta = (i-1)/4*pi;
    T = [cos(theta) 0 sin(theta) 0;0 1 0 0;-sin(theta) 0 cos(theta) 0;0 0 0 1];
    points = [V2 ones(size(V2,1),1)]*T;
    points = points(:,1:3);
    
    p = patch('Faces', FV.faces, 'Vertices', V2, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none');     
    %p = patch('Faces', FV.faces, 'Vertices', points, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none');     
    p.FaceColor = 'interp';

%      f = getframe(gcf);
%         imind = frame2im(f);
%         [imind,cm] = rgb2ind(imind,256);
%         if j == 1
%             imwrite(imind,cm,filenames,'gif', 'Loopcount',inf);
%         else
%             imwrite(imind,cm,filenames,'gif','WriteMode','overwrite');
%         end 
    delete(p);
end
%%
% p = patch('Faces', FV.faces, 'Vertices', V2, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); 
% axis equal; axis off; 
% p.FaceColor = 'interp';
% view3d rot; hold on;