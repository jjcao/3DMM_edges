% jjcao @ 2018

clc;clearvars;close all;
MYTOOLBOXROOT='../jjcao_code/toolbox/';
%MYTOOLBOXROOT='E:/jjcao_code/toolbox/';
addpath ([MYTOOLBOXROOT 'jjcao_mesh'])
addpath ([MYTOOLBOXROOT 'jjcao_mesh/feature'])
addpath ([MYTOOLBOXROOT 'jjcao_io'])
addpath ([MYTOOLBOXROOT 'kdtree'])
addpath ([MYTOOLBOXROOT 'jjcao_interact'])
addpath ([MYTOOLBOXROOT 'jjcao_common'])
addpath ([MYTOOLBOXROOT 'jjcao_math'])
DEBUG=1;

%% input
inputFile = 'output/fface1_texture';%test_LFW1,image_0018,fface1,sface1
load([inputFile '.mat']);

inputImg = 'inputImages/fface1.png';
rgb = imread(inputImg);
rgb = flip(rgb,1);
[I,map] = rgb2ind(rgb, 65536);

[X,Y] = meshgrid(1:size(I,2),1:size(I,1));
Z = -(X.^2 + Y.^2); Z(:) = 1;

%%
figure; warp(X,Y,Z,I,map);
axis equal; view3d rot; hold on;

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

%%
p = patch('Faces', FV.faces, 'Vertices', V2, 'FaceVertexCData', FV.facevertexcdata, 'EdgeColor', 'none'); 
axis equal; axis off; 
p.FaceColor = 'interp';
view3d rot; hold on;
