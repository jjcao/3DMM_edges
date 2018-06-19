% ex_render_texture_back
% the projected texture contains artifacts 暂不知道原因。
% 还是直接在3D做blending方便。

clc;clearvars;close all;
addpath utils;
addpath ../jjcao_code/toolbox/jjcao_interact;
addpath ../jjcao_code/toolbox/jjcao_mesh;
addpath ../jjcao_code/toolbox/jjcao_mesh/feature;
addpath ../jjcao_code/toolbox/kdtree;

inputFile = 'image_0018_texture.mat';
load(['output/' inputFile]);

oglp.height=size(im,1);
oglp.width=size(im,2);
oglp.i_amb_light = [1 1 1];
oglp.i_dir_light = [1 1 1];
oglp.i_dir_light = [0 0 0];

Rr = R;
Rr(4,4)=1;
Sr = eye(4).*s;
Tr = eye(4);
Tr(1:2,4)=t;
%T = Tr*Sr*Rr;
T = Tr*Sr;
clear Tr Sr Rr

renderim = render_texture_back(FV, T, oglp, zeros(size(im)));
imshow(renderim);


