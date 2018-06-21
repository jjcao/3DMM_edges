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

inputFile = 'inputImages/fface1.png';

rgb = imread(inputFile);
[I,map] = rgb2ind(rgb, 65536);
[X,Y] = meshgrid(-100:100,-80:80);
Z = -(X.^2 + Y.^2); Z(:) = 1;

figure; warp(X,Y,Z,I,map);
axis equal; view3d rot; hold on;