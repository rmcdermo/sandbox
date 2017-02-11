% McDermott
% 2-10-2017
% flameview.m

close all
clear all

% [A,map] = imread('burner_side_1.jpeg');
% [nx,ny,rgb] = size(A)
% image(A)

vidObj = VideoReader('9_psi_short3.mov');

k = 1;
while hasFrame(vidObj)
    s(k).cdata = readFrame(vidObj);
    image(s(k).cdata)
    pause
    k = k+1;
end

