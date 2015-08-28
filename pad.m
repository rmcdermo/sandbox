% McDermott
% 3-3-2013
% pad.m
%
% function [pp] = pad(p)
%
% For use with pcolor and shading faceted such that cell-centered data of
% size(p) may be correctly visualized without white space.

function [pp]=pad(p)

[nx,ny] = size(p);

pp = [p,p(:,ny)];
pp = [pp;pp(nx,:)];
pp(nx+1,ny+1) = p(nx,ny);
