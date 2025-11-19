function [source] = Omega_source(N, distance) 

%--------------------------------------------------------%

% Input arguments:
% N: number of voxels inside the brain. 
% a) 1925 (only inner voxels).
% b) 3294 (center of the head included) but external voxels do not have 19 neighbors.
% c) 3423 (the 1925 inner voxels have all max number of neighbors [19]). 

% distance: in cm^3.
% 1: 7 voxels.
% 1.5: 19 voxels.

% Output arguments:
% source: structure containing the following fields:
% forward: template of a source forward
% inverse: template of a source inverse.
% inside: index of the N inside voxels from the total number of 6804.
% dim: size of data.
% xx,yy,zz: spatial coordinates.
% connmat: logical connectivity matrix of NxN neighbors.
% dtempl: template to organize data.

%--------------------------------------------------------%

% Path
p.repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));


% Load source forward and inverse from a template to get positions
load([p.repo '\Files\Templates\source_forward_10mm_' num2str(N) '.mat'])
load([p.repo '\Files\Templates\source_inverse_10mm_' num2str(N) '.mat'])

% Get indices of the N inside voxels
inside = find(source_inverse.inside);

% Get x y z coordinates
for i=1:3
    pos(:,i) = source_inverse.pos(:,i)-min(source_inverse.pos(:,i))+1;
    dim(i) = max(pos(:,i));
end

% Initialize template of inside voxels (N) x total voxels (6804)
vtempl = zeros(dim);
dtempl = zeros(size(inside,1),prod(dim));

for v=1:size(inside,1)
    temp = zeros(dim);
    p=pos(inside(v),:);
    vtempl(p(1),p(2),p(3)) = v;
    temp(p(1),p(2),p(3)) = 1;
    dtempl(v,:) = temp(:);
end

% Generate connectivity matrix of neighboring voxels based on input distance 
connmat = zeros(size(inside,1),size(inside,1));

for v1=1:size(inside,1)
    vpos1 = source_inverse.pos(inside(v1),:);
    for v2=1:size(inside,1)
        vpos2 = source_inverse.pos(inside(v2),:);
        vdif = norm(vpos2-vpos1);
        if vdif<=distance
            connmat(v1,v2) = 1;
        end
    end
end

connmat = logical(connmat);

% Get x y z coordinates of the N inside voxels
for v=1:size(inside,1)
    ind = find(vtempl==v);
    [i1,i2,i3] = ind2sub(size(vtempl),ind);
    xx(v)=i1;
    yy(v)=i2;
    zz(v)=i3;
end

% Return:
source.forward = source_forward;
source.inverse = source_inverse;
source.inside   = inside;
source.dim            = dim;
source.connmat        = connmat;
source.dtempl         = dtempl;
source.xx             = xx;
source.yy             = yy;
source.zz             = zz;

end