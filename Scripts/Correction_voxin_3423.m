% Correction to make all inside voxels (the 1925 voxels) have 19 neighbors:
% total 3423 voxels.

% 1. Load the inside 1925 voxels
load('Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\correccion_vox_inside_10mm.mat')
load('Z:\OMEGA\OMEGA_data\sub-0001\ses-0001\source_inverse_10mm.mat')
voxel_inside_1925 = find(inside);

neigh_inside=[];
for vi=1:length(voxel_inside_1925)
        Dst = pdist2(source.pos(voxel_inside_1925(vi),:),source.pos);
        neigh = find(Dst<=1.5)';
        if length(neigh) ~= 19
            disp('error')
            break
        end
        neigh_inside = [neigh_inside;neigh];
end
        voxel_inside_3423 = unique(neigh_inside);
        inside = zeros(length(inside),1);
        inside(voxel_inside_3423) = 1;
        inside = logical(inside);

save('Z:\OMEGA\correccion_vox_inside_10mm_3423.mat', "inside")