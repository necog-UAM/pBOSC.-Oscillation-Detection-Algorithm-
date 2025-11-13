%% Script to find the distance from each voxel to MEG sensors and select 3 sources at a linearly increasing distance
% Input are voxel indices in 1925 grid
% Output are voxel indices of 3423 grid and their distances (mm) to nearest
% MEG sensor.

function simvoxs3423 = Sim_voxels(voxels)

[source] = Omega_source(1925, 1.5);

Nvoxin = length(source.inside);
% Minimun distance from each voxel to each sensor
mindist  = zeros(1,length(source.inside));
minidx  = zeros(1,length(source.inside));
for vx = 1:Nvoxin
    eucdist = pdist2(source.forward.grid.pos(source.inside(vx),:), source.forward.grad.chanpos);
    tmpidx(vx) = min(eucdist);
end

[mindist,distidx] = sort(tmpidx,'ascend');

% Find 3 voxels in a line with increasing distances

for vv = 1:length(voxels)
    vx(vv) = find(distidx == voxels(vv)); 
end

% Inices of simulated voxels in the 6804 grid
simvoxs = source.inside(distidx(vx)); 

% Corresponding indices in the 3423 grid
[source] = Omega_source(3423, 1.5);

for ii = 1:length(simvoxs)
    simvoxs3423(ii) = find(source.inside==simvoxs(ii));
end

figure, plot3(source.forward.grid.pos(source.inside,1), source.forward.grid.pos(source.inside,2), source.forward.grid.pos(source.inside,3),'.' )
hold on
plot3(source.forward.grad.chanpos(:,1), source.forward.grad.chanpos(:,2), source.forward.grad.chanpos(:,3),'.' )
plot3( source.forward.grid.pos(source.inside(simvoxs3423),1), source.forward.grid.pos(source.inside(simvoxs3423),2),...
    source.forward.grid.pos(source.inside(simvoxs3423),3),'.', 'MarkerSize',25)

for sc = 1:length(simvoxs3423)
display(['Source: ' num2str(simvoxs3423(sc)) ' Distance (mm): ' num2str(mindist(vx(sc)))]);
end

end