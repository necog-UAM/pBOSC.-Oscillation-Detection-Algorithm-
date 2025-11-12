  [~, ~,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5);
    Nvoxin = length(epis);
    Nfrex = length(frex);
    Ntp = size(dataclean.trial{1},2);
    time = dataclean.time{1};

    episodes = int8(zeros(Nvoxin,Nfrex,Ntp)); % Srate in episodes was 128, so now adapt to 256
    for v=1:Nvoxin
        for ep=1:length(epis{v})
            fm = epis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis{v}(ep).timeps;
          
            episodes(v,fbin,tpts) = 1;
        end
    end


epis5d = int8(zeros(1,1,1,Nfrex,Ntp));
epis5d = repmat(epis5d,[dim,1]);


for vx = 1:size(episodes,1)
    epis5d(xx(vx),yy(vx),zz(vx),:,:) = episodes(vx,:,:);
end

conn = conndef(2,'maximal');            % salen demasiados elementos conectados en el cluster 1; probar poniendo un thr arriba

CC = bwconncomp(epis5d,conn);

% Get the labeled matrix
% Preallocate output label matrix
labeled5d = zeros(CC.ImageSize, 'uint16'); % Use uint16 or uint32 depending on the number of objects

% Fill in the labels manually
for k = 1:CC.NumObjects
    labeled5d(CC.PixelIdxList{k}) = true;
end



