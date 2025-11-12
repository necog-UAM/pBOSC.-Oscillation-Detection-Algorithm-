vox_in=find(source.inside==1)
for vi=1:1925
    for vo=1:length(source.inside)
        D = pdist2(source.pos(vox_in(vi),:),source.pos(vo,:))
        if D<=1.5
            vo

        end
    end
end
unique