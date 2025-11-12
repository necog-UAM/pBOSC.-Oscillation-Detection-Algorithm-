function conepis = connect_episES(epis, Ntp)


[source_forward, source,voxel_inside, dim,xx,yy,zz,connmat,dtempl] = Omega_voxin(1925, 1.5);
[frex, Nfrex, fb_lab,fb1,fb2] = Omega_frex();

Nvoxin = length(voxin);

    Nvoxin = length(epis);
    Nfrex = length(frex);

    episodes = zeros(Nvoxin, Nfrex, Ntp);
    for v=1:Nvoxin
        for ep=1:length(epis{v})
            fm = epis{v}(ep).freq;
            fbin = dsearchn(frex',fm);
            tpts = epis{v}(ep).timeps;
            episodes(v,fbin,tpts) = 1;
        end
    end
for v = 1:size(episodes,1)
    for fx = 1:size(episodes,2)
        fepis = find(squeeze(episodes(v,fx,:)));
        if ~isempty(fepis)
            d = diff(fepis);
            breakidx = find(d>1); % when the episode lasts
            vxnb = find(connmat(v,:));
            fxnb = fx-1:fx+1;
            for brks = 1:length(breakidx)
                [voxrow, fxcol] = find(squeeze(episodes(vxnb,fxnb,breakidx(brks)+1))); % if break point is continued in nb chnn and fx
                if ~isempty(voxrow)
                    vxnb(voxrow)
                    fxnb(fxcol)
                else

                end
            end

        end
        





idx = find(connmat(vx,:))
tmpepis = squeeze(mean(mean(episodes(idx,18-1:18+1,:),2),1));
tmpepis(tmpepis>0) = 1;
figure, imagesc(tmpepis)