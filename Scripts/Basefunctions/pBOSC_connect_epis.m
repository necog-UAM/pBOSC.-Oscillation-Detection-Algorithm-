function epis2 = pBOSC_connect_epis(cfg, epis)

epis2 = cell(1,size(epis,2));
Nvox = size(epis,2);
fsample = cfg.fsample;
for v=1:Nvox   % only the 1925 voxels inside the cortex to save time
    if ~isempty(epis{v})

        listepis = [1:length(epis{v})];           % original episodes, 1st iteration
        checkepis = [];
        epis2{v}(1) = epis{v}(1);
        Nit = 0;             % number of iterations
        while ~isempty(listepis)
            if Nit>0
                listepis = checkepis;
                checkepis = [];
                epis{v} = epis2{v};
                epis2{v}(1) = epis{v}(1);
            end

            ep1 = 1;
            ctep = 1;
            while ep1 <= length(listepis)
                epis{v}(ep1).timeps = single(epis{v}(ep1).timeps);
                t0 = epis{v}(ep1).timeps(1);
                cycntps = 1/epis{v}(ep1).freq*fsample;       % number of time points corresponding to 1 cycle
                tf = round(epis{v}(ep1).timeps(end) + 1/2*cycntps);     % + 1/2 missing cycle
                tep1 = [t0:tf];
                if ep1+1 <= length(epis{v}) & epis{v}(ep1).freq ~= 100
                    for ep2 = ep1+1:length(epis{v})
                        epis{v}(ep2).timeps = single(epis{v}(ep2).timeps);
                        t0 = round(epis{v}(ep2).timeps(1) - 1/2*cycntps);   % - 1/2 missing cycle
                        cycntps = 1/epis{v}(ep2).freq*fsample;   % number of time points corresponding to 1 cycle
                        tf = epis{v}(ep2).timeps(end);
                        tep2 = [t0:tf];
                        if isempty(intersect(tep1,tep2))
                            ct = 1;
                            break;
                        elseif ~isempty(intersect(tep1,tep2))
                            f10 = (epis{v}(ep1).freq-0.1*epis{v}(ep1).freq);
                            f1f = (epis{v}(ep1).freq+0.1*epis{v}(ep1).freq);
                            f1 = [round(f10,1):0.1:round(f1f,1)];
                            f20 = (epis{v}(ep2).freq-0.1*epis{v}(ep2).freq);
                            f2f = (epis{v}(ep2).freq+0.1*epis{v}(ep2).freq);
                            f2 = [round(f20,1):0.1:round(f2f,1)];
                            if isempty(intersect(f1,f2))
                                ct = 1;
                            elseif ~isempty(intersect(f1,f2))
                                epis2{v}(ctep).freq = (epis{v}(ep1).freq*length(epis{v}(ep1).timeps) + epis{v}(ep2).freq*length(epis{v}(ep2).timeps)) ./ (length(epis{v}(ep1).timeps) + length(epis{v}(ep2).timeps));   % weighted average of freqs
                                epis2{v}(ctep).timeps = unique([tep1 tep2]);
                                epis2{v}(ctep).dur_sec = length(epis2{v}(ctep).timeps)./fsample;     % duration of the new connected episode
                                epis2{v}(ctep).dur_cyc = epis2{v}(ctep).freq .* epis2{v}(ctep).dur_sec;
                                epis2{v}(ctep).power = (epis{v}(ep1).power*length(epis{v}(ep1).timeps) + epis{v}(ep2).power*length(epis{v}(ep2).timeps)) ./ (length(epis{v}(ep1).timeps) + length(epis{v}(ep2).timeps));   % weighted average of power
                                if ep1 == 1 & length(listepis)==1      % only 1st episode is a candidate
                                    epis2{v}(ep2)= [] ;
                                end
                                ct = 0;
                                break;
                            end
                        end
                    end
                else
                    ct=1;
                end
                if ct==0
                    checkepis = [checkepis ep1];      % new episodes might have new connections; run the algorithm again
                    ctep = ctep + 1;
                    ep1 = ep1+1;
                    epis{v}(ep2).freq = 100;           % mark episodes that have already been connected
                elseif ct==1
                    epis2{v}(ctep) = epis{v}(ep1);
                    ctep = ctep + 1;
                    ep1 = ep1+1;
                end
            end
            Nit = Nit+1;
        end
    epis2{v} = epis2{v}([epis2{v}.freq] ~= 100);            % remove epis marked with freq = 100
    end
end


