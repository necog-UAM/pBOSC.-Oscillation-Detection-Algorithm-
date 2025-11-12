% Los subs 43 (0049) y 113 (0170) tienen un newepis que opcupa más de 2GB,
% por lo que singlepis no se guarda. 

for s = 1:Nsub
    s
    cd([dpath '\sub-' subs{s} '\ses-' sess{s}])
    load newepis.mat
    epis2 = cell(size(epis));  % Preallocate cell array for epis2
    for ep = 1:length(epis)
        fnams = fieldnames(epis{ep});
        epis2{ep} = repmat(struct(), size(epis{ep}));  % Preallocate structure array
        for sz = 1:length(epis{ep})
            for nams = 1:length(fnams)
                if nams == 4
                    epis2{ep}(sz).(fnams{nams}) = int32(epis{ep}(sz).(fnams{nams}));
                else
                    epis2{ep}(sz).(fnams{nams}) = round(single(epis{ep}(sz).(fnams{nams})),3);
                end
            end
        end
    end
    epis = epis2;
    tmp = whos("epis");
    if (tmp.bytes *1e-9) >= 2
        warning(['Sub ' num2str(subs{s}) ' > 2GB'])
        save singlepis epis Ntp
    else
        save singlepis epis
    end
end

