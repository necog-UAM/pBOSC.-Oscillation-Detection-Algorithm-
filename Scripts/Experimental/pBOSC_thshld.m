function thshld = Omega_thshld(signal, noise, conf)
tic
% define variables
delta_pc = 0.5; % percentil jumps
delta_smt = 2/delta_pc; %smooth window size
thshld = single(zeros(1,size(signal,2)));
thshld = repmat(thshld, [size(signal,1) 1]);

for vx=1:size(signal,1)

    sign_pc = prctile(squeeze(signal(vx,:,:)),[0:delta_pc:100],2);
    cntrs_tmp = prctile(squeeze(signal(vx,:,:)),[delta_pc/2:delta_pc:100],2);
    smth_signal = [];
    smth_noise = [];
    cntrs = [] ;
    ct=1;

    tmpsignal = sort(squeeze(signal(vx,:,:)),2);
    tmpnoise = sort(squeeze(noise(vx,:,:)),2);
  
    for ii = 1:size(cntrs_tmp,2)
        try
            smth_signal(:,ct) = sum(tmpsignal>=sign_pc(:,ii-delta_smt) & tmpsignal<=sign_pc(:,ii+delta_smt),2);
            smth_noise(:,ct)  = sum(tmpnoise>=sign_pc(:,ii-delta_smt) & tmpnoise<=sign_pc(:,ii+delta_smt),2);
            cntrs(:,ct) =  cntrs_tmp(:,ii);
            ct=ct+1;
        end
    end

    snratio = (smth_signal - smth_noise)./(smth_signal + smth_noise);
    [r,c] = find(snratio>=conf); % row and column
    [C,~,IC] = unique(r);

    %% To get idx
    idx = zeros(size(signal,2),1) + size(snratio,2);
    for row = 1:length(C)
        tmpidx = find(IC==row);
        idx(C(row)) = c(tmpidx(1)); %la 1º vez que aparece esa fx en la row, que col le corresponde
    end

    thshld(vx,:) = diag(cntrs(:,idx));

end

tend = toc;
display(['Total time was ' num2str(tend/60) ' minutes'])

end
