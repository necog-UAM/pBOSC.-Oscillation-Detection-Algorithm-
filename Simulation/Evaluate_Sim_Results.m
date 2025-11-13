%% Evaluate results of the simulations

p.scripts = fileparts(fileparts(matlab.desktop.editor.getActiveFilename)); % Store the path of this script
addpath(genpath(p.scripts))

p.pBOSC = fileparts(p.scripts);
p.data = fullfile(p.pBOSC, 'Simresults');
p.files = fullfile(p.pBOSC, 'Files');
p.simresults =  fullfile(p.pBOSC, 'Simresults');

% Toolbox paths
addpath('Z:\Toolbox\fieldtrip-20230118') 
ft_defaults
addpath('Z:\Toolbox\')

%% Simulation paramteres 

load([p.files '\voxelsdistance.mat'])% Depth is based on Omega_voxdepth.m
Nsources = simvoxs3423 
Nfreqs =  [5 10 20];
Ncycles = [3 10 20];
Nsnr = [0.5 1 1.5];
iterations = 5;

%% Check model results

nmods = length(Nsnr)^4; % Number of models (27*3)
resTPsim = zeros(nmods,iterations);
resTPalt = zeros(nmods,iterations);
resFPsim = zeros(nmods,iterations);
resFPalt = zeros(nmods,iterations);
distAlt = zeros(nmods,iterations);
resTPtotal = zeros(nmods,iterations);
resFPtotal = zeros(nmods,iterations);
misloc = NaN(nmods,iterations);
for md = 1:nmods
    for it = 1:iterations
        cd([p.simresults '\Repetition_' num2str(it)])
        tmpfil = ls('modelorig_*');
        load(tmpfil(md,:));
        resTPsim(md,it) = str2num(model.simvoxel{4,2});
        resTPalt(md,it) = str2num(model.altvoxel{4,2});
        resFPsim(md,it) = str2num(model.simvoxel{2,2});
        resFPalt(md,it) = str2num(model.altvoxel{2,2});
        % calculate distances from alt voxel if true positives are better
        if str2num(model.simvoxel{4,2}) <= str2num(model.altvoxel{4,2})
            distAlt(md,it) = str2num(model.altvoxel{9,2});
            if distAlt(md,it) <= 1.5
                resTPtotal(md,it) = str2num(model.altvoxel{4,2});
                resFPtotal(md,it) = str2num(model.altvoxel{2,2});
            else
                resTPtotal(md,it) = str2num(model.simvoxel{4,2});
                resFPtotal(md,it) = str2num(model.simvoxel{2,2});
                misloc(md,it) = str2num(model.altvoxel{4,2});
            end
        else
            distAlt(md,it) = 0;
            resTPtotal(md,it) = str2num(model.simvoxel{4,2});
            resFPtotal(md,it) = str2num(model.simvoxel{2,2});
        end
    end
    md
end

 save([p.simresults '\modelresults.mat'] ,'resTPtotal' ,'resFPtotal')
    
%% Results
idx = strfind(tmpfil(1,:),'_'); % idx to find the model numbers

firstvoxidx = find(tmpfil(:,idx+1)=='1');
secondvoxidx = find(tmpfil(:,idx+1)=='2');
thirdvoxidx = find(tmpfil(:,idx+1)=='3');

    firstvoxTP = mean(mean(resTPtotal(firstvoxidx,:)));
    secondvoxTP = mean(mean(resTPtotal(secondvoxidx,:)));
    thirdvoxTP = mean(mean(resTPtotal(thirdvoxidx,:)));
    firstvoxFP = mean(mean(resFPtotal(firstvoxidx,:)));
    secondvoxFP = mean(mean(resFPtotal(secondvoxidx,:)));
    thirdvoxFP = mean(mean(resFPtotal(thirdvoxidx,:)));

firstfrexidx = find(tmpfil(:,idx+2)=='1');
secondfrexidx = find(tmpfil(:,idx+2)=='2');
thirdfrexidx = find(tmpfil(:,idx+2)=='3');

    firstfrexTP = mean(mean(resTPtotal(firstfrexidx,:)));
    secondfrexTP = mean(mean(resTPtotal(secondfrexidx,:)));
    thirdfrexTP = mean(mean(resTPtotal(thirdfrexidx,:)));
    firstfrexFP = mean(mean(resFPtotal(firstfrexidx,:)));
    secondfrexFP = mean(mean(resFPtotal(secondfrexidx,:)));
    thirdfrexFP = mean(mean(resFPtotal(thirdfrexidx,:)));

firstcycidx = find(tmpfil(:,idx+3)=='1');
secondcycidx = find(tmpfil(:,idx+3)=='2');
thirdcycidx = find(tmpfil(:,idx+3)=='3');

    firstcycTP = mean(mean(resTPtotal(firstcycidx,:)));
    secondcycTP = mean(mean(resTPtotal(secondcycidx,:)));
    thirdcycTP = mean(mean(resTPtotal(thirdcycidx,:)));
    firstcycFP = mean(mean(resFPtotal(firstcycidx,:)));
    secondcycFP = mean(mean(resFPtotal(secondcycidx,:)));
    thirdcycFP = mean(mean(resFPtotal(thirdcycidx,:)));

firstsnridx = find(tmpfil(:,idx+4)=='1');
secondsnridx = find(tmpfil(:,idx+4)=='2');
thirdsnridx = find(tmpfil(:,idx+4)=='3');

    firstsnrTP = mean(mean(resTPtotal(firstsnridx,:)));
    secondsnrTP = mean(mean(resTPtotal(secondsnridx,:)));
    thirdsnrTP = mean(mean(resTPtotal(thirdsnridx,:)));
    firstsnrFP = mean(mean(resFPtotal(firstsnridx,:)));
    secondsnrFP = mean(mean(resFPtotal(secondsnridx,:)));
    thirdsnrFP = mean(mean(resFPtotal(thirdsnridx,:)));

    %% Found episodes at >1.5 cm
for k = 1:size(distAlt,1)
    for j = 1:size(distAlt,2)
        if distAlt(k,j)<=1.5 & resTPtotal(k,j)>=100
        correctloc(k,j) = 1;
        else
        end
    end
end

sum(correctloc(:)) ./ prod(size(correctloc)) * 100 % Hits 
sum(correctloc(:)==0) ./ prod(size(correctloc)) * 100 % Misses
sum(misloc(:)>0) ./ prod(size(misloc)) * 100 % Found episodes at <1.5cm

distances = distAlt(misloc==100);
Nbins = 1.5:12.5;
[N, edges] = histcounts(distances,Nbins)
percbins = N/sum(~isnan(distances)) *100 % en %

% Figure: Found episodes at >1.5 cm
figure, histogram(distances, Nbins, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', [.6 .6 .6]);
hold on;
[f, xi] = ksdensity(distances);
plot(xi, f, 'k-', 'LineWidth', 2);
xlabel(['Distance (cm)'])
xticks([])

% Localization error at each SNR scenario
sum(sum(distAlt(firstsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100
sum(sum(distAlt(secondsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100
sum(sum(distAlt(thirdsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100

% Localization error 
% 1st SNR
sum(sum(distAlt(firstsnridx,:)>1.5)) / prod(size(distAlt(firstsnridx,:)))*100
% 2nd SNR
sum(sum(distAlt(secondsnridx,:)>1.5)) / prod(size(distAlt(secondsnridx,:)))*100
% 3rd SNR
sum(sum(distAlt(thirdsnridx,:)>1.5)) / prod(size(distAlt(thirdsnridx,:)))*100
 
%% Analysis: Anovas

    firstvoxTP_all = mean(resTPtotal(firstvoxidx,:),2);
    secondvoxTP_all = mean(resTPtotal(secondvoxidx,:),2);
    thirdvoxTP_all = mean(resTPtotal(thirdvoxidx,:),2);

    firstfrexTP_all = mean(resTPtotal(firstfrexidx,:),2);
    secondfrexTP_all = mean(resTPtotal(secondfrexidx,:),2);
    thirdfrexTP_all = mean(resTPtotal(thirdfrexidx,:),2);

    firstcycTP_all = mean(resTPtotal(firstcycidx,:),2);
    secondcycTP_all = mean(resTPtotal(secondcycidx,:),2);
    thirdcycTP_all = mean(resTPtotal(thirdcycidx,:),2);

    firstsnrTP_all = mean(resTPtotal(firstsnridx,:),2);
    secondsnrTP_all = mean(resTPtotal(secondsnridx,:),2);
    thirdsnrTP_all = mean(resTPtotal(thirdsnridx,:),2);

% Voxels
    [p,t,voxStats] = kruskalwallis([firstvoxTP_all secondvoxTP_all thirdvoxTP_all])
    SS_between = t{2,2};
    SS_within = t{3,2}; 
    eta_p_squared = SS_between / (SS_between + SS_within)
    figure; boxplot([firstvoxTP_all secondvoxTP_all thirdvoxTP_all], [{'1'},{'2'},{'3'}])
    title('Voxels'), ylim([-10 110])
    nbins = 0:1.5:max(max(distAlt));
    [postvoxStats,MEANS,H,] = multcompare(voxStats)
% Frex
    [p,t,frexStats] = kruskalwallis([firstfrexTP_all secondfrexTP_all thirdfrexTP_all])
    SS_between = t{2,2};
    SS_within = t{3,2};
    eta_p_squared = SS_between / (SS_between + SS_within)
    figure; boxplot([firstfrexTP_all secondfrexTP_all thirdfrexTP_all], [{'1'},{'2'},{'3'}])
    title('Frequency'), ylim([-10 110])
    [postfrexStats] = multcompare(frexStats)
% Cycles
    [p,t,cycStats] = kruskalwallis([firstcycTP_all secondcycTP_all thirdcycTP_all])
    SS_between = t{2,2};
    SS_within = t{3,2};
    eta_p_squared = SS_between / (SS_between + SS_within)
    figure; boxplot([firstcycTP_all secondcycTP_all thirdcycTP_all], [{'1'},{'2'},{'3'}])
    title('Cycles'), ylim([-10 110])
    [postcycStats] = multcompare(cycStats)
% SNR
    [p,t,snrStats] = kruskalwallis([firstsnrTP_all secondsnrTP_all thirdsnrTP_all])
    SS_between = t{2,2};
    SS_within = t{3,2};
    eta_p_squared = SS_between / (SS_between + SS_within)
    figure; boxplot([firstsnrTP_all secondsnrTP_all thirdsnrTP_all], [{'1'},{'2'},{'3'}])
    title('SNR'), ylim([-10 110])
    [postsnrStats] = multcompare(snrStats)
% Repetitions
    [p,t,iterStats] = kruskalwallis([resTPtotal(:,1) resTPtotal(:,2) resTPtotal(:,3) resTPtotal(:,4) resTPtotal(:,5)])
    figure; boxplot([resTPtotal(:,1) resTPtotal(:,2) resTPtotal(:,3) resTPtotal(:,4) resTPtotal(:,5)], [{'1'},{'2'},{'3'}, {'4'},{'5'}])

    [postiterStats] = multcompare(iterStats)

    % figure; violinplot([resTPtotal(:,1) resTPtotal(:,2) resTPtotal(:,3) resTPtotal(:,4) resTPtotal(:,5)] )

    % Plots
figure, plot([firstvoxTP secondvoxTP thirdvoxTP], '.-', 'MarkerSize',20),
hold on,
plot([firstfrexTP secondfrexTP thirdfrexTP], '.-', 'MarkerSize',20),
plot([firstcycTP secondcycTP thirdcycTP], '.-','MarkerSize',20),
plot([firstsnrTP secondsnrTP thirdsnrTP], '.-','MarkerSize',20),
xlim([0.5 3.5]), ylim([0 100]), legend({'Voxel', 'Frex', 'Cycle', 'SNR'})
title('TP of 4 variables')

%iterations mean
figure, plot(mean(resTPtotal,1), '.-', 'MarkerSize',15)
hold on, yline(mean(mean(resTPtotal)), 'r-')
xlim([0.5 5.5]), ylim([0 100])
title('Repetitions'), xticks(1:5)

% % INdividual var plots
% figure,
% subplot(221),hold on
% tmp = resTPtotal(firstvoxidx,:); tmp = mean(tmp); firstvoxerrorsTP = std(tmp(:)) ./ length(tmp); 
% tmp = resTPtotal(secondvoxidx,:); tmp = mean(tmp); secondvoxerrorsTP = std(tmp(:))./ length(tmp);  
% tmp = resTPtotal(thirdvoxidx,:); tmp = mean(tmp); thirdvoxerrorsTP = std(tmp(:))./ length(tmp);  
% errorbar([firstvoxTP, secondvoxTP, thirdvoxTP], [firstvoxerrorsTP, secondvoxerrorsTP, thirdvoxerrorsTP], 'k-.'),
% xlim([0.5 3.5]), , title('Voxels (depth)')
% 
% subplot(222)
% tmp = resTPtotal(firstfrexidx,:); tmp = mean(tmp);firstfrexerrorsTP = std(tmp(:))./ length(tmp); 
% tmp = resTPtotal(secondfrexidx,:); tmp = mean(tmp);secondfrexerrorsTP = std(tmp(:))./ length(tmp);  
% tmp = resTPtotal(thirdfrexidx,:); tmp = mean(tmp);thirdfrexerrorsTP = std(tmp(:))./ length(tmp); 
% errorbar([firstfrexTP, secondfrexTP, thirdfrexTP], [firstfrexerrorsTP, secondfrexerrorsTP, thirdfrexerrorsTP], 'k-.'),
% xlim([0.5 3.5]), title('Frequencies')
% 
% subplot(223)
% tmp = resTPtotal(firstcycidx,:); tmp = mean(tmp);firstcycerrorsTP = std(tmp(:))./ length(tmp); 
% tmp = resTPtotal(secondcycidx,:); tmp = mean(tmp);secondcycerrorsTP = std(tmp(:))./ length(tmp);  
% tmp = resTPtotal(thirdcycidx,:); tmp = mean(tmp);thirdcycerrorsTP = std(tmp(:))./ length(tmp); 
% errorbar([firstcycTP, secondcycTP, thirdcycTP], [firstcycerrorsTP, secondcycerrorsTP, thirdcycerrorsTP], 'k-.'),
% xlim([0.5 3.5]), title('Cycles')
% 
% subplot(224)
% tmp = resTPtotal(firstsnridx,:); tmp = mean(tmp);firstsnrerrorsTP = std(tmp(:))./ length(tmp);  
% tmp = resTPtotal(secondsnridx,:); tmp = mean(tmp);secondsnrerrorsTP = std(tmp(:))./ length(tmp);  
% tmp = resTPtotal(thirdsnridx,:); tmp = mean(tmp);thirdsnrerrorsTP = std(tmp(:))./ length(tmp);  
% errorbar([firstsnrTP, secondsnrTP, thirdsnrTP], [firstsnrerrorsTP, secondsnrerrorsTP, thirdsnrerrorsTP], 'k-.'),
% xlim([0.5 3.5]), title('SNR')
% 
% % FP
% figure,
% subplot(221),hold on
% tmp = resFPtotal(firstvoxidx,:); tmp = mean(tmp);firstvoxerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(secondvoxidx,:); tmp = mean(tmp);secondvoxerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(thirdvoxidx,:); tmp = mean(tmp);thirdvoxerrorsFP = std(tmp(:))./ length(tmp);  
% errorbar([firstvoxFP, secondvoxFP, thirdvoxFP], [firstvoxerrorsFP, secondvoxerrorsFP, thirdvoxerrorsFP], 'r-.'),
% xlim([0.5 3.5]), title('Voxels (depth)')
% 
% subplot(222)
% tmp = resFPtotal(firstfrexidx,:);tmp = mean(tmp); firstfrexerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(secondfrexidx,:); tmp = mean(tmp);secondfrexerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(thirdfrexidx,:); tmp = mean(tmp);thirdfrexerrorsFP = std(tmp(:))./ length(tmp);  
% errorbar([firstfrexFP, secondfrexFP, thirdfrexFP], [firstfrexerrorsFP, secondfrexerrorsFP, thirdfrexerrorsFP], 'r-.'),
% xlim([0.5 3.5]),  title('Frequencies')
% 
% subplot(223)
% tmp = resFPtotal(firstcycidx,:);tmp = mean(tmp); firstcycerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(secondcycidx,:);tmp = mean(tmp); secondcycerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(thirdcycidx,:);tmp = mean(tmp); thirdcycerrorsFP = std(tmp(:))./ length(tmp);  
% errorbar([firstcycFP, secondcycFP, thirdcycFP], [firstcycerrorsFP, secondcycerrorsFP, thirdcycerrorsFP], 'r-.'),
% xlim([0.5 3.5]), title('Cycles')
% 
% subplot(224)
% tmp = resFPtotal(firstsnridx,:); tmp = mean(tmp);firstsnrerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(secondsnridx,:); tmp = mean(tmp);secondsnrerrorsFP = std(tmp(:))./ length(tmp);  
% tmp = resFPtotal(thirdsnridx,:); tmp = mean(tmp);thirdsnrerrorsFP = std(tmp(:))./ length(tmp);  
% errorbar([firstsnrFP, secondsnrFP, thirdsnrFP], [firstsnrerrorsFP, secondsnrerrorsFP, thirdsnrerrorsFP], 'r-.'),
% xlim([0.5 3.5]),   title('SNR')


%% Hits matrix
mtxsnr1 = zeros(3,3);
ct = [1,4,7];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr1(mtxr,mtxc) = round(mean(mean(resTPtotal(tmpidx,:),2)));
        init = init+9;
    end
end
mtxsnr2 = zeros(3,3);
ct = [2,5,8];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr2(mtxr,mtxc) = round(mean(mean(resTPtotal(tmpidx,:),2)));
        init = init+9;
    end
end
mtxsnr3 = zeros(3,3);
ct = [3,6,9];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr3(mtxr,mtxc) = round(mean(mean(resTPtotal(tmpidx,:),2)));
        init = init+9;
    end
end

load RdBu_R.mat % in Files

xlabels = [{'5'}; {'10'}; {'20'}];ylabels = [{'3'}; {'10'}; {'20'}];
figure('Position', [0 100 1000 250])
subplot(131),h =heatmap(xlabels,ylabels,mtxsnr1), h.ColorbarVisible = 'off',
colormap(RdBu_r), h.FontSize = 14
caxis([0 100])
subplot(132),h =heatmap(xlabels,ylabels,mtxsnr2), h.ColorbarVisible = 'off',
colormap(RdBu_r), h.FontSize = 14 
caxis([0 100])
subplot(133),h =heatmap(xlabels,ylabels,mtxsnr3), h.ColorbarVisible = 'off',
colormap(RdBu_r) ,h.FontSize = 14
caxis([0 100])

% print(gcf, ['Z:\OMEGA\Enrique\figures\TPorig.tiff'], '-dtiff', '-r300');

% MEANS
% total
mean(mean(mtxsnr1))
mean(mean(mtxsnr2))
mean(mean(mtxsnr3))
% ciclos
[mean(mtxsnr1,2) mean(mtxsnr2,2) mean(mtxsnr3,2)]
% frex
[mean(mtxsnr1,1)' mean(mtxsnr2,1)' mean(mtxsnr3,1)']




%% False Positives Matrix
mtxsnr1 = zeros(3,3);
ct = [1,4,7];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr1(mtxr,mtxc) = round(mean(mean(resFPtotal(tmpidx,:),2)),3);
        init = init+9;
    end
end
mtxsnr2 = zeros(3,3);
ct = [2,5,8];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr2(mtxr,mtxc) = round(mean(mean(resFPtotal(tmpidx,:),2)),3);
        init = init+9;
    end
end
mtxsnr3 = zeros(3,3);
ct = [3,6,9];
for mtxr = 1:3
    init = ct(mtxr);
    for mtxc = 1:3 
        tmpidx = init:27:81;
        mtxsnr3(mtxr,mtxc) = round(mean(mean(resFPtotal(tmpidx,:),2)),3);
        init = init+9;
    end
end

load RdBu_R.mat
xlabels = [{'5'}; {'10'}; {'20'}];ylabels = [{'3'}; {'10'}; {'20'}];
figure('Position', [0 100 1000 250])
subplot(131),h =heatmap(xlabels,ylabels,mtxsnr1), h.ColorbarVisible = 'off',
colormap(RdBu_r), h.FontSize = 14
caxis([0 100])
subplot(132),h =heatmap(xlabels,ylabels,mtxsnr2), h.ColorbarVisible = 'off',
colormap(RdBu_r), h.FontSize = 14
caxis([0 100])
subplot(133),h =heatmap(xlabels,ylabels,mtxsnr3), h.ColorbarVisible = 'off',
colormap(RdBu_r) ,h.FontSize = 14
caxis([0 100])
% print(gcf, ['Z:\OMEGA\Enrique\figures\FPorig.tiff'], '-dtiff', '-r300');

max(max(mtxsnr1))
max(max(mtxsnr2))
max(max(mtxsnr3))
