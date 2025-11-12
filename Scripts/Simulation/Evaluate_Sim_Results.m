%% Evaluate simulation results

% Prepare paths 
dpath = 'Z:\OMEGA\Enrique\Results';
addpath('Z:\Toolbox\fieldtrip-20230118')
ft_defaults
addpath('Z:\Toolbox\fBOSC-main')
start_fBOSC
addpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\BaseScripts')
addpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\Simulation')
addpath('Z:\Toolbox\NECOG')
cd(dpath)
addpath('Z:\Toolbox\BrewerMap-master\BrewerMap-master')

%% Simulation paramteres parameters
load Z:\OMEGA\Enrique\Results\voxelsdistance2.mat
Nsources = simvoxs3423%[2957 2405 745]; % Depth is based on Omega_voxdepth.m
Nfreqs =  [5 10 20];
Ncycles = [3 10 20];
Nsnr = [0.5 1 1.5];
iterations = 5;

%% Check model results
nmods = length(Nsnr)^4; % combinatoria. Num de valores de la variable elevado a num de vars
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
        % cd(['Z:\OMEGA\Enrique\Results\Sources\Iteration_' num2str(it)])
        cd(['Z:\OMEGA\Enrique\OMEGA-Repo\Simresults\Repetition_' num2str(it)])
        tmpfil = ls('modelorig2_*');
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

cd(dpath)
% save modelresultsorig2 resTPtotal resFPtotal
    % Results
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

% Anovas

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

%%
% Cuantas veces encuentra epis a <=1.5cm
for k = 1:size(distAlt,1)
    for j = 1:size(distAlt,2)
        if distAlt(k,j)<=1.5 & resTPtotal(k,j)>=100
        correctloc(k,j) = 1;
        else
        end
    end
end
sum(correctloc(:)) ./ prod(size(correctloc)) * 100 % cuantas veces hace bien
sum(correctloc(:)==0) ./ prod(size(correctloc)) * 100 % cuantas mal
sum(misloc(:)>0) ./ prod(size(misloc)) * 100 % cuantas veces encuentra el epis bien pero en otro voxel a mas de 1.5cm
sum(resTPtotal(:)<100) ./ prod(size(resTPtotal)) * 100

distances = distAlt(misloc==100);
Nbins = 1.5:12.5;
[N, edges] = histcounts(distances,Nbins)
percbins = N/sum(~isnan(distances)) *100 % en %

figure, histogram(distances, Nbins, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceAlpha', 0.6, 'FaceColor', [.6 .6 .6]);
hold on;
[f, xi] = ksdensity(distances);
plot(xi, f, 'k-', 'LineWidth', 2);
xlabel(['Distance (cm)'])
xticks([])

% Distances Peor y mejor escenario
sum(sum(distAlt(firstsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100
sum(sum(distAlt(secondsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100
sum(sum(distAlt(thirdsnridx,:)>1.5)) / sum(sum(distAlt>1.5)) * 100
% Escenario Ideal
bestidx = intersect(intersect(firstfrexidx,thirdcycidx), thirdsnridx);
sum(sum(distAlt(bestidx,:)>1.5)) / prod(size(distAlt(bestidx,:))) * 100
% Escenario Terror
worstidx = intersect(intersect(thirdfrexidx,firstcycidx), firstsnridx);
sum(sum(distAlt(worstidx,:)>1.5)) / prod(size(distAlt(worstidx,:))) * 100


%1st SNR
sum(sum(distAlt(firstsnridx,:)>1.5)) / prod(size(distAlt(firstsnridx,:)))*100
%2nd SNR
sum(sum(distAlt(secondsnridx,:)>1.5)) / prod(size(distAlt(secondsnridx,:)))*100
%3rd SNR
sum(sum(distAlt(thirdsnridx,:)>1.5)) / prod(size(distAlt(thirdsnridx,:)))*100
 
%%
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

    figure; violinplot([resTPtotal(:,1) resTPtotal(:,2) resTPtotal(:,3) resTPtotal(:,4) resTPtotal(:,5)] )

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

% INdividual var plots
figure,
subplot(221),hold on
tmp = resTPtotal(firstvoxidx,:); tmp = mean(tmp); firstvoxerrorsTP = std(tmp(:)) ./ length(tmp); 
tmp = resTPtotal(secondvoxidx,:); tmp = mean(tmp); secondvoxerrorsTP = std(tmp(:))./ length(tmp);  
tmp = resTPtotal(thirdvoxidx,:); tmp = mean(tmp); thirdvoxerrorsTP = std(tmp(:))./ length(tmp);  
errorbar([firstvoxTP, secondvoxTP, thirdvoxTP], [firstvoxerrorsTP, secondvoxerrorsTP, thirdvoxerrorsTP], 'k-.'),
xlim([0.5 3.5]), , title('Voxels (depth)')

subplot(222)
tmp = resTPtotal(firstfrexidx,:); tmp = mean(tmp);firstfrexerrorsTP = std(tmp(:))./ length(tmp); 
tmp = resTPtotal(secondfrexidx,:); tmp = mean(tmp);secondfrexerrorsTP = std(tmp(:))./ length(tmp);  
tmp = resTPtotal(thirdfrexidx,:); tmp = mean(tmp);thirdfrexerrorsTP = std(tmp(:))./ length(tmp); 
errorbar([firstfrexTP, secondfrexTP, thirdfrexTP], [firstfrexerrorsTP, secondfrexerrorsTP, thirdfrexerrorsTP], 'k-.'),
xlim([0.5 3.5]), title('Frequencies')

subplot(223)
tmp = resTPtotal(firstcycidx,:); tmp = mean(tmp);firstcycerrorsTP = std(tmp(:))./ length(tmp); 
tmp = resTPtotal(secondcycidx,:); tmp = mean(tmp);secondcycerrorsTP = std(tmp(:))./ length(tmp);  
tmp = resTPtotal(thirdcycidx,:); tmp = mean(tmp);thirdcycerrorsTP = std(tmp(:))./ length(tmp); 
errorbar([firstcycTP, secondcycTP, thirdcycTP], [firstcycerrorsTP, secondcycerrorsTP, thirdcycerrorsTP], 'k-.'),
xlim([0.5 3.5]), title('Cycles')

subplot(224)
tmp = resTPtotal(firstsnridx,:); tmp = mean(tmp);firstsnrerrorsTP = std(tmp(:))./ length(tmp);  
tmp = resTPtotal(secondsnridx,:); tmp = mean(tmp);secondsnrerrorsTP = std(tmp(:))./ length(tmp);  
tmp = resTPtotal(thirdsnridx,:); tmp = mean(tmp);thirdsnrerrorsTP = std(tmp(:))./ length(tmp);  
errorbar([firstsnrTP, secondsnrTP, thirdsnrTP], [firstsnrerrorsTP, secondsnrerrorsTP, thirdsnrerrorsTP], 'k-.'),
xlim([0.5 3.5]), title('SNR')

% FP
figure,
subplot(221),hold on
tmp = resFPtotal(firstvoxidx,:); tmp = mean(tmp);firstvoxerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(secondvoxidx,:); tmp = mean(tmp);secondvoxerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(thirdvoxidx,:); tmp = mean(tmp);thirdvoxerrorsFP = std(tmp(:))./ length(tmp);  
errorbar([firstvoxFP, secondvoxFP, thirdvoxFP], [firstvoxerrorsFP, secondvoxerrorsFP, thirdvoxerrorsFP], 'r-.'),
xlim([0.5 3.5]), title('Voxels (depth)')

subplot(222)
tmp = resFPtotal(firstfrexidx,:);tmp = mean(tmp); firstfrexerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(secondfrexidx,:); tmp = mean(tmp);secondfrexerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(thirdfrexidx,:); tmp = mean(tmp);thirdfrexerrorsFP = std(tmp(:))./ length(tmp);  
errorbar([firstfrexFP, secondfrexFP, thirdfrexFP], [firstfrexerrorsFP, secondfrexerrorsFP, thirdfrexerrorsFP], 'r-.'),
xlim([0.5 3.5]),  title('Frequencies')

subplot(223)
tmp = resFPtotal(firstcycidx,:);tmp = mean(tmp); firstcycerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(secondcycidx,:);tmp = mean(tmp); secondcycerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(thirdcycidx,:);tmp = mean(tmp); thirdcycerrorsFP = std(tmp(:))./ length(tmp);  
errorbar([firstcycFP, secondcycFP, thirdcycFP], [firstcycerrorsFP, secondcycerrorsFP, thirdcycerrorsFP], 'r-.'),
xlim([0.5 3.5]), title('Cycles')

subplot(224)
tmp = resFPtotal(firstsnridx,:); tmp = mean(tmp);firstsnrerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(secondsnridx,:); tmp = mean(tmp);secondsnrerrorsFP = std(tmp(:))./ length(tmp);  
tmp = resFPtotal(thirdsnridx,:); tmp = mean(tmp);thirdsnrerrorsFP = std(tmp(:))./ length(tmp);  
errorbar([firstsnrFP, secondsnrFP, thirdsnrFP], [firstsnrerrorsFP, secondsnrerrorsFP, thirdsnrerrorsFP], 'r-.'),
xlim([0.5 3.5]),   title('SNR')

%% Best conditions: model 9,36 and 63 (all voxels, freqs 2, cycles 3, SNR 2 and 3)secondfrexidx

% 1 freq
    
    % 3 cycles
    fx1cyc3snr3 = mean(mean(resTPtotal(firstfrexidx(9:9:end),:))); % snr 3
    fx1cyc3snr2 = mean(mean(resTPtotal(firstfrexidx(8:9:end),:))); % snr 2
    fx1cyc3snr1 = mean(mean(resTPtotal(firstfrexidx(7:9:end),:))); % snr 1
    
    % 2 cycles
    fx1cyc2snr3 = mean(mean(resTPtotal(firstfrexidx(6:9:end),:))); % snr 3
    fx1cyc2snr2 = mean(mean(resTPtotal(firstfrexidx(5:9:end),:))); % snr 2
    fx1cyc2snr1 = mean(mean(resTPtotal(firstfrexidx(4:9:end),:))); % snr 1
    
    % 1 cycle
    fx1cyc1snr3 = mean(mean(resTPtotal(firstfrexidx(3:9:end),:))); % snr 3
    fx1cyc1snr2 = mean(mean(resTPtotal(firstfrexidx(2:9:end),:))); % snr 2
    fx1cyc1snr1 = mean(mean(resTPtotal(firstfrexidx(1:9:end),:))); % snr 1

% 2 freq
    
    % 3 cycles
    fx2cyc3snr3 = mean(mean(resTPtotal(secondfrexidx(9:9:end),:))); % snr 3
    fx2cyc3snr2 = mean(mean(resTPtotal(secondfrexidx(8:9:end),:))); % snr 2
    fx2cyc3snr1 = mean(mean(resTPtotal(secondfrexidx(7:9:end),:))); % snr 1
    
    % 2 cycles
    fx2cyc2snr3 = mean(mean(resTPtotal(secondfrexidx(6:9:end),:))); % snr 3
    fx2cyc2snr2 = mean(mean(resTPtotal(secondfrexidx(5:9:end),:))); % snr 2
    fx2cyc2snr1 = mean(mean(resTPtotal(secondfrexidx(4:9:end),:))); % snr 1
    
    % 1 cycle
    fx2cyc1snr3 = mean(mean(resTPtotal(secondfrexidx(3:9:end),:))); % snr 3
    fx2cyc1snr2 = mean(mean(resTPtotal(secondfrexidx(2:9:end),:))); % snr 2
    fx2cyc1snr1 = mean(mean(resTPtotal(secondfrexidx(1:9:end),:))); % snr 1

% 3 freq
    
    % 3 cycles
    fx3cyc3snr3 = mean(mean(resTPtotal(thirdfrexidx(9:9:end),:))); % snr 3
    fx3cyc3snr2 = mean(mean(resTPtotal(thirdfrexidx(8:9:end),:))); % snr 2
    fx3cyc3snr1 = mean(mean(resTPtotal(thirdfrexidx(7:9:end),:))); % snr 1
    
    % 2 cycles
    fx3cyc2snr3 = mean(mean(resTPtotal(thirdfrexidx(6:9:end),:))); % snr 3
    fx3cyc2snr2 = mean(mean(resTPtotal(thirdfrexidx(5:9:end),:))); % snr 2
    fx3cyc2snr1 = mean(mean(resTPtotal(thirdfrexidx(4:9:end),:))); % snr 1
    
    % 1 cycle
    fx3cyc1snr3 = mean(mean(resTPtotal(thirdfrexidx(3:9:end),:))); % snr 3
    fx3cyc1snr2 = mean(mean(resTPtotal(thirdfrexidx(2:9:end),:))); % snr 2
    fx3cyc1snr1 = mean(mean(resTPtotal(thirdfrexidx(1:9:end),:))); % snr 1

xlabls = {'3', '10', '20', '3', '10', '20', '3', '10', '20'};

    figure,
    subplot(311)
    b = bar([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 0 fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1 0  fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1])
    ylim([0 100]), title('SNR 0.5'), set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(2,:) = [0.4667 0.8902 0.4196];
    b.CData(3,:) = [0.0824 0.7412 0.0078]; b.CData(5,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(7,:) = [1.0000 0.9686 0];b.CData(9,:) = [0.9098 0.6000 0.5725];b.CData(10,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on
   % Get standard errors
        init = 1; step = 9;
        for lp = 1:3
            tmp1snr = resTPtotal(firstfrexidx(init:step:end),:); snr1error(lp) = std(tmp1snr(:))./ length(tmp1snr(:));
            tmp2snr = resTPtotal(secondfrexidx(init:step:end),:); snr2error(lp) = std(tmp2snr(:))./ length(tmp2snr(:));
            tmp3snr = resTPtotal(thirdfrexidx(init:step:end),:); snr3error(lp) = std(tmp3snr(:))./ length(tmp3snr(:));
            init=init+3;
        end 
    er = errorbar([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 0 fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1 0  fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1], [snr1error 0 snr2error 0 snr3error]);      
    er.Color = [0 0 0]; er.LineStyle = 'none';  
  
    subplot(312)
    b = bar([fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2 0 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2 0 fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2])
    ylim([0 100]), title('SNR 1'), set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(2,:) = [0.4667 0.8902 0.4196];
    b.CData(3,:) = [0.0824 0.7412 0.0078]; b.CData(5,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(7,:) = [1.0000 0.9686 0];b.CData(9,:) = [0.9098 0.6000 0.5725];b.CData(10,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on
    % Get standard errors
        init = 2; step = 9;
        for lp = 1:3
            tmp1snr = resTPtotal(firstfrexidx(init:step:end),:); snr1error(lp) = std(tmp1snr(:))./ length(tmp1snr(:));
            tmp2snr = resTPtotal(secondfrexidx(init:step:end),:); snr2error(lp) = std(tmp2snr(:))./ length(tmp2snr(:));
            tmp3snr = resTPtotal(thirdfrexidx(init:step:end),:); snr3error(lp) = std(tmp3snr(:))./ length(tmp3snr(:));
            init=init+3;
        end 
    er = errorbar([fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2 0 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2 0 fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2], [snr1error 0 snr2error 0 snr3error]);      
    er.Color = [0 0 0]; er.LineStyle = 'none';  
    
    subplot(313)
    b = bar([fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3 0 fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3 0 fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3]);
    ylim([0 100]), title('SNR 1.5'),set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(2,:) = [0.4667 0.8902 0.4196];
    b.CData(3,:) = [0.0824 0.7412 0.0078]; b.CData(5,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(7,:) = [1.0000 0.9686 0];b.CData(9,:) = [0.9098 0.6000 0.5725];b.CData(10,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on
    % Get standard errors
        init = 3; step = 9;
        for lp = 1:3
            tmp1snr = resTPtotal(firstfrexidx(init:step:end),:); snr1error(lp) = std(tmp1snr(:))./ length(tmp1snr(:));
            tmp2snr = resTPtotal(secondfrexidx(init:step:end),:); snr2error(lp) = std(tmp2snr(:)) ./ length(tmp2snr(:));
            tmp3snr = resTPtotal(thirdfrexidx(init:step:end),:); snr3error(lp) = std(tmp3snr(:)) ./ length(tmp3snr(:));
            init=init+3;
        end 
    er = errorbar([fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3 0 fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3 0 fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3], [snr1error 0 snr2error 0 snr3error]);      
    er.Color = [0 0 0]; er.LineStyle = 'none';  
  



    % Same but another axis
xlabls = {'θ', 'α', 'β', 'θ', 'α', 'β','θ', 'α', 'β'};

    figure,
    subplot(311)
    b = bar([fx1cyc1snr1 fx2cyc1snr1 fx3cyc1snr1 0 fx1cyc2snr1 fx2cyc2snr1 fx3cyc2snr1 0  fx1cyc3snr1 fx2cyc3snr1 fx3cyc3snr1])
    ylim([0 100]), title('SNR 0.5'), set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(5,:) = [0.4667 0.8902 0.4196];
    b.CData(9,:) = [0.0824 0.7412 0.0078]; b.CData(2,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(10,:) = [1.0000 0.9686 0];b.CData(3,:) = [0.9098 0.6000 0.5725];b.CData(7,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on
  
    subplot(312)
    b = bar([fx1cyc1snr2 fx2cyc1snr2 fx3cyc1snr2 0 fx1cyc2snr2 fx2cyc2snr2 fx3cyc2snr2 0  fx1cyc3snr2 fx2cyc3snr2 fx3cyc3snr2])
    ylim([0 100]), title('SNR 1'), set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(5,:) = [0.4667 0.8902 0.4196];
    b.CData(9,:) = [0.0824 0.7412 0.0078]; b.CData(2,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(10,:) = [1.0000 0.9686 0];b.CData(3,:) = [0.9098 0.6000 0.5725];b.CData(7,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on
    
    subplot(313)
    b = bar([fx1cyc1snr3 fx2cyc1snr3 fx3cyc1snr3 0 fx1cyc2snr3 fx2cyc2snr3 fx3cyc2snr3 0  fx1cyc3snr3 fx2cyc3snr3 fx3cyc3snr3])
    ylim([0 100]), title('SNR 1.5'), set(gca,'XTick',[]), xticks([1:3 5:7 9:11]), xticklabels(xlabls)
    b.FaceColor = 'flat'; b.CData(1,:) = [0.8118 0.9294 0.6784]; b.CData(5,:) = [0.4667 0.8902 0.4196];
    b.CData(9,:) = [0.0824 0.7412 0.0078]; b.CData(2,:) = [0.9608 0.9490 0.5961]; b.CData(6,:) = [1.0000 0.9882 0.3882];
    b.CData(10,:) = [1.0000 0.9686 0];b.CData(3,:) = [0.9098 0.6000 0.5725];b.CData(7,:) = [0.8902 0.3490 0.3020];
    b.CData(11,:) = [0.9020 0.1059 0.0353];
    hold on


%% Otra forma de representar
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

load('Z:\OMEGA\Enrique\RdBu_R.mat')
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
print(gcf, ['Z:\OMEGA\Enrique\figures\TPorig.tiff'], '-dtiff', '-r300');
% MEANS
% total
mean(mean(mtxsnr1))
mean(mean(mtxsnr2))
mean(mean(mtxsnr3))
% ciclos
[mean(mtxsnr1,2) mean(mtxsnr2,2) mean(mtxsnr3,2)]
% frex
[mean(mtxsnr1,1)' mean(mtxsnr2,1)' mean(mtxsnr3,1)']




%% FP
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

load('Z:\OMEGA\Enrique\RdBu_R.mat')
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
print(gcf, ['Z:\OMEGA\Enrique\figures\FPorig.tiff'], '-dtiff', '-r300');

max(max(mtxsnr1))
max(max(mtxsnr2))
max(max(mtxsnr3))
%%

% Means
snr1perf = mean([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1  fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1])
snr2perf = mean([fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2  fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2])
snr3perf = mean([fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3 fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3  fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3])
snr1sd = std([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1  fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1])
snr2sd = std([fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2  fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2])
snr3sd = std([fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3 fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3  fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3])

frex1perf = mean([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2  fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3])
frex2perf = mean([fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2  fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3])
frex3perf = mean([fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1 fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2  fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3])
frex1sd = std([fx1cyc1snr1 fx1cyc2snr1 fx1cyc3snr1 fx1cyc1snr2 fx1cyc2snr2 fx1cyc3snr2  fx1cyc1snr3 fx1cyc2snr3 fx1cyc3snr3])
frex2sd = std([fx2cyc1snr1 fx2cyc2snr1 fx2cyc3snr1 fx2cyc1snr2 fx2cyc2snr2 fx2cyc3snr2  fx2cyc1snr3 fx2cyc2snr3 fx2cyc3snr3])
frex3sd = std([fx3cyc1snr1 fx3cyc2snr1 fx3cyc3snr1 fx3cyc1snr2 fx3cyc2snr2 fx3cyc3snr2  fx3cyc1snr3 fx3cyc2snr3 fx3cyc3snr3])

cyc1perf = mean([fx1cyc1snr1 fx2cyc1snr1 fx3cyc1snr1 fx1cyc1snr2 fx2cyc1snr2 fx3cyc1snr2  fx1cyc1snr3 fx2cyc1snr3 fx3cyc1snr3])
cyc2perf = mean([fx1cyc2snr1 fx2cyc2snr1 fx3cyc2snr1 fx1cyc2snr2 fx2cyc2snr2 fx3cyc2snr2  fx1cyc2snr3 fx2cyc2snr3 fx3cyc2snr3])
cyc3perf = mean([fx1cyc3snr1 fx2cyc3snr1 fx3cyc3snr1 fx1cyc3snr2 fx2cyc3snr2 fx3cyc3snr2  fx1cyc3snr3 fx2cyc3snr3 fx3cyc3snr3])
cyc1sd = std([fx1cyc1snr1 fx2cyc1snr1 fx3cyc1snr1 fx1cyc1snr2 fx2cyc1snr2 fx3cyc1snr2  fx1cyc1snr3 fx2cyc1snr3 fx3cyc1snr3])
cyc2sd = std([fx1cyc2snr1 fx2cyc2snr1 fx3cyc2snr1 fx1cyc2snr2 fx2cyc2snr2 fx3cyc2snr2  fx1cyc2snr3 fx2cyc2snr3 fx3cyc2snr3])
cyc3sd = std([fx1cyc3snr1 fx2cyc3snr1 fx3cyc3snr1 fx1cyc3snr2 fx2cyc3snr2 fx3cyc3snr2  fx1cyc3snr3 fx2cyc3snr3 fx3cyc3snr3])

