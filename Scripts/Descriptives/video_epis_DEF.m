addpath('Z:\Toolbox\fieldtrip-20230118\')
addpath(genpath('Z:\OMEGA\OMEGA-NaturalFrequencies-main'))
cd('Z:\OMEGA\OMEGA_data\sub-0001\ses-0001')
load epis_conctd
load datasource
load source_forward_10mm 
load source_inverse_10mm 

load ('Z:\OMEGA\scripts\splot_spdup_left')
connmat_left = double(connmat);
surf_left = surf;
maskval_left = maskval;

load ('Z:\OMEGA\scripts\splot_spdup_right')
connmat_right = double(connmat);
surf_right = surf;
maskval_right = maskval;

voxel_inside = find(source.inside==1);

J = jet_omega_mod;
f0 = size(J,1);
foi = exp(0.6:0.05:3.75); % de 1.8 Hz a 40 Hz
Nf = length(foi);
Nvox = length(epis);

cfg = [];
cfg.resamplefs = 128;
dataclean = ft_resampledata(cfg,dataclean);

data = dataclean.trial{1};
time = dataclean.time{1};
sr =  dataclean.fsample;

oscactiv = zeros(Nvox,length(time));
for vx = 1:Nvox
    for ep = 1:length(epis{vx})
        if epis{vx}(ep).dur_cyc>=3
            oscactiv(vx,epis{vx}(ep).timeps) = round(epis{vx}(ep).freq);
        end
    end
end


%% Prepare figures for video

cd ('Z:\OMEGA\video_epis')
Ncol = 12;
tt = sr*1;    % number of seconds per figure
ctv=1;
for t0 = 1:3:3000%length(time)
    tic
    t0
    figure('WindowState','maximized','Color',[1 1 1]);
    vx=0;
    tf = t0+tt;
    for q = 1:Ncol
        if q==7   
            title(['              Time: ' num2str(round(time(t0)+(time(tf)-time(t0))/2,2)) ' s (' num2str((t0+tf)/2) ')'], 'FontSize', 18)
        end
        subplot(4,Ncol,[12+q:Ncol:4*Ncol])
        vx0 = floor(Nvox/12)*(q-1)+2;
        vxf = floor(Nvox/12)*q-2;
        ytck = [];
        ytcklab = {};
        ct = 1;                 
        for vx = vx0:4:vxf          % plots for 1/4 voxels
            sc =  .7*vx*10e-13;     % scale
            plot(time(t0:tf),data(vx,t0:tf) - sc,'Color',[0.7 0.7 0.7]);
            set(gca,'Box','off','GridAlpha',1,'GridColor',[0 0 0],'XColor','none','YColor',[0.7 0.7 0.7],'XLim',[time(t0) time(tf)],...
                'XTick',[time(t0)+(time(tf)-time(t0))/2],'XTickLabel',{});

            ytck(ct) = mean(data(vx,t0:tf) - sc);
            ytcklab{ct}=num2str(vx);
            ct=ct+1;

            if vx == vx0
                yf = max(data(vx,t0:tf) - sc);
            elseif vx == vxf
                y0 = min(data(vx,t0:tf) - sc);
            end
            hold on

            osctp = find(oscactiv(vx,t0:tf)>0);
            tm = t0:tf;
            df = diff(osctp);
            oscactiv2 = oscactiv(vx,t0:tf);
            dfosc = diff(oscactiv2(osctp));
            df(abs(dfosc)>0) = 2;
            if ~isempty(df)
                t1 = 1;
                tpf = find(df>1);
                if isempty(tpf)
                    t2=length(df);
                    fcol = findbin(foi,mean(oscactiv(vx,tm(osctp(t1:t2)))));
                    plot(time(tm(osctp(t1:t2))),data(vx,tm(osctp(t1:t2))) - sc,'Color',J(fcol,:),'LineWidth',1.2);
                    set(gca,'XGrid','on','GridLineWidth',1.2)
                else
                    for i=1:length(tpf)+1
                        if i==length(tpf)+1
                            t2=length(df)+1;
                        else
                            t2 = tpf(i);
                        end
                        fcol = findbin(foi,mean(oscactiv(vx,tm(osctp(t1:t2)))));
                        plot(time(tm(osctp(t1:t2))),data(vx,tm(osctp(t1:t2))) - sc,'Color',J(fcol,:),'LineWidth',1.2);
                        set(gca,'XGrid','on','GridLineWidth',1.2)
                        t1 = t2+1;
                    end
                end
            end
        end
        set(gca, 'YTick',flip(ytck),'YTickLabel',flip(ytcklab))
        ylim([y0 yf])
    end
        
    % To speed up ft_sourceplot (see splot_spdup): 

    sloc = oscactiv(:,(t0+tf)/2);

    cfg=[];
    cfg.edgecolor = 'none';
    cfg.facecolor = [];
    cfg.vertexcolor = 'curv';
    cfg.funcolorlim   = [0.7 3.4];      % exponential scale   2-30 Hz
    cfg.opaclim = [0 1];
    cfg.opacitymap = 0.8*ones(1,64);
    cfg.funcolormap = J;

    val = sloc'*connmat_left;
    val(val==0)=NaN;           % Modified by Almu 10/07/2024 to avoid plotting results in subcortical regions
    val = log(val);
    
    subplot(4,Ncol,[5 6])
    ft_plot_mesh(surf_left, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
    ft_plot_mesh(surf_left, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val', 'facealpha', maskval_left, 'clim', [cfg.funcolorlim(1) cfg.funcolorlim(2)], 'alphalim', [cfg.opaclim(1) cfg.opaclim(2)], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
    view([90 0]),  camlight('left') 

    subplot(4,Ncol,[9 10])
    ft_plot_mesh(surf_left, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
    ft_plot_mesh(surf_left, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val', 'facealpha', maskval_left, 'clim', [cfg.funcolorlim(1) cfg.funcolorlim(2)], 'alphalim', [cfg.opaclim(1) cfg.opaclim(2)], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
    view([-90 0]), camlight('left')

    val = sloc'*connmat_right;
    val(val==0)=NaN;           % Modified by Almu 10/07/2024 to avoid plotting results in subcortical regions
    val = log(val);

    subplot(4,Ncol,[3 4])
    ft_plot_mesh(surf_right, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
    ft_plot_mesh(surf_right, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val', 'facealpha', maskval_right, 'clim', [cfg.funcolorlim(1) cfg.funcolorlim(2)], 'alphalim', [cfg.opaclim(1) cfg.opaclim(2)], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
    view([90 0]),  camlight('right')
    
    subplot(4,Ncol,[7 8])
    ft_plot_mesh(surf_right, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
    ft_plot_mesh(surf_right, 'edgecolor', cfg.edgecolor, 'facecolor', cfg.facecolor, 'vertexcolor', val', 'facealpha', maskval_right, 'clim', [cfg.funcolorlim(1) cfg.funcolorlim(2)], 'alphalim', [cfg.opaclim(1) cfg.opaclim(2)], 'alphamap', cfg.opacitymap, 'colormap', cfg.funcolormap, 'maskstyle', 'opacity');
    view([-90 0]), camlight('right')

    % colorbar
    subplot(4,Ncol,12)
    cbar = rep4+
    +mat(linspace(0, 1, 256)', 1, 20);  % Un vector de 256 valores repetidos en 20 columnas
    imagesc(cbar);
    colormap('jet_omega_mod');
    set(gca, 'Position',[0.89 0.785 0.014 0.12])    % [left, bottom, width, height]
    set(gca, 'XTick', []);
    set(gca, 'YTick', [35:50:240]);
    set(gca, 'YTickLabel', {'Delta','Theta','Alpha','Low-Beta','High-Beta'},'FontSize',13);

    print('-djpeg','-r100',['Sub1-Ses1-Frame' num2str(ctv) '.jpeg']);
    ctv = ctv+1;

    close
    toc
end


%% Make video

video = VideoWriter('OscEpisodes.avi'); % Nombre del archivo y formato
video.FrameRate = 2; % Configurar la tasa (frames per second)
open(video); % Abrir el archivo de video

cd ('Z:\OMEGA\video_epis')

for i = 1:878
    i
    img = imread(['Sub1-Ses1-Frame' num2str(i) '.jpeg ']);
    writeVideo(video, img);
end
close(video);



