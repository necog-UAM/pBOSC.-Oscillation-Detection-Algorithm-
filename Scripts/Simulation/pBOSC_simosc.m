function oscillation = pBOSC_simosc(cfg)
 
%% Generate an oscilltion with frequency and cycles centered 
dur = 0:1/cfg.fsample:(1/cfg.freq*cfg.cycles);
sim = sin(2*pi*cfg.freq*dur); 

zpad = zeros(1,(length(cfg.time) - length(sim)));

if mod(length(zpad),2)==0
    oscillation = [zpad(1:length(zpad)/2) sim zpad(1:length(zpad)/2)];
else
    oscillation = [zpad(1:ceil(length(zpad)/2)) sim zpad(1:floor(length(zpad)/2))];
end
  
end