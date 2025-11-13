% Aperiodic generators. Randomly selected voxels are at a distance (5mm) of
% oscillation sources.

function apgens = Sim_apgenerators(cfg)

% Load source models
[source] = Omega_source(3423, 5); % 2nd paramter is the radius distance (mm)

rng(cfg.seed)

Nvox = length(source.inside);

for it=1:cfg.repetitions
    gen = randperm(Nvox,cfg.numgens);
    while ~isempty(intersect(gen,find(sum(source.connmat(cfg.Nsources,:))))) % If near neighbor of any sine generator, repeat
            gen = randperm(Nvox,cfg.numgens);
    end
    apgens(it,:) = gen;
end

end