% ---- 2. Helper function: generate random compositions with given bounds --
function frac = constrained_composition(nSamples, nComp, low, high)
    % Generate nSamples compositions with each component in [low, high] and sum = 1.
    % Uses rejection sampling (efficient for nComp=4 and feasible range).
    frac = zeros(nSamples, nComp);
    for i = 1:nSamples
        ok = false;
        while ~ok
            x = low + (high-low)*rand(1,nComp);
            if abs(sum(x)-1) < 1e-6
                frac(i,:) = x;
                ok = true;
            end
        end
    end
end