function effects = transmissionBOmega(from, B, Omega, varAnd, varNot, multiplier)
    effects = cell(1, length(varAnd));
    for ii = 1:length(varAnd)
        vAnd = varAnd{ii};
        vNot = varNot{ii};
        m = multiplier(ii);

        BTilde = B;
        OmegaTilde = Omega;
        
        for v = vAnd
            [BTilde, OmegaTilde] = applyAndToB(BTilde, OmegaTilde, from, v);
        end
        
        for v = vNot
            [BTilde, OmegaTilde] = applyNotToB(BTilde, OmegaTilde, from, v);
        end
        
        effects{ii} = (eye(size(BTilde)) - BTilde) \ OmegaTilde(:, from);
        effects{ii} = m * effects{ii};
        
        % if isempty(vAnd) && isempty(vNot)
        %     return;
        % end
        
        if ~isempty(vAnd)
            effects{ii}(1:max(vAnd)) = 0;
        end
    end

    effects = sum(cat(2, effects{:}), 2);
    if ~isempty(varAnd)
        effects(1:max(cat(2, varAnd{:}))) = 0;
    end

end

