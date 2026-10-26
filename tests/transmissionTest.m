function tests = utilsTest
    tests = functiontests(localfunctions);
end

function testTransmission1(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    cond = makeCondition("!x2");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    % Doing it all manually
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(2,1) = 0;
    BTilde(2, :) = 0;

    manualBOmega = (eye(size(B)) - BTilde) \ OmegaTilde;
    manualBOmega = manualBOmega(:, 1);

    manualIrfs = irfs(:, 1) - irfs(2,1)*irfsOrtho(:,2)/irfsOrtho(2,2);

    assert(all(max(abs(effect - manualBOmega), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect - effectIrfs), [], 'all') < sqrt(eps())));
end

function testTransmission2(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    manualIrfs = irfs(3,1) .* irfsOrtho(:,3) ./ irfsOrtho(3, 3) -  ...
        irfs(2,1) .* irfsOrtho(3,2) ./ irfsOrtho(2,2) .* irfsOrtho(:,3) ./ irfsOrtho(3,3);

    BTilde = B;
    OmegaTilde = Omega; 
    OmegaTilde(4,1) = 0;
    BTilde(4, [1, 2, 4]) = 0;
    manualBOmegaPart1 = (eye(size(B)) - BTilde) \ OmegaTilde; 

    BTilde = B; 
    OmegaTilde = Omega; 
    OmegaTilde([3, 4], 1) = 0; 
    BTilde(3, 1) = 0;
    BTilde(4, [1, 2, 4]) = 0;
    manualBOmegaPart2 = (eye(size(B)) - BTilde) \ OmegaTilde; 
    manualBOmega = manualBOmegaPart1 - manualBOmegaPart2;
    manualBOmega = manualBOmega(:, 1);

    cond = makeCondition("x3 & !x2");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    % all paths must go through x3, so x1 to x3 are masked for both methods
    assert(all(isnan(effect(1:3))));
    assert(all(isnan(effectIrfs(1:3))));
    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(4:end) - effectIrfs(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs(4:end) - manualIrfs(4:end)), [], 'all') < sqrt(eps())));
end

function testTransmission3(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    cond = makeCondition("x2 & !x2");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    assert(all(effect == 0));
    assert(all(effectIrfs == 0));
end

function testTransmission4(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    manualIrfs = irfs(2, 1) .* irfsOrtho(:, 2) ./ irfsOrtho(2, 2) + ...
        irfs(3, 1) .* irfsOrtho(:, 3) ./ irfsOrtho(3, 3) -  ...
        2 * irfs(2, 1) .* irfsOrtho(3, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(:, 3) ./ irfsOrtho(3, 3);
    
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde([3, 4], 1) = 0;
    BTilde(3:end, 1) = 0;
    manualBOmegaPart1 = (eye(size(B)) - BTilde) \ OmegaTilde;

    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(4,1) = 0;
    BTilde(4:end, [1, 2, 4]) = 0;
    manualBOmegaPart2 = (eye(size(B)) - BTilde) \ OmegaTilde;


    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde([3, 4], 1) = 0;
    BTilde([3, 4], 1) = 0;
    BTilde(4:end, [1, 2]) = 0;
    manualBOmegaPart3 = (eye(size(B)) - BTilde) \ OmegaTilde;

    manualBOmega = manualBOmegaPart1 + manualBOmegaPart2 - 2 * manualBOmegaPart3;
    manualBOmega = manualBOmega(:, 1);

    
    cond = makeCondition("((x2 & !x3) | (!x2 & x3))");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");
   
    % the highest AND variable is x3, so x1 to x3 are masked for both methods
    assert(all(isnan(effect(1:3))));
    assert(all(isnan(effectIrfs(1:3))));
    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs(4:end) - manualIrfs(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(4:end) - effectIrfs(4:end)), [], 'all') < sqrt(eps())));
end


function testTransmission5(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';
    
    manualIrfs = irfs(2, 1) .* irfsOrtho(:, 2) ./ irfsOrtho(2, 2) + ...
        irfs(5, 1) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5) - ...
        irfs(2, 1) .* irfsOrtho(5, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5);
    
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(3:end, 1) = 0;
    BTilde(3:end, 1) = 0;
    manualBOmegaPart1 = (eye(size(B)) - BTilde) \ OmegaTilde;
    
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(6:end, 1) = 0;
    BTilde(6:end, 1:4) = 0;
    manualBOmegaPart2 = (eye(size(B)) - BTilde) \ OmegaTilde; 
    
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(3:end, 1) = 0;
    BTilde(3:end, 1) = 0;
    BTilde(6:end, 1:4) = 0;
    manualBOmegaPart3 = (eye(size(B)) - BTilde) \ OmegaTilde; 
    
    manualBOmega = manualBOmegaPart1 + manualBOmegaPart2 - manualBOmegaPart3;
    manualBOmega = manualBOmega(:, 1); 
    
    cond = makeCondition("x2 | x5");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf"); 
    
    % the highest AND variable is x5, so x1 to x5 are masked for both methods
    assert(all(isnan(effect(1:5))));
    assert(all(isnan(effectIrfs(1:5))));
    assert(all(max(abs(effect(6:end) - manualBOmega(6:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs(6:end) - manualIrfs(6:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(6:end) - effectIrfs(6:end)), [], 'all') < sqrt(eps())));
end


function testTransmission6(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';
    
    manualIrfs = irfs(2, 1) .* irfsOrtho(:, 2) ./ irfsOrtho(2, 2) - ...
        irfs(2, 1) .* irfsOrtho(3, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(:, 3) ./ irfsOrtho(3, 3) -  ...
        irfs(2, 1) .* irfsOrtho(4, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(:, 4) ./ irfsOrtho(4, 4) +  ...
        irfs(2, 1) .* irfsOrtho(3, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(4, 3) ./ irfsOrtho(3, 3) .* irfsOrtho(:, 4) ./ irfsOrtho(4, 4) -  ...
        irfs(2, 1) .* irfsOrtho(5, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5) +  ...
        irfs(2, 1) .* irfsOrtho(3, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(5, 3) ./ irfsOrtho(3, 3) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5) +  ...
        irfs(2, 1) .* irfsOrtho(4, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(5, 4) ./ irfsOrtho(4, 4) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5) -  ...
        irfs(2, 1) .* irfsOrtho(3, 2) ./ irfsOrtho(2, 2) .* irfsOrtho(4, 3) ./ irfsOrtho(3, 3) .* irfsOrtho(5, 4) ./ irfsOrtho(4, 4) .* irfsOrtho(:, 5) ./ irfsOrtho(5, 5);
    
    % too many terms to calculate it using the second method
    
    cond = makeCondition("x2 & !x3 & !x4 & !x5");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    % the only AND variable is x2, so x1 and x2 are masked for both methods
    assert(all(isnan(effect(1:2))));
    assert(all(isnan(effectIrfs(1:2))));
    assert(all(max(abs(effectIrfs(3:end) - manualIrfs(3:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(3:end) - effectIrfs(3:end)), [], 'all') < sqrt(eps())));
end


function testTransmission7(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    manualIrfs = irfs(:, 1);
    manualBOmega = (eye(size(B)) - B) \ Omega; 
    manualBOmega = manualBOmega(:, 1);

    cond = makeCondition("(x1 | x2 | x3) | !(x1 | x2 | x3)");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf"); 

    % the condition simplifies to the total effect, so nothing is masked
    assert(~any(isnan(effect)));
    assert(~any(isnan(effectIrfs)));
    assert(all(max(abs(effect - manualBOmega), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect - effectIrfs), [], 'all') < sqrt(eps())));
end

function testTransmission8(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    % algorithms were created in the way that the variables in the condition
    % must be before the outcome variable. To prevent misinterpretation, the
    % effects on all variables up to and including the highest AND variable
    % are set to NaN.

    cond = makeCondition("x1 & x2");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    % Doing is manually
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(2:end, :) = 0;
    BTilde(3:end, 1) = 0;

    manualBOmega = (eye(size(B)) - BTilde) \ OmegaTilde;
    manualBOmega = manualBOmega(:, 1);

    manualIrfs = irfs(1, 1) * irfsOrtho(2,1) / irfsOrtho(1,1) * irfsOrtho(:, 2) / irfsOrtho(2, 2);

    assert(all(isnan(effect(1:2))));
    assert(all(isnan(effectIrfs(1:2))));

    assert(all(max(abs(effect(3:end) - manualBOmega(3:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs(3:end) - manualIrfs(3:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(3:end) - effectIrfs(3:end)), [], 'all') < sqrt(eps())));

end

function testTransmission9(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    % algorithms were created in the way that the variables in the condition
    % must be before the outcome variable. To prevent misinterpretation, the
    % effects on all variables up to and including the highest AND variable
    % of any term are set to NaN.

    cond = makeCondition("(x1 & x2) | x3");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    % Doing it manually. Q((x1 & x2) | x3) = Q(x1 & x2) + Q(x3) - Q(x1 & x2 & x3)
    % Term 1: Q(x1 & x2)
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(2:end, :) = 0;
    BTilde(3:end, 1) = 0;
    manualBOmega1 = (eye(size(B)) - BTilde) \ OmegaTilde;
    manualBOmega1 = manualBOmega1(:, 1);

    manualIrfs1 = irfs(1, 1) * irfsOrtho(2, 1) / irfsOrtho(1, 1) * irfsOrtho(:, 2) / irfsOrtho(2, 2);

    % Term 2: Q(x3)
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(4:end, :) = 0;
    BTilde(4:end, 1:2) = 0;
    manualBOmega2 = (eye(size(B)) - BTilde) \ OmegaTilde;
    manualBOmega2 = manualBOmega2(:, 1);

    manualIrfs2 = irfs(3, 1) * irfsOrtho(:, 3) / irfsOrtho(3, 3);

    % Term 3: Q(x1 & x2 & x3)
    BTilde = B;
    OmegaTilde = Omega;
    OmegaTilde(2:end, :) = 0;
    BTilde(3:end, 1) = 0;
    BTilde(4:end, 2) = 0;
    manualBOmega3 = (eye(size(B)) - BTilde) \ OmegaTilde;
    manualBOmega3 = manualBOmega3(:, 1);

    manualIrfs3 = irfs(1, 1) * irfsOrtho(2, 1) / irfsOrtho(1, 1) * ...
        irfsOrtho(3, 2) / irfsOrtho(2, 2) * irfsOrtho(:, 3) / irfsOrtho(3, 3);

    manualBOmega = manualBOmega1 + manualBOmega2 - manualBOmega3;
    manualIrfs = manualIrfs1 + manualIrfs2 - manualIrfs3;

    assert(all(isnan(effect(1:3))));
    assert(all(isnan(effectIrfs(1:3))));

    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs(4:end) - manualIrfs(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(4:end) - effectIrfs(4:end)), [], 'all') < sqrt(eps())));
end

function testTransmission10(testCase)

    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    % NOT conditions never lead to masking. Paths to x1 cannot go through x2,
    % so the effect on x1 is the total effect, and the effect on x2 is zero
    % since all paths to x2 go through x2.

    cond = makeCondition("!x2");
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");

    assert(~any(isnan(effect)));
    assert(~any(isnan(effectIrfs)));
    assert(abs(effect(1) - irfs(1, 1)) < sqrt(eps()));
    assert(abs(effectIrfs(1) - irfs(1, 1)) < sqrt(eps()));
    assert(abs(effect(2)) < sqrt(eps()));
    assert(abs(effectIrfs(2)) < sqrt(eps()));
end

%-------------------------------------------------------------------------------
% Masking of uninterpretable effects.
% All checks use a simulated system with four variables and are run for both
% the BOmega and the irf method. Effects that are not masked are compared to
% the effect along all paths satisfying the condition, which is computed by
% explicitly enumerating the paths.
%-------------------------------------------------------------------------------

function testMaskingNotX1(testCase)
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("!x1");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(abs(e(1)) < sqrt(eps()));
        assert(abs(e(2) - pathSpecificEffect(B, Omega, 1, 2, [], [1])) < sqrt(eps()));
    end
end

function testMaskingX1(testCase)
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("x1");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(isnan(e(1)));
        assert(abs(e(2) - pathSpecificEffect(B, Omega, 1, 2, [1], [])) < sqrt(eps()));
    end
end

function testMaskingX2(testCase)
    [effect, effectIrfs] = maskingEffects("x2");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(isnan(e(1)));
    end
end

function testMaskingNotX2(testCase)
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("!x2");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(abs(e(1) - pathSpecificEffect(B, Omega, 1, 1, [], [2])) < sqrt(eps()));
    end
end

function testMaskingX2AndNotX3(testCase)
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("x2 & !x3");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(isnan(e(1)));
        assert(isnan(e(2)));
        assert(abs(e(3)) < sqrt(eps()));
        assert(abs(e(4) - pathSpecificEffect(B, Omega, 1, 4, [2], [3])) < sqrt(eps()));
    end
end

function testMaskingNotX2AndX3(testCase)
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("!x2 & x3");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(all(isnan(e(1:3))));
        assert(abs(e(4) - pathSpecificEffect(B, Omega, 1, 4, [3], [2])) < sqrt(eps()));
    end
end

function testMaskingNotX4(testCase)
    % A NOT on the last variable never masks earlier variables (this is the
    % situation of channels such as notThrough('x', 0:maxHorizon, order)).
    [~, ~, irfs, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("!x4");
    for e = {effect, effectIrfs}
        e = e{1};
        assert(max(abs(e(1:3) - irfs(1:3, 1))) < sqrt(eps()));
        assert(abs(e(4)) < sqrt(eps()));
    end
end

function testMaskingDisjunction(testCase)
    % Masking is determined by the highest AND variable of any term of the
    % (simplified) condition. x1 | x3 = x1 + x3 - x1 & x3, so x1 to x3 are
    % masked.
    [B, Omega, ~, ~] = simulateMaskingSystem(4);
    [effect, effectIrfs] = maskingEffects("x1 | x3");
    manual = pathSpecificEffect(B, Omega, 1, 4, [1], []) + ...
        pathSpecificEffect(B, Omega, 1, 4, [3], []) - ...
        pathSpecificEffect(B, Omega, 1, 4, [1, 3], []);
    for e = {effect, effectIrfs}
        e = e{1};
        assert(all(isnan(e(1:3))));
        assert(abs(e(4) - manual) < sqrt(eps()));
    end
end

function testMaskingMethodsAgree(testCase)
    % Both methods must mask the same entries and agree on all others.
    conds = ["!x1", "x1", "x2", "!x2", "x2 & !x3", "!x2 & x3", "!x4", "x1 | x3"];
    for c = conds
        [effect, effectIrfs] = maskingEffects(c);
        assert(isequal(isnan(effect), isnan(effectIrfs)));
        keep = ~isnan(effect);
        assert(max(abs(effect(keep) - effectIrfs(keep)), [], 'all') < sqrt(eps()));
    end
end

function [B, Omega, irfs, irfsOrtho] = simulateMaskingSystem(k)
    % Simulate a structural system with `k` variables and a single horizon.
    % Starting from a Cholesky factor `L` (lower triangular, positive diagonal)
    % and a rotation `Q`, the structural IRFs are `L * Q'` and the
    % orthogonalised IRFs are `L`. The systems form is then
    % `B = I - D * inv(L)` (strictly lower triangular) and `Omega = D * Q'`
    % with `D = diag(diag(L))`, so that `(I - B) \ Omega == L * Q'`.
    s = rng(6150533);
    [Q, ~] = qr(randn(k));
    L = tril(randn(k), -1) + diag(1 + rand(k, 1));
    rng(s);
    D = diag(diag(L));
    B = eye(k) - D / L;
    Omega = D * Q';
    irfs = L * Q';
    irfsOrtho = L;
end

function [effect, effectIrfs] = maskingEffects(s)
    % Effects of condition `s` on the simulated system using both methods.
    [B, Omega, irfs, irfsOrtho] = simulateMaskingSystem(4);
    cond = makeCondition(s);
    effect = transmission(1, B, Omega, cond, "BOmega");
    effectIrfs = transmission(1, irfs, irfsOrtho, cond, "irf");
end

function effect = pathSpecificEffect(B, Omega, from, outcome, varAnd, varNot)
    % Effect of shock `from` on `outcome` along all paths that go through all
    % variables in `varAnd` and through none of the variables in `varNot`.
    % Paths are enumerated explicitly, making this independent of both methods.
    effect = 0;
    for mask = 0:(2^(outcome-1) - 1)
        nodes = [find(bitget(mask, 1:(outcome-1))), outcome];
        if ~all(ismember(varAnd, nodes)) || any(ismember(varNot, nodes))
            continue;
        end
        e = Omega(nodes(1), from);
        for t = 1:(length(nodes) - 1)
            e = e * B(nodes(t+1), nodes(t));
        end
        effect = effect + e;
    end
end
