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

    % only checking elements after 3, because all must go through x3
    % no problem for irf method
    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(4:end) - effectIrfs(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
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
   
    % all paths must go through x3 so we can only compare from 4 onwards 
    % for the BOmega method
    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
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
    
    % all paths must go through x5 so we can only compare from 6 onwards 
    % for the BOmega method
    assert(all(max(abs(effect(6:end) - manualBOmega(6:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
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

    % all paths must go through x5 so we can only compare from 6 onwards 
    % for the BOmega method
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(5:end) - effectIrfs(5:end)), [], 'all') < sqrt(eps())));
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

    % all paths must go through x3 so we can only compare from 4 onwards 
    % for the BOmega method
    assert(all(max(abs(effect(4:end) - manualBOmega(4:end)), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect(4:end) - effectIrfs(4:end)), [], 'all') < sqrt(eps())));
end

function testTransmission8(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    % If the condition is (x1 AND x2), then the effect on x1 should be zero 
    % because the path does not include x2 yet, but the effect on x2 can 
    % be anything because x2 is on the path. 

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
    manualBOmega(1, 1) = 0;   % Because x2 is not on the path and user might interpret it.

    manualIrfs = irfs(1, 1) * irfsOrtho(2,1) / irfsOrtho(1,1) * irfsOrtho(:, 2) / irfsOrtho(2, 2);

    assert(all(max(abs(effect - manualBOmega), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect - effectIrfs), [], 'all') < sqrt(eps())));

end

function testTransmission9(testCase)
    B = jsondecode(fileread('./tests/simulated-svar-k3-p1/B.json'))';
    Omega = jsondecode(fileread('./tests/simulated-svar-k3-p1/Omega.json'))';
    irfs = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs.json'))';
    irfsOrtho = jsondecode(fileread('./tests/simulated-svar-k3-p1/irfs_ortho.json'))';

    % If the condition is (x1 AND x2) OR x3, then it makes not much sense 
    % to interpret the effect on x1 and x2 even though the effect on x2 could 
    % be nonzero. Thus, to prevent the user from making mistakes, we should 
    % simply set these to zero. However, what does make sense is to interpret 
    % the effect on x3, since x3 is on the path if it is the outcome variable. 
    % Thus, the user might want to provide the outcome as the last variable that 
    % must lie on the path which means we shouldn't just set this value to zero 
    % to protect the user. 

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
    manualBOmega1(1, 1) = 0;  % To prevent misinterpretation

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
    manualBOmega3(1:2, 1) = 0;  % To prevent misinterpretation

    manualIrfs3 = irfs(1, 1) * irfsOrtho(2, 1) / irfsOrtho(1, 1) * ...
        irfsOrtho(3, 2) / irfsOrtho(2, 2) * irfsOrtho(:, 3) / irfsOrtho(3, 3);

    manualBOmega = manualBOmega1 + manualBOmega2 - manualBOmega3;
    manualIrfs = manualIrfs1 + manualIrfs2 - manualIrfs3;

    assert(all(max(abs(effect - manualBOmega), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effectIrfs - manualIrfs), [], 'all') < sqrt(eps())));
    assert(all(max(abs(effect - effectIrfs), [], 'all') < sqrt(eps())));
end
