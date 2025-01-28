function tests = simplifyingTest
    tests = functiontests(localfunctions);
end

function testConstructorQ(testCase)
    q1 = Q('x1');
    q2 = Q('x2', -1.5);
    q3 = Q({'x1', 'x2'});
    q4 = Q({'x1', 'x2'}, [-1, 2]);
    q5 = Q([1, 2, 3]);

    assert(isequal(q1.vars{1}, 'x1'));
    assert(isequal(q2.vars{1}, 'x2'));
    assert(q2.multiplier(1) == -1.5);
    assert(length(q3.vars) == 2);
    assert(all(q4.multiplier == [-1, 2]));
    assert(isequal(q5.vars,{'x1', 'x2', 'x3'}));
end

function testCollectTerms(testCase)
    q = Q({'x1', 'x1'}, [1, 1]);
    q = collectTerms(q);

    assert(length(q.vars) == 1);
    assert(isequal(q.vars{1}, 'x1'));
    assert(q.multiplier(1) == 2);

    % In the matlab version, T denotes true
    q = Q({'x1', 'T', 'x1'}, [1, 1, -1]);
    q = collectTerms(q);

    assert(length(q.vars) == 1);
    assert(isequal(q.vars{1}, 'T'));
    assert(q.multiplier(1) == 1);
end

function testStringAnd(testCase)
    s = stringAnd('T', 'x1');
    assert(isequal(s, 'x1'));

    s = stringAnd('x1', 'x1');
    assert(isequal(s, 'x1'));

    s = stringAnd('x1 & x2', 'x1');
    assert(isequal(s, 'x2 & x1'));

    s = stringAnd('!x1 & x2', 'x1');
    assert(isequal(s, 'x2 & x1 & !x1'));

    s = stringAnd('x1', 'x2');
    assert(isequal(s, 'x2 & x1'));
end

function testCheckContradiction(testCase)
    [hasContradictions, contradictions] = checkContradiction([1], [2]);
    assert(~hasContradictions);

    [hasContradictions, contradictions] = checkContradiction([4, 3, 1], [5, 2]);
    assert(~hasContradictions);

    [hasContradictions, contradictions] = checkContradiction([1], [2, 1]);
    assert(hasContradictions);

    [hasContradictions, contradictions] = checkContradiction([3, 1], [2, 1])
    assert(hasContradictions);
    assert(contradictions(1) == 0);
    assert(contradictions(2) == 1);

    [hasContradictions, contradictions] = checkContradiction({[3, 1], [2, 1]}, {[5, 4], [5, 4]});
    assert(~hasContradictions);

    [hasContradictions, contradictions] = checkContradiction({[3, 1], [2, 1]}, {[5, 1], [2, 4]});
    assert(hasContradictions);
    assert(contradictions(1) == 1);
    assert(contradictions(2) == 1);

    [hasContradictions, contradictions] = checkContradiction({[3, 1], [2, 1]}, {[5, 1], [5, 4]});
    assert(hasContradictions);
    assert(contradictions(1) == 1);
    assert(contradictions(2) == 0);
end

function testGetVarNumsAndMultiplier(testCase)
    q = Q({'x1'}, 1);
    [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q);
    assert(all(andNums{1} == 1));
    assert(length(andNotNums{1}) == 0);

    q = Q({'!x1'}, 2);
    [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q);
    assert(all(andNotNums{1} == 1));
    assert(length(andNums{1}) == 0);

    q = Q({'x2 & !x1'}, 2);
    [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q);
    assert(all(andNums{1} == 2));
    assert(all(andNotNums{1} == 1));

    q = Q({'x2 & !x1', 'x2 & x1'}, [2, -1]);
    [andNums, andNotNums, multiplier] = getVarNumsAndMultiplier(q);
    assert(all(andNums{1} == 2));
    assert(all(andNums{2} == [2, 1]));
    assert(all(andNotNums{1} == 1));
    assert(length(andNotNums{2}) == 0);
    assert(all(multiplier == [2, -1]));
end

function testRemoveContradictions(testCase)
    q = Q('x1', 1);
    qRemoved = removeContradictions(q);
    assert(isequal(q.vars, qRemoved.vars));
    assert(isequal(q.multiplier, qRemoved.multiplier));

    q = Q('x1 & !x1', 1);
    qRemoved = removeContradictions(q);
    assert(isequal(qRemoved.vars{1}, 'T'));
    assert(isequal(qRemoved.multiplier, 0));

    q = Q({'x1 & !x1', 'x2 & x1'}, [1, 2]);
    qRemoved = removeContradictions(q);
    assert(isequal(qRemoved.vars{1}, 'x2 & x1'));
    assert(isequal(qRemoved.multiplier(1), 2));
    assert(length(qRemoved.vars) == 1);
end
