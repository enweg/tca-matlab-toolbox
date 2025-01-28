function tests = utilsTest
    tests = functiontests(localfunctions);
end

function testMapX2YAndBack(testCase)
    order = [3, 1, 2];
    strY = "y_{1, 2} & y_{3, 1}";
    strX = mapY2X(strY, order);
    assert(isequal(strY, mapX2Y(strX, order)));

    strY2 = "y_{1, 2} & y_{3, 1}";
    strX = mapY2X(strY2, order);
    assert(isequal(strY, mapX2Y(strX, order)));

    order = [3, 1, 2];
    strY = "!y_{1, 2} & y_{3, 1}";
    strX = mapY2X(strY, order);
    assert(isequal(strY, mapX2Y(strX, order)));

    order = [3, 1, 2];
    strY = "!y_{1, 0} & y_{3, 1}";
    strX = mapY2X(strY, order);
    assert(isequal(strY, mapX2Y(strX, order)));

    order = [3, 1, 2];
    strY = "!(y_{1, 0} & y_{3, 1})";
    strX = mapY2X(strY, order);
    assert(isequal(strY, mapX2Y(strX, order)));
end

function testSlideIn(testCase)
    A = [];
    for i = 3:-1:1
        A = [A, repmat(i, 3, 3)];
    end
    B = zeros(12, 12);
    B = slideIn(B, A);
    assert(all(B(1:3, 1:3) == 1, 'all'));
    assert(all(B(1:3, 4:end) == 0, 'all'));
    assert(all(B(4:6, 1:3) == 2, 'all'));
    assert(all(B(4:6, 4:6) == 1, 'all'));
    assert(all(B(4:6, 7:end) == 0, 'all'));
    assert(all(B(7:9, 1:3) == 3, 'all'));
    assert(all(B(7:9, 4:6) == 2, 'all'));
    assert(all(B(7:9, 7:9) == 1, 'all'));
    assert(all(B(7:9, 10:end) == 0, 'all'));
    assert(all(B(10:end, 1:3) == 0, 'all'));
    assert(all(B(10:end, 4:6) == 3, 'all'));
    assert(all(B(10:end, 7:9) == 2, 'all'));
    assert(all(B(10:end, 10:end) == 1, 'all'));
end
