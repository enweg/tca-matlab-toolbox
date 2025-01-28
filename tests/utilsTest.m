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
