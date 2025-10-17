function minValue = findMinValue(f, a, b)
    options = optimset('fminbnd');
    minValue = fminbnd(f, a, b, options);
end
