function v = dynareCellArrayToVec(ca)
  v = repelem("", length(ca));
  for i=1:length(ca)
    v(i) = string(cell2mat(ca(i)));
  end
  v = vec(v);
end
