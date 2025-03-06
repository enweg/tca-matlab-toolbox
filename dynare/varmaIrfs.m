function irfs=varmaIrfs(Phi0, As, Psis, horizon)
  p = length(As);
  q = length(Psis);
  n = size(Phi0, 1);

  % calculating irfs
  irfs = zeros(n, n, horizon+1);
  irfs(:, :, 1) = Phi0;
  for h=1:horizon
    for i=1:min(p, h)
      irfs(:, :, h+1) = irfs(:, :, h+1) + As{i}*irfs(:, :, h-i+1);
    end
    if h <= q
      irfs(:, :, h+1) = irfs(:, :, h+1) + Psis{h} * Phi0;
    end
  end

end
