
function apodize, x, frac

  nx = x
  n = n_elements(nx)
  n1 = round(n*frac)
  avg = average([nx[0:n1-1], nx[n-n1:n-1]])
  sinu = sin(0.5*!dpi*findgen(n1)/n1)
  nx[0:n1-1] = avg + (nx[0:n1-1] - avg)*sinu
  cosi = cos(0.5*!dpi*findgen(n1)/n1)
  nx[n-n1:n-1] = avg + (nx[n-n1:n-1] - avg)*cosi

  return, nx

end
