pathnodes = zeros(Int,l+1)
pathnodes[l+1] = curnode
pathedges = zeros(Int,l)
pathprods = cell(1,l+1)
pathprods[l+1] = eye(n)
while i <= l
  if i == 0
    npaths = npaths + 1
    pathprod = pathprods[1]
    cur = dot(lyap.dual[pathedges[1]], p_k(pathprod*x, x))
    if cur > best
      best = cur
      best_e = pathedges
    end
    i = i + 1
  else
    pathedges[i] = pathedges[i] + 1
    if pathedges[i] > nEdges
      pathedges[i] = 0
      i = i + 1
    elseif edges(pathedges[i], 2) == pathnodes[i+1]
      pathnodes[i] = edges(pathedges[i],1)
      pathprods[i] = pathprods[i+1] * getmat(s, pathedges[i])
      i = i - 1
    end
  end
end

