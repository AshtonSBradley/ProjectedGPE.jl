δ(j::Integer,k::Integer) = j == k ? 1 : 0

function ladderops(M,a)
  N=0:M-1
  X  = [δ(n,m+1)*√n + δ(n,m-1)*√(n+1) for n in N, m in N]a/√2
  Px = [δ(n,m+1)*√n - δ(n,m-1)*√(n+1) for n in N, m in N]im/(√2a)
  return X,Px
end
