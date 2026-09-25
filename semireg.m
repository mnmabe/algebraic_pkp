function SemiReg(M,D,N)
  P<t>:=PowerSeriesRing(Integers(),40000);
  f:=(&*[(1-t^D[i])^M[i] : i in [1..#M]])*(1-t)^-N;
  b:=0;
  while Coefficient(f,b) gt 0 do
      b:=b+1;
  end while;
  return b;
end function;

//SemiReg([Q,L],[2,1],N);