function HilbertSeriesQL(n, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);
  
  Q1, Q2, Q3 := QuadraticPermEquations(P);
  L1, L2 := LinearPermEquations(P);
  F := Q1 cat Q2 cat Q3 cat L1 cat L2;

  HS := HilbertSeriesIdeal(Top(F));
    
  return HS;
end function;

function DReg(HS)
   d:= Degree(HS)+1; 

   return d;
end function;