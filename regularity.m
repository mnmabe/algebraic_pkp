// Compute Hilbert series and degree of regularity under semi-regularity assumption.
// Test semi-regularity.
// Test generic coordinates.
// Compute Castelnuovo-Mumford regularity < Macaulay bound.
// [1]: The complexity of algebraic algorithms for LWE. M. Steiner.

// only quadratic polynomials
function HilbertSeriesQtop(n, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);
  
  Q1, Q2, Q3 := QuadraticPermEquations(P);
  Q := Q1 cat Q2 cat Q3;

  return HilbertSeriesIdeal(Top(Q));
end function;

// only linear polynomials. Note the hilbert series is not a polynomial
// in this case.
function HilbertSeriesLtop(n, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);
  
  L1, L2 := QuadraticPermEquations(P);
  L := L1 cat L2;

  G := Top(L);
  if not IsGroebner(G) then
    print "Ltop is not a GB, computing GB.";
    G := GroebnerBasis(G);
  end if;

  return HilbertSeriesIdeal(G);
end function;
  
// quadratic and linear
function HilbertSeriesQtopLtop(n, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);

  G := Top(PermEquations(P));
  if not IsGroebner(G) then
    print "Qtop union Ltop is not a GB, computing GB.";
    G := GroebnerBasis(G);
  end if;

  return HilbertSeriesIdeal(G);
end function;

// Complete system
function HilbertSeriesSystem(n, m, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);

  A, c, t, P0 := GenPKPInstance(n, m, q);
  H := VerticalJoin(A, Matrix(GF(q), 1, Ncols(A), [1: i in [1.. Ncols(A)]]));
  s := ZeroMatrix(GF(q), 1, m+1);
  s[1,m+1] := t;

  S := System(P, H, c, s);
  if not CheckSystem(S, P0) then
    print "Something went wrong!";
  end if;
  
  G := Top(S);
  if not IsGroebner(G) then
    print "S^top is not a GB, computing GB.";
    G := GroebnerBasis(G);
  end if;

  return HilbertSeriesIdeal(G);
end function;

function DegreeOfRegularityQtop(n, q)
  HS := HilbertSeriesQtop(n, q);
  return Degree(HS) + 1;
end function;

function DegreeOfRegularityQtopLtop(n, q)
  HS := HilbertSeriesQtopLtop(n, q);
  return Degree(HS) + 1;
end function;




// Test if S is in generic coordinates. This implies the solving degree is bounded by
// the CM regularity. S does not need to be groebner basis.
//
// This DOES require the system is defined over an algebraically closed field. We will ignore
// that for now...
function IsInGenericCoordinates(S)
  Fqx := Parent(S[1]);
  x := [Fqx.i : i in [1..Rank(Fqx)]];
  
  I := Ideal(S);

  // Check conditions of [1] Theorem 7:
  if I eq Fqx then
    print "Condition 1 failed: I = (1).";
    return 0;
  elif not Dimension(I) eq 0 then
    print "Condition 2 failed: Dim(I) != 0.";
    return 0;
  end if;

  // these are conditions 2 and 3, which are equivalent to being in generic coordinates.
  // Don't really need both, just for completeness
  if Radical(I) eq Ideal(x) then
    return true;
  elif Dimension(Ideal(Top(S))) eq 0 then
    print "Something went wrong! Condition 2 failed but condition 3 didn't, but they should be equivalent.";
    return true;
  end if;
    
  return false;  
end function;

// Bound solving degree by testing if system is in generic coordinates, then bounding by
// CM regularity which is bounded by the Macaulay bound.
// Corollary 6 of [1].
function BoundSolvingDegreeByCMRegularity(S)
  // Requires S^top is in generic coordinates over the algebraic closure.
  Fqx := Parent(S[1]);
  n := Rank(Fqx);
  Fq := BaseRing(Fqx);
  Fq_bar := AlgebraicClosure(Fq);
  Fq_bary<[y]> := PolynomialRing(Fq_bar, n, "grevlex");

  F_hom := [Fq_bary!f : f in Top(S)];
  if not IsInGenericCoordinates(F_hom) then
    print "Input not in generic coordinates over algebraic closure.";
    return 0;
  end if;

  degrees := Sort(Setseq(Set([TotalDegree(f) : f in S])));

  m := #S;
  l := Min(n+1, m);

  bound := 1;
  // Sort(S) sorts based on ordering (?). A graded ordering means it is ordered with respect to
  // degree first. Reverse it to start with high degree polynomials.
  for f in Reverse(Sort(S))[1..l] do
    bound +:= TotalDegree(f);
  end for;
  bound -:= l;

  d := Max([TotalDegree(f) : f in S]);
  print "Computed bound = ", bound;
  print "(n + 1)(d - 1) + 1 = ", (n+1)*(d-1)+1;

  return bound;
end function;


function test_bound_permeqns(n, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);

  S := PermEquations(P);
  
  bound := BoundSolvingDegreeByCMRegularity(S);
  return bound;
end function;

function test_bound(n, m, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");
  P := Matrix(Fqx, n, n, x);

  A, c, t, P0 := GenPKPInstance(n, m, q);
  H := VerticalJoin(A, Matrix(GF(q), 1, Ncols(A), [1: i in [1.. Ncols(A)]]));
  s := ZeroMatrix(GF(q), 1, m+1);
  s[1,m+1] := t;
  
  S := System(P, H, c, s);
  
  print "System has size =", #S;
  if not CheckSystem(S, P0) then
    print "Something went wrong!";
  end if;
  
  bound := BoundSolvingDegreeByCMRegularity(S);
  return bound;
end function;
