// Compute system of polynomial equations in F_q[x_1, .. x_{n^2}] whose solution
// is a permutation matrix P


function ParityCheckEquations(P, H, c, s)
  Fqx := Parent(P[1,1]);
  H := ChangeRing(H, Fqx);
  c := ChangeRing(c, Fqx);
  s := ChangeRing(s, Fqx);

  eqs := H*Transpose(c*P) - Transpose(s);
  return [eqs[i][1] : i in [1..Nrows(eqs)]];
end function;

function FieldEquations(P)
  return [x^2 - x : x in Eltseq(P)];
end function;

// Separate linear and quadratic eqns for testing
function LinearPermEquations(P)
  n := Nrows(P);
  L1 := [];
  L2 := [];
  
  for i in [1..n] do
    // the first n linear equations of the 1st perm equations
    f := &+[P[i][j] : j in [1..n]] - 1;
    L1 := Append(L1, f);
    
    // the other half of n linear equations of the 1st perm equations
    f := &+[P[j][i] : j in [1..n]] - 1;
    L2 := Append(L2, f);
  end for;

  return L1, L2;
end function;

function QuadraticPermEquations(P)
  n := Nrows(P);
  Q1 := [];
  Q2 := [];

  for i in [1..n] do
    for j in [1..n] do
      for k in [1..n] do
        if j ne k then
          f := P[i][j]*P[i][k];
          if not f in Q1 then
            Q1 := Append(Q1, f);
          end if;
        end if;

        if i ne k then    
          f := P[i][j]*P[k][j];
          if not f in Q2 then
            Q2 := Append(Q2, f);
          end if;
        end if;
      end for;
    end for;
  end for;
  Q3 := FieldEquations(P);

  return Q1, Q2, Q3;
end function;

function PermEquations(P)
  L1, L2 := LinearPermEquations(P);
  Q1, Q2, Q3 := QuadraticPermEquations(P);
  S := L1 cat L2 cat Q1 cat Q2 cat Q3;
  //S := L cat Q;
  return S;
end function;

function System(P, H, c, s)
  S := ParityCheckEquations(P, H, c, s);
  S cat:= PermEquations(P);
  return S;
end function;

function GenSystem(n, m, q)
  Fqx<[x]> := PolynomialRing(GF(q), n^2, "grevlex");

  A, c, t, P0 := GenPKPInstance(n, m, q); //this gives us our givens here for the instance
  H:= VerticalJoin(A, Matrix(GF(q), 1, Ncols(A), [1: i in [1.. Ncols(A)]]));
  s:=HorizontalJoin(ZeroMatrix(GF(q),1,m), Matrix(GF(q),1,1,[t])); // PKP s here

  P := Matrix(Fqx, n, n, x);
  S := System(P, H, c, s);
  return S;
end function;
