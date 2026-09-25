// ------------------------------------------------------------
// Construct Hhat and sstar
//
// Input:
//   A : n x m matrix over F_q
//   C : ell x n matrix over F_q, where C[k,r] = c_{k,r}
//
// Output:
//   Hhat  : (n + ell*m) x n^2 matrix over F_q
//   sstar : (n + ell*m) x 1 column matrix over F_q
// ------------------------------------------------------------

function GenParityCheck(n,m,l,A,C)
  F := BaseRing(A);
  numRows := n + l*m;
  numCols := n^2;

  Hhat := ZeroMatrix(F, numRows, numCols);

  // --------------------------------------------------------
  // U in F_q^{n x n^2} 
  //
  // u_{i,j} = 1 if (i-1)n < j <= i n,
  //           0 otherwise.
  //U is NOT needed for RegSDP 
  // --------------------------------------------------------

  //for i in [1..n] do
      //for j in [((i-1)*n + 1)..(i*n)] do
          //Hhat[i,j] := 1;
      //end for;
  //end for;

  // --------------------------------------------------------
  // V = (I_n || I_n || ... || I_n) in F_q^{n x n^2}
  //
  // There are n copies of I_n.
  // These rows are placed in Hhat from n+1 to 2n.
  // --------------------------------------------------------

  for block in [1..n] do //each block of an identity matrix 
      for i in [1..n] do //help iterate rows of each identity matrix 
          j := (block-1)*n + i; //columns of idenity matrix 
          Hhat[i,j] := 1;
      end for;
  end for;


  // --------------------------------------------------------
  // W_k in F_q^{m x n^2}
  //
  // w_{i,j}^{(k)} = c_{k, floor((j-1)/n)+1}
  //                 * A_{((j-1) mod n)+1, i}.
  //
  // The W_k blocks are stacked below U and V.
  // --------------------------------------------------------

  for k in [1..l] do
    rowOffset := n + (k-1)*m;//making sure we get ml rows but start at n+1 
      for i in [1..m] do
          for j in [1..numCols] do

              cIndex := Floor((j-1)/n) + 1;
              aRow   := ((j-1) mod n) + 1; // check indexing here with j v j-1 

              Hhat[rowOffset + i, j] := C[k,cIndex] * A[aRow,i];

          end for;
      end for;
  end for;


    // --------------------------------------------------------
    // sstar = (1,...,1,0,...,0)^T
    //
    // First 2n entries are 1.
    // Last ell*m entries are 0.
    // --------------------------------------------------------

    sEntries := [(r le n) select 1 else 0 : r in [1..numRows]];

    sstar := Matrix(F, numRows, 1, sEntries);

  return Hhat, sstar;
end function; 


function ParityCheckEquations(P, H, C, s, n)
  Fqx := Parent(P[1,1]);
  H := ChangeRing(H, Fqx);
  C := ChangeRing(C, Fqx);
  s := ChangeRing(s, Fqx);
  pvec := Matrix(BaseRing(P), 1, n^2, Eltseq(P));

  eqs := pvec*Transpose(H) - Transpose(s);

  //eqs := H*Transpose(C*P) - Transpose(s);

  return [eqs[1][i] : i in [1..Ncols(eqs)]];
end function;

function FieldEquations(P)
  return [x^2 - x : x in Eltseq(P)];
end function;


function RegularFormEquations(P)
n := Nrows(P);
R1 := [];

   for i in [1..n] do
    for j in [1..n] do
      for k in [1..n] do
        if j lt k then
          f := P[i][j]*P[i][k];
          if not f in R1 then
            R1 := Append(R1, f);
          end if;
        end if;
      end for;
    end for;
  end for;
return R1;
end function; 

//----------------------------
//we have n^2 x vals 
//we need to randomly assign u coordinates as 0 every n 

function ErrorFreeIndices(n,u)

I:=[];

for k in [1..n] do //how many blocks we have and would like to do this for 
  R := RandomSubset({((k-1)*n+1)..(k*n)}, u); 
  I cat:= [j : j in [((k-1)*n + 1)..(k*n)] | j in R];
  //I := Append(I,R); 
end for; 

return I; 
end function; 

function ErrorFreeSystem(P,u,I)
pvec := Eltseq(P);
S0:= [pvec[i] : i in I];

return S0; 
end function; 


function HBOSystem(P,H,C,s,n)
  L1:= ParityCheckEquations(P, H, C, s, n); 
  Q1:= FieldEquations(P);
  Q2:= RegularFormEquations(P);

  S:= L1 cat Q1 cat Q2;

  return S;
end function; 


function HilbertSeriesQL(S)

  HS := HilbertSeriesIdeal(Top(S));
    
  return HS;
end function;

function DReg(HS)
   d:= Degree(HS)+1; 

   return d;
end function;