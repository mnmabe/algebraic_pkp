// Generate PKP problem instance
function GenPKPInstance(n, m, q, l)
  A := RandomMatrix(GF(q), m, n);
  N := Nullspace(Transpose(A)); //Kernel of A^T
  d := Dimension(N);

// --------------------------------------------------------
// Build an ell x n matrix C0 whose rows are in the nullspace.
// This version chooses ell linearly independent rows.
// --------------------------------------------------------

 B := Basis(N);

 C := Matrix(GF(q), l , n, &cat[Eltseq(B[i]) : i in [1..l]]);

  S := SymmetricGroup(n);
  // Construct the permutation matrix corresponding to p
  P := PermutationMatrix(GF(q), Random(S));
  while P eq IdentityMatrix(GF(q),n) do
    P := PermutationMatrix(GF(q), Random(S));
  end while;

  C_P := C*P;
  
  // Over F_q we don't want Weight(c) but the sum of the c_i.
  //t := &+[C_P[1,i] : i in [1..n]];
  t := Vector(GF(q), [&+[C_P[i,j] : j in [1..n]]: i in [1..l]]);
  
  // Also return the secret permutation for verification
  return A, C_P, t, Transpose(P), C;
end function;
