// Table 1 of PKP_as_SDP_and_MQ (September 2026 draft).
// Run: magma -b table1.m
// To update both tables with official Python TII estimates:
//   magma -b table1.m > table1.txt
//   python scripts/update-table1-tii.py /path/to/CryptographicEstimators
// CCMMP2 is valid on n coordinates only for q=2, ell=1; otherwise use
// the permutation-vector ESDP (length n^2, weight n) and CCMMP4.
// L2 is the actual linear rank; Q2 counts the original quadratic equations.
// For CCMMP4, Q2=N+(N-1)*r2 and L2=rank(H)+2*r2: the first recurrence
// is linear.
// d_sol is the observed maximum input/F4 step degree, using grevlex with
// x first, then y in row order; it can vary with the instance and Magma.

function Table1System(n, m, ell, q)
    K := GF(q);
    repeat A := RandomMatrix(K, m, n); until Rank(A) eq m;
    B := BasisMatrix(Nullspace(Transpose(A)));
    repeat U := RandomMatrix(K, ell, n-m); until Rank(U) eq ell;
    P := PermutationMatrix(K, Random(SymmetricGroup(n)));
    C := U*B*P;
    secret := Transpose(P);
    assert IsZero(C*secret*Transpose(A));

    if q eq 2 and ell eq 1 then
        N := n;
        witness := Eltseq(C*secret);
        t := Weight(C[1]); // Integer Hamming weight, not the sum in GF(2).
        b := #Intseq(t, 2);
        r := 0; // CCMMP2 has no r_1,r_2.
        R := PolynomialRing(K, N*(b+1), "grevlex");
        x := [R.i : i in [1..N]];
        y := Matrix(R, N, b, [R.i : i in [N+1..Rank(R)]]);
        H := VerticalJoin(A, Matrix(K, 1, N, [1 : i in [1..N]]));
        syndrome := [K!0 : i in [1..m]] cat [K!t];
        F := [R.i^2-R.i : i in [1..Rank(R)]];
        F cat:= [y[1,1]-x[1]] cat [y[1,j] : j in [2..b]];
        F cat:= [x[i]+y[i,1]+y[i-1,1] : i in [2..N]];
        F cat:= [y[i,j-1]*y[i-1,j-1]+y[i,j]+y[i-1,j-1]+y[i-1,j]
                 : i in [2..N], j in [2..b]];
        F cat:= [y[i,b]*y[i-1,b]+y[i-1,b] : i in [2..N]];
        F cat:= [y[N,j]-K!((t div 2^(j-1)) mod 2) : j in [1..b]];
        w := 0;
        for i in [1..N] do
            w +:= Integers()!witness[i];
            witness cat:= [K!((w div 2^(j-1)) mod 2) : j in [1..b]];
        end for;
        bound := (N div 2) + (b-1)*((N-1) div 2) + 1;
    else
        N := n^2;
        t := n;
        r := 1;
        while q^r le Max(t, N-t)+1 do r +:= 1; end while;
        R := PolynomialRing(K, N*(r+1), "grevlex");
        x := [R.i : i in [1..N]];
        y := Matrix(R, N, r, [R.i : i in [N+1..Rank(R)]]);
        ones := Matrix(K, 1, n, [1 : i in [1..n]]);
        H := VerticalJoin(KroneckerProduct(IdentityMatrix(K,n), ones),
                          KroneckerProduct(ones, IdentityMatrix(K,n)));
        for a in [1..ell] do
            H := VerticalJoin(H, KroneckerProduct(Matrix(K,1,n,Eltseq(C[a])), A));
        end for;
        syndrome := [K!1 : i in [1..2*n]] cat [K!0 : i in [1..m*ell]];
        // A primitive companion matrix gives period q^r-1 > max(t,N-t).
        M := Transpose(CompanionMatrix(PrimitivePolynomial(K,r)));
        y0 := Matrix(K, r, 1, [1] cat [0 : j in [2..r]]);
        prev := ChangeRing(y0, R);
        F := [u^2-u : u in x];
        for i in [1..N] do
            curr := Transpose(Submatrix(y, i, 1, 1, r));
            F cat:= Eltseq(curr-(1-x[i])*prev-x[i]*ChangeRing(M,R)*prev);
            prev := curr;
        end for;
        F cat:= Eltseq(prev-ChangeRing(M^t*y0,R));
        witness := Eltseq(secret);
        state := y0;
        for i in [1..N] do
            if witness[i] eq 1 then state := M*state; end if;
            witness cat:= Eltseq(state);
        end for;
        bound := 0; // No CCMMP2 bound for CCMMP4.
    end if;
    // Reduce the syndrome equations together with their right-hand side.
    aug := EchelonForm(HorizontalJoin(H, Matrix(K,#syndrome,1,syndrome)));
    rho := Rank(H);
    assert Rank(aug) eq rho;
    F cat:= [&+[aug[i,j]*x[j] : j in [1..N]]-aug[i,N+1] : i in [1..rho]];
    assert &and[Evaluate(f,witness) eq 0 : f in F];
    return F, t, r, bound, A, C;
end function;

function Table1Degrees(F, r : Solve := true)
    R := Universe(F);
    nv := Rank(R);
    Q := #[f : f in F | TotalDegree(f) eq 2];
    linear := [f : f in F | TotalDegree(f) eq 1];
    // L2 counts independent linear equations, including weight constraints.
    L := Rank(Matrix(BaseRing(R),
                     [[MonomialCoefficient(f,R.i) : i in [1..nv]] : f in linear]));
    // Formal semi-regular prediction: first nonpositive coefficient of
    // (1-z^2)^Q/(1-z)^(nv-L). It is a heuristic, not a bound.
    S<z> := PowerSeriesRing(Rationals(), nv+2);
    h := (1-z^2)^Q/(1-z)^(nv-L);
    sr := 0;
    while Coefficient(h,sr) gt 0 do sr +:= 1; end while;

    if r eq 0 then
        I := ideal<R | [HomogeneousComponent(f,TotalDegree(f)) : f in F]>;
        assert Dimension(I) eq 0;
        hs := HilbertSeries(quo<R|I>);
        assert Degree(Denominator(hs)) eq 0;
        dreg := Sprint(Degree(Numerator(hs))+1);
    else
        // In CCMMP4, setting x=0, y_1=y_N=0 annihilates every top form.
        // The (N-2)*r interior y variables are free, so d_reg is not finite.
        dreg := "inf";
    end if;
    // Magma handbook: GroebnerBasis(S), second output = F4 step degrees.
    // https://magma.maths.usyd.edu.au/magma/handbook/groebner_bases
    dsol := 0; // Solve=false lets a timed caller save the other columns first.
    if Solve then
        G, steps := GroebnerBasis(F);
        dsol := Max([TotalDegree(f) : f in F] cat steps);
    end if;
    return nv, Q, L, sr, dreg, dsol;
end function;

procedure Table1(: Seed := 1)
    print "n m ell q | model t r1 r2 | n2 Q2 L2 | d_bound d_SR d_reg d_sol";
    for q in [2,5] do
        for ell in [1,2] do
            for n in [3..5] do
                m := n div 2;
                SetSeed(Seed+10000*q+100*n+10*m+ell);
                F, t, r, bound := Table1System(n,m,ell,q);
                nv, Q, L, sr, dreg, dsol := Table1Degrees(F,r);
                printf "%o %o %o %o | %o %o %o %o | %o %o %o | %o %o %o %o\n",
                    n,m,ell,q, r eq 0 select "CCMMP2" else "CCMMP4", t,
                    r eq 0 select "--" else "1", r eq 0 select "--" else Sprint(r),
                    nv,Q,L, bound eq 0 select "--" else Sprint(bound), sr,dreg,dsol;
            end for;
        end for;
    end for;
end procedure;

if not assigned Table1Library then
    Table1();
end if;
