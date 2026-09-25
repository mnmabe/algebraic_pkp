"---- Parameters ------";
//n:= Code length
//k:= Code dimension
//h:= number of blocks
N:=h*Floor(n/h); //n for truncation method 
K:=k+N-n;//k for truncation method 
b:=Floor(n/h); //length of blocks
l:=Floor((N-K)/b);

"---- Complexity Pieces -----";

function Prob(f,u,g,b,h)
  P:=((1-(u+1)/b)^f)*((1-(u/b))^(h-g-f));
  return P;
end function;

function Comp1(N,K,l,b,u,n0)
  T1:=(N-K)*(N-K-(l*(b-u-1)))*(n0+N-K-(l*(b-u-1)));
  return T1;
end function;

function Comp2(n0)
  T2:=((n0^2 + 3*n0)/2)^2.81;
  return T2;
end function;

function Comp3(z,v,w)
  T3:=((5*z)/2)^2.81 + ((5*z)/2)*(v + 2*w + w*v + Binomial(w,2));
  return T3;
end function;

function Comp4(b,g)
  T4:=b^g;
  return T4;
end function;


"----- Finding Best f,u,g for CompTotal (time) ---------";

function timecomp(N,K,h,b,l)
    min_val := Infinity();
    best_u := 0;
    best_f := 0;
    best_g := 0;

    for u in [0..(b-2)] do
        for f in [0..h] do
            for g in [0..(h-f)] do

              n0:= K-f-((h-g)*u);
              if n0 lt 0 then
                    continue;
              end if;

              w:= g*b;
              v:=n0-w;

              if v lt 0 then
                    continue;
              end if;

              // number of quadratic equations
              m1 := f*Binomial(b-u-1,2) + (h-g-f)*Binomial(b-u,2);

              // number of monomials after partial enumeration
              z:=v + 2*g + g*v + Binomial(g,2);

              if z le 0 then
                    continue;
              end if;

              //condition based on m1
              // m1 - v(v+1)/2 >= 5z/2
              if 2*m1 - v*(v+1) lt 5*z then
                  continue;
              end if;

              y:=(N-K-(l*(b-u-1)));

              if y le 0 then
                    continue;
              end if;

              C:= (1.0*Prob(f,u,g,b,h))^-1*(Comp1(N,K,l,b,u,n0) + Comp2(n0) + Comp3(z,v,w)*Comp4(b,g));

              val := Log(2,C);

              if val lt min_val then
                  min_val := val;
                  best_u := u;
                  best_f := f;
                  best_g := g;
                  end if;
            end for;
        end for;
    end for;
    return min_val, best_u, best_f, best_g;
end function;

d, u, f, g := timecomp(N,K,h,b,l);