"------ Hybrid BO Eq (9) and (10) Degree of Regularity --------";

function SemiReg9(n,k,B,u,h,f)
  P<t>:=PowerSeriesRing(Integers(),4000);
  H:=(((1-t)^(n-k))*((1+(B-u)*t*(1-t)^-1)^f)*((1+B*t*(1-t)^-1)^(h-f)))*(1-t)^-1;
  b:=0;
  while Coefficient(H,b) gt 0 do
      b:=b+1;
  end while;
  return b;
end function;

//---- Monomial count series used in the complexity estimate ---

function MonomialCount9(B,u,h,f,d)
    P<t> := PowerSeriesRing(Integers(), 4000);

    Mons := (1 + (B-u)*t*(1-t)^-1)^f*(1 + B*t*(1-t)^-1)^(h-f);

    mc := &+[Coefficient(Mons,i) : i in [0..d]];
    
    return mc;
end function;


"------ finding the best cost with corresponding u and f -----";
//n:=length of code word
//k:=dimension of code
//h:=number of blocks
//B:=Floor(n/h) length of block

function dreg9(n,k,B,h)
    best_cost := Infinity();
    best_d := 0;
    best_u := 0;
    best_f := 0;
    best_mc := 0; 

    for u in [1..B-1] do
        for f in [1..h] do
            if f*u lt k then

                d := SemiReg9(n,k,B,u,h,f);

                // Ignore non-useful degree-1 cases
                if d ne 1 then

                mc := MonomialCount9(B,u,h,f,d);

                // Success probability
                Prob := (1 - u/B)^f;

                    if Prob gt 0 then 
                        raw_cost := 3 * (1/Prob) * (k + 1 - f*u) * (mc)^2;
                        cost := Log(raw_cost)/Log(2);


                        if cost lt best_cost then
                            best_cost := cost;
                            best_d := d;
                            best_u := u;
                            best_f := f;
                            best_mc := mc;
                        end if; 
                    end if; 
                end if; 
            end if;
        end for;
    end for;
    return best_cost, best_d, best_u, best_f, best_mc;
end function;

C, d, u, f, mc := dreg9(n,k,B,h);


"------ Hybrid BO Eq (9) and (10) for Binary Semi Reg Degree of Regularity --------";

function SemiReg2(n,k,B,u,h,f)
  P<t>:=PowerSeriesRing(Integers(),4000);
  H:=(((1+t)^(-(n-k)))*((1+(B-1-u)*t)^f)*((1+(B-1)*t)^(h-f)))*(1-t)^-1;
  b:=0;
  while Coefficient(H,b) gt 0 do
      b:=b+1;
  end while;
  return b;
end function;

"---- Monomial-count series used in the complexity estimate ---";

function MonomialCount2(B,u,h,f,d)
    P<t> := PowerSeriesRing(Integers(), 4000);

    Mons := (1+(B-1-u)*t)^f*(1+(B-1)*t)^(h-f);

    mc := &+[Coefficient(Mons,i) : i in [0..d]];
    
    return mc;
end function;


"------ finding the minimum with corresponding u and f -----";
//n:=length of code word
//k:=dimension of code
//h:=number of blocks
//B:=Floor(n/h) length of block

function dreg2(n,k,B,h)
    best_cost := Infinity();
    best_d := 0;
    best_u := 0;
    best_f := 0;
    best_mc := 0; 

    for u in [1..B-1] do
        for f in [1..h] do
            if f*u lt k then

                d := SemiReg2(n,k,B,u,h,f);

                // Ignore non-useful degree-1 cases
                if d ne 1 then

                mc := MonomialCount2(B,u,h,f,d);

                // Success probability
                Prob := (1 - u/B)^f;

                    if Prob gt 0 then 
                        raw_cost := 3 * (1/Prob) * (k + 1 - f*u) * (mc)^2;
                        cost := Log(raw_cost)/Log(2);


                        if cost lt best_cost then
                            best_cost := cost;
                            best_d := d;
                            best_u := u;
                            best_f := f;
                            best_mc := mc;
                        end if; 
                    end if; 
                end if; 
            end if;
        end for;
    end for;
    return best_cost, best_d, best_u, best_f, best_mc;
end function;

C, d, u, f, mc := dreg2(n,k,B,h);


"------ Hybrid BO Eq (8) and (10) Semi Reg Degree of Regularity --------";

function SemiReg8(n,k,B,u,h)
  P<t>:=PowerSeriesRing(Integers(),4000);
  H:= ((1-t)^(n-k)*(1+(B-u)*t*(1-t)^-1)^h)*(1-t)^-1;
  b:=0;
  while Coefficient(H,b) gt 0 do
      b:=b+1;
  end while;
  return b;
end function;

"---- Monomial-count series used in the complexity estimate ---";

function MonomialCount8(B,u,h,f,d)
    P<t> := PowerSeriesRing(Integers(), 4000);

    Mons := (1 + (B-u)*t*(1-t)^-1)^f*(1 + B*t*(1-t)^-1)^(h-f);

    mc := &+[Coefficient(Mons,i) : i in [0..d]];
    
    return mc;
end function;


"------ finding the minimum with corresponding u and f -----";
//n:=length of code word
//k:=dimension of code
//h:=number of blocks
//B:=Floor(n/h) length of block

function dreg8(n,k,B,h)
    f:=h;
    best_cost := Infinity();
    best_d := 0;
    best_u := 0;
    best_mc := 0; 

    for u in [1..B-1] do

        if f*u lt k then
            
        d := SemiReg8(n,k,B,u,h);
            
            // Ignore non-useful degree-1 cases
            if d ne 1 then

            mc := MonomialCount8(B,u,h,f,d);

            // Success probability
             Prob := (1 - u/B)^f;

                if Prob gt 0 then 
                    raw_cost := 3 * (1/Prob) * (k + 1 - f*u) * (mc)^2;
                    cost := Log(raw_cost)/Log(2);


                    if cost lt best_cost then
                        best_cost := cost;
                        best_d := d;
                        best_u := u;
                        best_mc := mc;
                    end if; 
                end if; 
            end if; 
         end if;
    end for;
    return best_cost, best_d, best_u, best_mc;
end function;

C, d, u, mc := dreg8(n,k,B,h);
