"--------- #Polys and #Vars for ESDP CCMMP2 --------------";
l := Floor(Log(2,t))+1;

function ESDPnumpolys(n,m,l)
  return 2*n*l+2*n+l+m;
end function;

function ESDPQuadPolys(n,l)
  return n+2*n*l-l;
end function;

function ESDPLinPolys(n,l,m)
  return n+2*l+m-1;
end function;

function ESDPnumvars(n,l)
  return n*(l+1);
end function;


"--------- #Polys and #Vars for ESDP CCMMP4 --------------";

r:=Ceiling(Log(q,n^2-n))+1;

function ESDPnumpolys(n,m,l,r)
  return n^2+2*n+m*l+n^2*r+r;
end function;

function ESDPQuadPolys(n,r)
  return n^2 + r*(n^2-1);
end function;

function ESDPLinPolys(n,m,l,r)
  return 2*n+m*l+2*r;
end function;

function ESDPnumvars(n,r)
  return n^2*(r+1);
end function;

function ESDPDim(n,m,l)
  return n^2-(2*n-1+m*l);
end function;

function ESDPrank(n,m,l)
  return (2*n-1+m*l);
end function;


"--------- #Polys and #Vars for RegSDP BO --------------";

function RegSDPnumpolys(n,m,l)
  return n^2*(n+1)/2 + n + m*l;
end function;

function RegSDPQuadPolys(n)
  return n^2*(n+1)/2
end function;

function RegSDPLinPolys(n,m,l)
  return n+m*l;
end function;

function RegSDPnumvars(n)
  return n^2;
end function;

function RegSDPDim(n,m,l)
  return n^2-(n+m*l);
end function;

function RegSDPrank(n,m,l)
  return (n+m*l);
end function;

"--------- #Polys and #Vars for PKP Direct Modeling --------------";

function PKPnumPolys(n,m,l)
  return n^3 + 2*n + m*l;
end function;

function PKPQuadPolys(n)
  return n^3;
end function;

function PKPLinPolys(n,m,l)
  return 2*n + m*l ;
end function;

function PKPnumvars(n)
  return n^2;
end function;

