# algebraic_pkp
This is a public repository for Algebraic Modelings of the Permuted Kernel Problem. 

For the tables below, start by generaing an instance of PKP using gen_pkp.m. This will generate a matrix A that is mxn and a codeword C that is lxn with weight t. We also return the permutation P and the permuted codeword C_P for verification purposes.

For each system, we can also find the number of quadratic and linear polynomials as well as the number of variables using the numpolys.m. Note here that the l parameter given for the CCMMP2 modeling is only used for that modeling and is different from the l used in the other modelings which comes from the codeword C. 

Compute the degree of regularity using the semi-regular assumption using semireg.m.

Compute the complexity for some degree estimate using comp.m. 

Generating Table 1: 


Generating Table 2: 
After generating a PKP instance, find the transpose of A, and use gen_BO.m to generate the polynomial system resulting from reducing PKP to a RegSDP instance. This file also gives the actual degree of regularity and solving degree of the system. To find the hybrid witness degrees, use the BO_hdreg.m file. Equations correspond to (8) or (9) respectively, where (8) considers guessing u positions in every vector block, while (9) is more general guessing error free positions in a number of error free blocks. Checking solving degree of hybrid system, solve the generated BO system combined with the generated Error Free System. 

Generating Table 3: 
After generating a PKP instance, use gen_direct.m to generate the polynomial system resulting from a direct algebraic modeling of PKP. This file also gives the actual degree of regularity and solving degree of the system.  

Generating Table 4: 
Table 4 uses the TII estimator for PKP and ESDP. To find the state-of-the-art complexities for RegSDP, use SOA_regSDP_comp.m which will give the best complexity with corresponding parameters (f,u,g). 

Generating Table 5: 
Table 5 then uses the PKP parameters seen in Table 4 to find the number of polynomials and variables, degree estimates, and complexities for the 3 different algebraic approaches. 
