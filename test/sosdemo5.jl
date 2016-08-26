# Adapted from:
# SOSDEMO5 --- Upper bound for the structured singular value mu
# Section 3.5 of SOSTOOLS User's Manual

facts("SOSDEMO5") do
for solver in sdp_solvers
context("With solver $(typeof(solver))") do
  @polyvar x1 x2 x3 x4 x5 x6 x7 x8
  vartable = [x1, x2, x3, x4, x5, x6, x7, x8]

  # The matrix under consideration
  alpha = 3 + sqrt(3);
  beta = sqrt(3) - 1;
  a = sqrt(2/alpha);
  b = 1/sqrt(alpha);
  c = b;
  d = -sqrt(beta/alpha);
  f = (1 + im)*sqrt(1/(alpha*beta));
  U = [a 0; b b; c im*c; d f];
  V = [0 a; b -b; c -im*c; -im*f -d];
  M = U*V';

  # Constructing A(x)'s
  A = Vector{Matrix{Float64}}(4)
  gam = 0.8724;

  Z = monomials(vartable, 1)
  for i = 1:4
      H = M[i,:]*M[i,:]' - (gam^2)*sparse([i],[i],[1],4,4)
      H = [real(H) -imag(H); imag(H) real(H)]
      A[i] = dot(Z, H*Z)
  end


  m0 = JuMP.Model(solver = solver)

  # -- Q(x)'s -- : sums of squares
  # Monomial vector: [x1; ... x8]
  Q = Vector{MatPolynomial{JuMP.Variable}}(4)
  for i = 1:4
    @SOSvariable m tmp >= 0 Z
    Q[i] = tmp
  end

  # -- r's -- : constant sum of squares
  Z = monomials(vartable, 0)
  r = Matrix{MatPolynomial{JuMP.Variable}}(4,4)
  for i = 1:4
      for j = (i+1):4
        @SOSvariable m tmp >= 0 Z
        r[i,j] = tmp
      end
  end

  # Constraint : -sum(Qi(x)*Ai(x)) - sum(rij*Ai(x)*Aj(x)) + I(x) >= 0
  expr = 0
  # Adding term
  for i = 1:4
      expr -= A[i]*Q[i]
  end
  for i = 1:4
      for j = (i+1):4
          expr -= A[i]*A[j]*r[i,j]
      end
  end
  # Constant term: I(x) = -(x1^4 + ... + x8^4)
  I = -sum(vartable.^4)
  expr = expr + I

  @SOSconstraint m expr >= 0

  status = solve(m)

  # Program is feasible, thus 0.8724 is an upper bound for mu.
  @fact status --> :Optimal
end; end; end
