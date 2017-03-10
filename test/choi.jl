@testset "Choi with $solver" for solver in sdp_solvers
  @polyvar x y z

  m = Model(solver = solver)

  C = [x^2+2y^2 -x*y -x*z;
       -x*y y^2+2z^2 -y*z;
       -x*z -y*z z^2+2x^2]

  @polyconstraint m C >= 0

  status = solve(m)
  @test status == :Infeasible
end
