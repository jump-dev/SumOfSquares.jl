facts("Non-symmetric matrix SOS constraint") do
    @polyvar x
    m = Model()
    @fact_throws ArgumentError addpolynonnegativeconstraint(m, [1 x; -x 0], BasicSemialgebraicSet())
end
