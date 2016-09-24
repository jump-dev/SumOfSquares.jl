facts("Substitution") do
  @polyvar x[1:3]
  @fact (x[1])([x[2]], [x[1]]) == x[2] --> true
  p = x[1]*x[2]*x[3]
  @fact Int(p([1, 2, 3], x)) --> 6
  p = x[1]^2 + x[1]*x[3] - 3
  @fact Int(p([5, x[1]*x[2], 4], x)) --> 42
  p = x[1]^2 + x[2]^2
  q = p([1 -1; 1 1] * x[1:2], x[1:2])
  @fact q == 2p --> true
end
