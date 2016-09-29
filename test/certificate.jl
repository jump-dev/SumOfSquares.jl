facts("Monomial selection for certificate") do
  @polyvar x y
  @fact_throws ErrorException getmonomialsforcertificate([x*y, y^2], :Sparse)
  @fact getmonomialsforcertificate([x*y, y^2]) --> MonomialVector([y])
end
