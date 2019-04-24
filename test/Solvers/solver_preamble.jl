# In `solver_tests.jl` we want to load the preamble only once but if we load
# several time, say `csdp_tests.jl`, from the REPL it may be because we modify
# the tests in-between so we want to reload it
if !(@isdefined shared_preamble) || !shared_preamble
    include("../Tests/Tests.jl")
    using Test, JuMP
end
