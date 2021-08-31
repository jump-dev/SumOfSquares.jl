const EXAMPLES = filter(ex -> endswith(ex, ".jl") && ex != "run_examples.jl" && ex != "all_examples.jl" && ex != "goldsteinprice.jl" && ex != "sparse_polynomials.jl" && ex != "scaled_perm.jl" && ex != "symmetry.jl" && ex != "discourse_65377.jl" && ex != "discourse_65377_1.jl" && ex != "prajna.jl",
                        readdir(@__DIR__))
