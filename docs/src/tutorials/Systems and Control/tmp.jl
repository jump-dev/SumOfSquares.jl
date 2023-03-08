Vs = [V - γ]
using ImplicitPlots
using Plots
function plot_lyapunovs(Vs, J; resolution = 1000)
    p = plot()
    eliminated = x[setdiff(1:n_x, J)]
    for i in eachindex(Vs)
        V = subs(Vs[i], eliminated => zeros(length(eliminated)))
        @show V
        @show variables(V)
        implicit_plot!(p, V; resolution, label="V$i")
    end
    return p
end
plot_lyapunovs(Vs, [5, 6])
