import Pajarito
using JuMP
using DynamicPolynomials

model = SOSModel(Pajarito.Optimizer)

@polyvar(t)
monos = monomials([t], 0:r)

struct Box
    xl::Float64
    yl::Float64
    xu::Float64
    yu::Float64
end

boxes = [
    Box(0, 10, 6, 11),
    Box(0,  9, 3, 10),
]

N = 10

T = 1:(N+1)

@variable(model, H[1:N, boxes], Bin)

Mxl = Myl = -20
Mxu = Myu = 20

p = Dict()
for j in 1:N
    @constraint(model, sum(H[j, box] for box in boxes) == 1)
    p[(:x,j)] = @variable(model, variable_type=Poly(monos))
    p[(:y,j)] = @variable(model, variable_type=Poly(monos))
    for box in boxes
        xl, xu, yl, yu = box.xl, box.xu, box.yl, box.yu
        @constraint(model, p[(:x,j)] >= Mxl + (xl-Mxl)*H[j,box], domain = (t >= T[j] && t <= T[j+1]))
        @constraint(model, p[(:x,j)] <= Mxu + (xu-Mxu)*H[j,box], domain = (t >= T[j] && t <= T[j+1]))
        @constraint(model, p[(:y,j)] >= Myl + (yl-Myl)*H[j,box], domain = (t >= T[j] && t <= T[j+1]))
        @constraint(model, p[(:y,j)] <= Myu + (yu-Myu)*H[j,box], domain = (t >= T[j] && t <= T[j+1]))
    end
end
