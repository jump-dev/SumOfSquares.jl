@testset "Simple DSOS/SDSOS example with $(factory.constructor)" for factory in sdp_factories
    function sdsos_simple(cone)
        @testset "with $cone" begin
            model = SOSModel(factory)
            @variable(model, α)
            @objective(model, Max, α)
            @polyvar x y
            @constraint(model, x^2 + α*x*y + y^2 in SDSOSCone())
            JuMP.optimize!(model)
            @test JuMP.value(α) ≈ 2.0 rtol=1e-5
            @test JuMP.objective_value(model) ≈ 2.0 rtol=1e-5
        end
    end
    sdsos_simple(DSOSCone())
    sdsos_simple(SDSOSCone())
    sdsos_simple(SOSCone())
end
