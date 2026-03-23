using Test
include(joinpath(@__DIR__, "..", "src", "expr.jl"))
include(joinpath(@__DIR__, "..", "src", "series.jl"))

@testset "Series expansion" begin
    ε = Sym(:ε)
    x = Sym(:x)

    @testset "collect_order basics" begin
        # 3 + 2ε + ε²
        e = add(Num(3), mul(Num(2), ε), pow(ε, Num(2)))
        @test collect_order(e, ε, 0) == Num(3)
        @test collect_order(e, ε, 1) == Num(2)
        @test collect_order(e, ε, 2) == Num(1)
        @test collect_order(e, ε, 3) == Num(0)
    end

    @testset "expand (1 + εx)^{-1} to O(ε³)" begin
        # (1 + εx)^{-1} = 1 - εx + ε²x² - ε³x³ + O(ε⁴)
        base = add(Num(1), mul(ε, x))
        expr = pow(base, Num(-1))
        expanded = expand_in(expr, ε, 3)

        c0 = collect_order(expanded, ε, 0)
        c1 = collect_order(expanded, ε, 1)
        c2 = collect_order(expanded, ε, 2)
        c3 = collect_order(expanded, ε, 3)

        @test c0 == Num(1)
        @test c1 == neg(x)
        @test c2 == mul(x, x)
        @test c3 == neg(mul(x, mul(x, x)))
    end

    @testset "expand (1 + ε)^2 to O(ε²)" begin
        base = add(Num(1), ε)
        expr = pow(base, Num(2))
        expanded = expand_in(expr, ε, 2)

        @test collect_order(expanded, ε, 0) == Num(1)
        @test collect_order(expanded, ε, 1) == Num(2)
        @test collect_order(expanded, ε, 2) == Num(1)
    end

    @testset "expand product with truncation" begin
        # (1 + ε)(1 + ε) = 1 + 2ε + ε²
        a = add(Num(1), ε)
        expanded = expand_in(mul(a, a), ε, 1)
        # At O(ε): should have 1 + 2ε, dropping ε² term
        @test collect_order(expanded, ε, 0) == Num(1)
        @test collect_order(expanded, ε, 1) == Num(2)
        @test collect_order(expanded, ε, 2) == Num(0)  # truncated
    end

    @testset "expand (1 + εx)^{-2} to O(ε²)" begin
        base = add(Num(1), mul(ε, x))
        expr = pow(base, Num(-2))
        expanded = expand_in(expr, ε, 2)

        # (1+εx)^{-2} = 1 - 2εx + 3ε²x² + ...
        @test collect_order(expanded, ε, 0) == Num(1)
        @test collect_order(expanded, ε, 1) == mul(Num(-2), x)
        @test collect_order(expanded, ε, 2) == mul(Num(3), mul(x, x))
    end

    @testset "expand constant" begin
        expanded = expand_in(Num(5), ε, 3)
        @test expanded == Num(5)
        @test collect_order(expanded, ε, 0) == Num(5)
        @test collect_order(expanded, ε, 1) == Num(0)
    end
end
