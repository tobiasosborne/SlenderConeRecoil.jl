using Test
using SlenderConeRecoil

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

    @testset "expand Laurent powers with epsilon-leading base" begin
        y = Sym(:y)

        leading_only = expand_in(pow(mul(ε, x), Num(-2)), ε, 2)
        @test collect_order(leading_only, ε, -2) == pow(x, Num(-2))
        @test collect_order(leading_only, ε, -1) == Num(0)
        @test collect_order(leading_only, ε, 0) == Num(0)
        @test !occursin("1//0", sprint(show, leading_only))

        base = add(mul(ε, x), mul(pow(ε, Num(3)), y))
        expanded = expand_in(pow(base, Num(-2)), ε, 4)
        bindings = Dict(:x => 2.0, :y => 3.0)

        @test eval_sexpr(collect_order(expanded, ε, -2), bindings) ≈ 0.25
        @test eval_sexpr(collect_order(expanded, ε, 0), bindings) ≈ -0.75
        @test eval_sexpr(collect_order(expanded, ε, 2), bindings) ≈ 27 / 16
        @test eval_sexpr(collect_order(expanded, ε, 4), bindings) ≈ -27 / 8
        @test !occursin("1//0", sprint(show, expanded))
    end

    @testset "explicit integer/Laurent grammar with ε-free coefficients" begin
        y = Sym(:y)
        sqrtx = pow(x, Num(1//2))
        sinx = Func(:sin, SExpr[x])
        expr = add(mul(pow(ε, Num(-1)), x),
                   sinx,
                   mul(ε, sqrtx),
                   mul(pow(ε, Num(2)), y))
        expanded = expand_in(expr, ε, 2)

        @test collect_order(expanded, ε, -1) == x
        @test collect_order(expanded, ε, 0) == sinx
        @test collect_order(expanded, ε, 1) == sqrtx
        @test collect_order(expanded, ε, 2) == y
    end

    @testset "unsupported half-integer and multiple-scale grammar fails clearly" begin
        ξ = Sym(:ξ)
        unsupported = (
            pow(ε, Num(1//2)),
            pow(add(Num(1), ε), Num(1//2)),
            Func(:sin, SExpr[mul(ξ, pow(ε, Num(-1//2)))]),
        )

        for expr in unsupported
            err = try
                expand_in(expr, ε, 3)
                nothing
            catch caught
                caught
            end
            @test err isa ArgumentError
            @test occursin("unsupported ε-expansion grammar",
                           sprint(showerror, err))
        end

        @test_throws ArgumentError expand_in(Num(1), ε, -1)
    end

    @testset "negative powers of zero fail clearly" begin
        @test_throws ArgumentError pow(Num(0), Num(-1))
    end

    @testset "expand constant" begin
        expanded = expand_in(Num(5), ε, 3)
        @test expanded == Num(5)
        @test collect_order(expanded, ε, 0) == Num(5)
        @test collect_order(expanded, ε, 1) == Num(0)
    end
end
