using Test
include(joinpath(@__DIR__, "..", "src", "expr.jl"))

@testset "Expression types" begin
    r = Sym(:r)
    z = Sym(:z)

    @testset "Leaf construction" begin
        @test r isa Sym
        @test Num(3) isa Num
        @test Num(3).val == 3//1
        @test string(r) == "r"
        @test string(Num(3)) == "3"
        @test string(Num(1//2)) == "1//2"
    end

    @testset "Smart constructors" begin
        # add flattens and drops zeros
        @test add(r, Num(0)) == r
        @test add(Num(0), Num(0)) == Num(0)
        @test add(r, add(z, r)) isa Add
        @test length(add(r, add(z, r)).terms) == 3  # flattened

        # mul absorbs numerics and drops zeros
        @test mul(Num(0), r) == Num(0)
        @test mul(Num(1), r) == r
        @test mul(Num(2), Num(3)) == Num(6)
        @test mul(Num(2), r) isa Mul
        @test mul(Num(2), r).coeff == 2//1

        # pow simplifies
        @test pow(r, Num(0)) == Num(1)
        @test pow(r, Num(1)) == r
        @test pow(Num(2), Num(3)) == Num(8)
    end

    @testset "Operator overloads" begin
        e = r + z
        @test e isa Add
        e2 = 2 * r
        @test e2 isa Mul && e2.coeff == 2//1
        e3 = r^2
        @test e3 isa Pow
        e4 = r / z
        @test e4 isa Mul  # r * z^(-1)
    end

    @testset "Substitution" begin
        e = r + z
        result = substitute(e, r, Num(1))
        @test result isa Add
        # After substitution: Num(1) + z
        @test result.terms[1] == Num(1)
        @test result.terms[2] == z

        # Multi-pair substitution
        e2 = r * z
        result2 = substitute(e2, r => Num(2), z => Num(3))
        @test result2 == Num(6)
    end

    @testset "Differentiate r^2 w.r.t. r gives 2r" begin
        # r^2
        e = pow(r, Num(2))
        de = differentiate(e, r)
        # Should be 2*r^1 * 1 = 2*r
        @test de == mul(Num(2), r)
    end

    @testset "Differentiation rules" begin
        # d/dr(r) = 1
        @test differentiate(r, r) == Num(1)
        # d/dr(z) = 0
        @test differentiate(z, r) == Num(0)
        # d/dr(5) = 0
        @test differentiate(Num(5), r) == Num(0)

        # d/dr(r + z) = 1 + 0 = 1
        @test differentiate(r + z, r) == Num(1)

        # d/dr(r*z) = z  (product rule: 1*z + r*0 = z)
        @test differentiate(r * z, r) == z

        # d/dr(sin(r)) = cos(r)
        e_sin = Func(:sin, [r])
        de_sin = differentiate(e_sin, r)
        @test de_sin == Func(:cos, [r])

        # d/dr(r^3) = 3*r^2
        e3 = r^3
        de3 = differentiate(e3, r)
        @test de3 == mul(Num(3), pow(r, Num(2)))
    end

    @testset "Pretty-print" begin
        e = add(mul(Num(2), r), neg(z))
        s = string(e)
        @test occursin("2", s)
        @test occursin("r", s)
    end

    @testset "Equality and hashing" begin
        @test Sym(:r) == Sym(:r)
        @test Sym(:r) != Sym(:z)
        @test hash(Sym(:r)) == hash(Sym(:r))
        @test Num(3) == Num(3//1)

        d = Dict{SExpr,Int}()
        d[Sym(:r)] = 1
        @test d[Sym(:r)] == 1
    end
end
