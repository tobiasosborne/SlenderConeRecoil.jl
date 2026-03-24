# Lightweight symbolic expression tree — mirrors TensorGR.jl AST patterns.
# Types: Sym, Num, Add, Mul, Pow, Func
# Operations: substitute, differentiate, pretty-print

export SExpr, Sym, Num, Add, Mul, Pow, Func
export substitute, differentiate, walk, add, mul, pow, neg

# ── Abstract type ──────────────────────────────────────────────────────
abstract type SExpr end

# ── Leaf types ─────────────────────────────────────────────────────────
struct Sym <: SExpr
    name::Symbol
end
Base.show(io::IO, s::Sym) = print(io, s.name)

struct Num <: SExpr
    val::Rational{Int}
end
Num(x::Int) = Num(Rational{Int}(x))
Num(x::Float64) = Num(rationalize(Int, x))
Base.show(io::IO, n::Num) = isone(denominator(n.val)) ? print(io, numerator(n.val)) : print(io, n.val)

# ── Branch types ───────────────────────────────────────────────────────
struct Add <: SExpr
    terms::Vector{SExpr}
end

struct Mul <: SExpr
    coeff::Rational{Int}
    factors::Vector{SExpr}
end

struct Pow <: SExpr
    base::SExpr
    exp::SExpr
end

struct Func <: SExpr
    name::Symbol      # :sin, :cos, :besselj, etc.
    args::Vector{SExpr}
end

# ── Smart constructors (normalize on construction) ─────────────────────
function add(args::SExpr...)
    terms = SExpr[]
    numsum = Rational{Int}(0)
    for a in args
        if a isa Add
            for t in a.terms
                if t isa Num
                    numsum += t.val
                else
                    push!(terms, t)
                end
            end
        elseif a isa Num
            numsum += a.val
        else
            push!(terms, a)
        end
    end
    if !iszero(numsum)
        pushfirst!(terms, Num(numsum))
    end
    isempty(terms) && return Num(0)
    length(terms) == 1 && return terms[1]
    sort!(terms, by=_sort_key)
    Add(terms)
end

function mul(args::SExpr...)
    coeff = Rational{Int}(1)
    factors = SExpr[]
    for a in args
        if a isa Num
            coeff *= a.val
        elseif a isa Mul
            coeff *= a.coeff
            append!(factors, a.factors)
        else
            push!(factors, a)
        end
    end
    iszero(coeff) && return Num(0)
    isempty(factors) && return Num(coeff)
    # If all remaining factors are Num, fold them in
    if all(f -> f isa Num, factors)
        for f in factors
            coeff *= f.val
        end
        return Num(coeff)
    end
    isone(coeff) && length(factors) == 1 && return factors[1]
    sort!(factors, by=_sort_key)
    Mul(coeff, factors)
end

# Canonical ordering key for deterministic equality
_sort_key(e::Sym) = (0, string(e.name), 0)
_sort_key(e::Num) = (1, "", Float64(e.val))
_sort_key(e::Pow) = (2, string(_sort_key(e.base)), 0)
_sort_key(e::Func) = (3, string(e.name), 0)
_sort_key(e::Add) = (4, "", length(e.terms))
_sort_key(e::Mul) = (5, "", Float64(e.coeff))

function pow(base::SExpr, exp::SExpr)
    exp isa Num && iszero(exp.val) && return Num(1)
    exp isa Num && isone(exp.val) && return base
    base isa Num && exp isa Num && isinteger(exp.val) &&
        return Num(base.val ^ Int(exp.val))
    Pow(base, exp)
end

function neg(a::SExpr)
    mul(Num(-1), a)
end

# ── Operator overloads ─────────────────────────────────────────────────
Base.:+(a::SExpr, b::SExpr) = add(a, b)
Base.:-(a::SExpr, b::SExpr) = add(a, neg(b))
Base.:-(a::SExpr) = neg(a)
Base.:*(a::SExpr, b::SExpr) = mul(a, b)
Base.:*(c::Integer, a::SExpr) = mul(Num(c), a)
Base.:*(a::SExpr, c::Integer) = mul(Num(c), a)
Base.:/(a::SExpr, b::SExpr) = mul(a, pow(b, Num(-1)))
Base.:^(a::SExpr, b::SExpr) = pow(a, b)
Base.:^(a::SExpr, n::Integer) = pow(a, Num(n))

# ── Tree walking (bottom-up) ──────────────────────────────────────────
function children(e::SExpr)
    e isa Sym && return SExpr[]
    e isa Num && return SExpr[]
    e isa Add && return e.terms
    e isa Mul && return e.factors
    e isa Pow && return SExpr[e.base, e.exp]
    e isa Func && return e.args
    SExpr[]
end

function walk(f, e::SExpr)
    if e isa Sym || e isa Num
        return f(e)
    elseif e isa Add
        f(add([walk(f, t) for t in e.terms]...))
    elseif e isa Mul
        f(mul(Num(e.coeff), [walk(f, t) for t in e.factors]...))
    elseif e isa Pow
        f(pow(walk(f, e.base), walk(f, e.exp)))
    elseif e isa Func
        f(Func(e.name, [walk(f, a) for a in e.args]))
    else
        error("walk: unhandled expression type $(typeof(e))")
    end
end

# ── Substitution ───────────────────────────────────────────────────────
function substitute(expr::SExpr, old::SExpr, new::SExpr)
    walk(e -> e == old ? new : e, expr)
end

function substitute(expr::SExpr, pairs::Pair{<:SExpr,<:SExpr}...)
    result = expr
    for (old, new) in pairs
        result = substitute(result, old, new)
    end
    result
end

include("expr_ops.jl")
