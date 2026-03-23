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

# ── Differentiation (structural) ──────────────────────────────────────
function differentiate(e::SExpr, x::Sym)::SExpr
    if e isa Num
        Num(0)
    elseif e isa Sym
        e == x ? Num(1) : Num(0)
    elseif e isa Add
        add([differentiate(t, x) for t in e.terms]...)
    elseif e isa Mul
        # Product rule: d/dx (c * f1*f2*...) = c * Σ_i (f1*...*f_i'*...*fn)
        fs = e.factors
        terms = SExpr[]
        for i in eachindex(fs)
            di = differentiate(fs[i], x)
            if !(di isa Num && iszero(di.val))
                others = [j == i ? di : fs[j] for j in eachindex(fs)]
                push!(terms, mul(Num(e.coeff), others...))
            end
        end
        isempty(terms) ? Num(0) : add(terms...)
    elseif e isa Pow
        # d/dx (base^exp): handle exp=const case and general case
        if e.exp isa Num
            # n * base^(n-1) * base'
            n = e.exp
            mul(n, pow(e.base, add(n, Num(-1))), differentiate(e.base, x))
        else
            # General: base^exp * (exp' * ln(base) + exp * base'/base)
            # For now, only handle const exponent — flag others
            error("differentiate: non-constant exponent not yet supported")
        end
    elseif e isa Func
        # Chain rule for known functions
        if length(e.args) == 1
            inner = e.args[1]
            di = differentiate(inner, x)
            outer_deriv = _func_deriv(e.name, inner)
            mul(outer_deriv, di)
        else
            error("differentiate: multi-arg Func not yet supported")
        end
    else
        error("differentiate: unknown expression type")
    end
end

function _func_deriv(name::Symbol, arg::SExpr)::SExpr
    if name == :sin
        Func(:cos, [arg])
    elseif name == :cos
        neg(Func(:sin, [arg]))
    elseif name == :exp
        Func(:exp, [arg])
    elseif name == :log
        pow(arg, Num(-1))
    elseif name == :sqrt
        mul(Num(1//2), pow(arg, Num(-1//2)))
    else
        error("differentiate: unknown function $name")
    end
end

# ── Pretty-print ───────────────────────────────────────────────────────
function Base.show(io::IO, e::Add)
    for (i, t) in enumerate(e.terms)
        if i > 1
            # Check if term has negative coefficient
            if t isa Mul && t.coeff < 0
                print(io, " - ")
                show(io, Mul(-t.coeff, t.factors))
            elseif t isa Num && t.val < 0
                print(io, " - ")
                show(io, Num(-t.val))
            else
                print(io, " + ")
                show(io, t)
            end
        else
            show(io, t)
        end
    end
end

function Base.show(io::IO, e::Mul)
    if e.coeff == -1
        print(io, "-")
    elseif e.coeff != 1
        show(io, Num(e.coeff))
        print(io, "*")
    end
    for (i, f) in enumerate(e.factors)
        i > 1 && print(io, "*")
        needs_parens = f isa Add
        needs_parens && print(io, "(")
        show(io, f)
        needs_parens && print(io, ")")
    end
end

function Base.show(io::IO, e::Pow)
    needs_parens = !(e.base isa Sym || e.base isa Num || e.base isa Func)
    needs_parens && print(io, "(")
    show(io, e.base)
    needs_parens && print(io, ")")
    print(io, "^")
    needs_parens_exp = !(e.exp isa Sym || e.exp isa Num)
    needs_parens_exp && print(io, "(")
    show(io, e.exp)
    needs_parens_exp && print(io, ")")
end

function Base.show(io::IO, e::Func)
    print(io, e.name, "(")
    for (i, a) in enumerate(e.args)
        i > 1 && print(io, ", ")
        show(io, a)
    end
    print(io, ")")
end

# ── Equality and hashing ──────────────────────────────────────────────
Base.:(==)(a::Sym, b::Sym) = a.name == b.name
Base.:(==)(a::Num, b::Num) = a.val == b.val
Base.:(==)(a::Add, b::Add) = a.terms == b.terms
Base.:(==)(a::Mul, b::Mul) = a.coeff == b.coeff && a.factors == b.factors
Base.:(==)(a::Pow, b::Pow) = a.base == b.base && a.exp == b.exp
Base.:(==)(a::Func, b::Func) = a.name == b.name && a.args == b.args

Base.hash(s::Sym, h::UInt) = hash(s.name, h)
Base.hash(n::Num, h::UInt) = hash(n.val, h)
Base.hash(a::Add, h::UInt) = hash(a.terms, hash(:Add, h))
Base.hash(m::Mul, h::UInt) = hash(m.factors, hash(m.coeff, hash(:Mul, h)))
Base.hash(p::Pow, h::UInt) = hash(p.exp, hash(p.base, hash(:Pow, h)))
Base.hash(f::Func, h::UInt) = hash(f.args, hash(f.name, hash(:Func, h)))
