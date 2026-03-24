# Differentiation, pretty-printing, equality/hashing for SExpr types.
# Split from expr.jl to keep files under 200 LOC.

export differentiate

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
        if e.exp isa Num
            n = e.exp
            mul(n, pow(e.base, add(n, Num(-1))), differentiate(e.base, x))
        else
            error("differentiate: non-constant exponent not yet supported")
        end
    elseif e isa Func
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
