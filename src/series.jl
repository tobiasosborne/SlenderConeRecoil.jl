# Series expansion in a small parameter — mirrors TensorGR perturbation pattern.
# expand_in(expr, param, order) — Taylor expand expression tree in param to given order
# collect_order(expr, param, n) — extract coefficient of param^n

export expand_in, collect_order

# ── Expand expression in small parameter to given order ────────────────
"""
    expand_in(expr, ε, order)

Symbolically expand `expr` as a power series in `ε` up to `O(ε^order)`.
Returns an Add of terms, each a coefficient times ε^k.
"""
function expand_in(expr::SExpr, ε::Sym, order::Int)::SExpr
    # Strategy: recursively expand each node type, truncating at `order`.
    _expand(expr, ε, order)
end

function _expand(e::Num, ε::Sym, order::Int)::SExpr
    e
end

function _expand(e::Sym, ε::Sym, order::Int)::SExpr
    e
end

function _expand(e::Add, ε::Sym, order::Int)::SExpr
    add([_expand(t, ε, order) for t in e.terms]...)
end

function _expand(e::Mul, ε::Sym, order::Int)::SExpr
    # Expand each factor, then multiply and truncate
    expanded_factors = [_expand(f, ε, order) for f in e.factors]
    result = Num(e.coeff)
    for f in expanded_factors
        result = _mul_and_truncate(result, f, ε, order)
    end
    result
end

function _expand(e::Pow, ε::Sym, order::Int)::SExpr
    if e.exp isa Num && isinteger(e.exp.val)
        n = Int(e.exp.val)
        base_expanded = _expand(e.base, ε, order)
        if n == 0
            return Num(1)
        elseif n > 0
            result = base_expanded
            for _ in 2:n
                result = _mul_and_truncate(result, base_expanded, ε, order)
            end
            return result
        else
            # Negative integer power: factor the lowest ε order, then use
            # the binomial expansion around a unit leading term.
            return _expand_negative_pow(base_expanded, -n, ε, order)
        end
    end
    # Non-integer exponent: leave as-is
    pow(_expand(e.base, ε, order), _expand(e.exp, ε, order))
end

function _expand(e::Func, ε::Sym, order::Int)::SExpr
    # For now, leave functions unexpanded unless their argument contains ε
    Func(e.name, [_expand(a, ε, order) for a in e.args])
end

# ── Expand (base)^{-n} as Laurent/binomial series ──────────────────────
"""
Expand base^{-n} where the base may have a nonzero lowest ε order:
base = ε^m a_m * (1 + r). Then
base^{-n} = ε^{-mn} a_m^{-n} * (1 + r)^{-n}.
"""
function _expand_negative_pow(base::SExpr, n::Int, ε::Sym, order::Int)::SExpr
    m, a_m = _lowest_nonzero_order(base, ε)
    leading_order = -n * m
    leading_order > order && return Num(0)

    unit_order = order - leading_order
    r = _relative_tail(base, ε, m, a_m, unit_order)
    series = _negative_binomial_series(r, n, ε, unit_order)
    coeff_series = _mul_and_truncate(pow(a_m, Num(-n)), series, ε, unit_order)

    leading_order == 0 &&
        return _mul_and_truncate(Num(1), coeff_series, ε, order)

    _mul_and_truncate(pow(ε, Num(leading_order)), coeff_series, ε, order)
end

function _lowest_nonzero_order(e::SExpr, ε::Sym)::Tuple{Int,SExpr}
    orders = sort!(unique(k for (k, _) in _to_order_terms(e, ε)))
    for k in orders
        coeff = collect_order(e, ε, k)
        !_iszero_expr(coeff) && return (k, coeff)
    end
    throw(ArgumentError("cannot expand negative power in $(ε.name): base has no nonzero ε-series terms"))
end

function _relative_tail(base::SExpr, ε::Sym, leading_order::Int,
                        leading_coeff::SExpr, order::Int)::SExpr
    terms = SExpr[]
    orders = sort!(unique(k for (k, _) in _to_order_terms(base, ε)))
    for k in orders
        k <= leading_order && continue
        relative_order = k - leading_order
        relative_order > order && continue

        coeff = collect_order(base, ε, k)
        _iszero_expr(coeff) && continue

        scaled_coeff = _mul_and_truncate(coeff, pow(leading_coeff, Num(-1)), ε, 0)
        term = relative_order == 0 ?
            scaled_coeff :
            _mul_and_truncate(scaled_coeff, pow(ε, Num(relative_order)), ε, order)
        push!(terms, term)
    end

    isempty(terms) ? Num(0) : add(terms...)
end

function _negative_binomial_series(x::SExpr, n::Int, ε::Sym, order::Int)::SExpr
    _iszero_expr(x) && return Num(1)

    x_expanded = _expand(x, ε, order)
    series = Num(1)
    x_power = Num(1)  # x^0
    binom_coeff = Rational{Int}(1)

    for k in 1:order
        binom_coeff *= (-n - k + 1) // k
        x_power = _mul_and_truncate(x_power, x_expanded, ε, order)
        _iszero_expr(x_power) && break
        term = _mul_and_truncate(Num(binom_coeff), x_power, ε, order)
        series = add(series, term)
    end

    series
end

function _iszero_expr(e::SExpr)::Bool
    e isa Num && iszero(e.val)
end

# ── Multiply two expanded expressions, truncating at order ─────────────
"""
Multiply two expressions that are already expanded in ε, dropping terms
above `order`.
"""
function _mul_and_truncate(a::SExpr, b::SExpr, ε::Sym, order::Int)::SExpr
    a_terms = _to_order_terms(a, ε)
    b_terms = _to_order_terms(b, ε)

    result_terms = SExpr[]
    for (ka, ca) in a_terms
        for (kb, cb) in b_terms
            k = ka + kb
            k > order && continue
            term = mul(ca, cb)
            if !(term isa Num && iszero(term.val))
                if k == 0
                    push!(result_terms, term)
                else
                    push!(result_terms, mul(term, pow(ε, Num(k))))
                end
            end
        end
    end
    isempty(result_terms) ? Num(0) : add(result_terms...)
end

# ── Decompose expression into (order, coefficient) pairs ───────────────
"""
Break an expression into pairs (k, coeff) where the term is coeff * ε^k.
"""
function _to_order_terms(e::SExpr, ε::Sym)::Vector{Tuple{Int,SExpr}}
    if e isa Add
        result = Tuple{Int,SExpr}[]
        for t in e.terms
            append!(result, _to_order_terms(t, ε))
        end
        return result
    end
    k, c = _extract_eps_order(e, ε)
    [(k, c)]
end

"""
Given a term, extract (k, coeff) such that term = coeff * ε^k.
"""
function _extract_eps_order(e::SExpr, ε::Sym)::Tuple{Int,SExpr}
    # e = ε^k  →  (k, 1)
    if e == ε
        return (1, Num(1))
    end
    if e isa Pow && e.base == ε && e.exp isa Num && isinteger(e.exp.val)
        return (Int(e.exp.val), Num(1))
    end
    # e = c * f1 * f2 * ...  →  collect ε powers from factors
    if e isa Mul
        total_k = 0
        remaining = SExpr[]
        for f in e.factors
            if f == ε
                total_k += 1
            elseif f isa Pow && f.base == ε && f.exp isa Num && isinteger(f.exp.val)
                total_k += Int(f.exp.val)
            else
                push!(remaining, f)
            end
        end
        coeff = isempty(remaining) ? Num(e.coeff) : mul(Num(e.coeff), remaining...)
        return (total_k, coeff)
    end
    # Constant w.r.t. ε
    (0, e)
end

# ── Collect coefficient of ε^n ─────────────────────────────────────────
"""
    collect_order(expr, ε, n)

Extract the coefficient of ε^n from `expr`.
"""
function collect_order(expr::SExpr, ε::Sym, n::Int)::SExpr
    terms = _to_order_terms(expr, ε)
    coeffs = [c for (k, c) in terms if k == n]
    isempty(coeffs) ? Num(0) : add(coeffs...)
end
