# TensorGR.jl Design Patterns for Reuse

Summary from studying `../TensorGR.jl/src/`.

## Expression type hierarchy

Five concrete types under `TensorExpr`:

| Type | Role | Fields |
|------|------|--------|
| `Tensor` | named symbol | name, indices |
| `TScalar` | scalar value | value (Any) |
| `TProduct` | product node | Rational coeff + Vector of factors |
| `TSum` | sum node | Vector of terms |
| `TDeriv` | derivative | index + child + identity |

## Key patterns to reuse

1. **Rational coefficient in product node** — eliminates floating-point drift in symbolic work. `TProduct` stores `coeff::Rational{Int}` times a list of factors.

2. **Smart constructors** — `tproduct()` and `tsum()` normalize on construction: flatten nested sums/products, absorb scalars into coefficients, propagate zeros. This prevents an entire class of bugs.

3. **Bottom-up `walk(f, expr)`** — single tree-walker that rebuilds the tree. Substitution is a one-liner:
   ```julia
   substitute(expr, old, new) = walk(e -> e == old ? new : e, expr)
   ```

4. **Fixed-point simplification** — apply simplification passes in a loop, checking `hash(expr)` for convergence. Each pass is a pure function.

5. **Perturbation expansion** — uses `all_compositions(n, k)` (integer partitions) to distribute perturbation order among factors via Leibniz rule. Memoized.

## What to simplify

- **Drop** index machinery (TIndex, Up/Down, canonicalization, dummy tracking) — not needed for scalar PDEs
- **Replace** `TScalar` wrapping `Any` with typed `Sym` (symbol) and `Num` (number) leaves
- **Drop** rewrite-rule engine with pattern unification — direct `walk` + `substitute` suffices
- **Drop** dummy-index avoidance in perturbation — scalar theory has no indices

## Files to reference

Core patterns live in ~600 LOC across 4 files:
- `types.jl` (114 LOC) — type definitions
- `ast/walk.jl` (125 LOC) — tree walking
- `algebra/arithmetic.jl` (98 LOC) — smart constructors
- `algebra/simplify.jl` (~100 LOC of pipeline logic) — simplification loop
