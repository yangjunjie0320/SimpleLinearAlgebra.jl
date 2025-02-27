# SimpleLinearAlgebra

[![Build Status](https://github.com/yangjunjie0320/SimpleLinearAlgebra.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yangjunjie0320/SimpleLinearAlgebra.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yangjunjie0320/SimpleLinearAlgebra.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yangjunjie0320/SimpleLinearAlgebra.jl)

This is a simple linear algebra library that builds upon `julia` language.

## Installation

```bash
git clone https://github.com/yangjunjie0320/SimpleLinearAlgebra.jl
cd SimpleLinearAlgebra.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

```julia
using LinearAlgebra
using SimpleLinearAlgebra

# forward substitution
l = tril(rand(10, 10))
b = rand(10)
prob = ForwardSubstitution(l, b)
soln = kernel(prob)

x = soln.x
@assert isapprox(l * x, b, atol=1e-10)
```

