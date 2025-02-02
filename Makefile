jl:=julia

init: src/* examples/*
	$(jl) --project=${PWD} -e 'using Pkg; Pkg.instantiate(); Pkg.precompile();'
	$(jl) --project=${PWD}/examples -e 'using Pkg; Pkg.develop(path="${PWD}"); Pkg.instantiate(); Pkg.precompile();'
	echo "Environment initialized at: ${PWD}"

update: src/* examples/*
	$(jl) --project=${PWD} -e 'using Pkg; Pkg.update()'
	$(jl) --project=${PWD}/examples -e 'using Pkg; Pkg.update(); Pkg.precompile()'
	echo "Environment updated at: ${PWD}"

test: src/* test/*
	echo "Running tests at: ${PWD}"
	$(jl) --project=${PWD} -e 'using Pkg; Pkg.test();'

# examples: src/* examples/*
# 	echo "Running examples at: ${PWD}"
# 	$(jl) --project=${PWD}/examples -e 'using Pkg; include(joinpath(pwd(), "examples", "main.jl"));'
