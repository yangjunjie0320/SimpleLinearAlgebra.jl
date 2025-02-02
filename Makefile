jl:=julia

init: src/*jl
	$(jl) -e 'd=pwd(); @assert isdir(d);\
	using Pkg: activate, instantiate, develop, precompile; \
	activate(d); instantiate(); activate(joinpath(d, "examples")); \
	develop(path=d); instantiate(); precompile();'
	@echo "Environment initialized at: $$PWD"

update: src/*jl
	$(jl) -e 'd=pwd(); @assert isdir(d);\
	using Pkg: activate, update, develop, precompile; \
	activate(d); update(); activate(joinpath(d, "examples")); \
	update(); precompile();'
	@echo "Environment updated at: $$PWD"

test: test/*jl
	@echo "Running tests at: $$PWD"
	$(jl) -e 'd=pwd(); @assert isdir(d);\
	using Pkg: activate, test; \
	activate(d); test();'
