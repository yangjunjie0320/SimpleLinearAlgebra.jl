# setup a new julia package

using PkgTemplates

p = [License(; name="MIT"), Git(; ssh=true)]
append!(p, [GitHubActions(; x86=true)])
append!(p, [Codecov()])
# append!(p, [Documenter{GitHubActions}()])

println(p)

t = Template(;
    user = "yangjunjie0320",
    name = "SimpleLinearAlgebra",
    julia = v"1.0.0",
    plugins = p,
)

t("SimpleLinearAlgebra")