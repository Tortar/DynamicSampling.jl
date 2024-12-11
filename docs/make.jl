
using Documenter
using DynamicSampling

println("Documentation Build")
makedocs(
    modules = [DynamicSampling],
    sitename = "DynamicSampling.jl",
    pages = [
        "API" => "index.md",
    ],
    warnonly = [:doctest, :missing_docs, :cross_references],
)

@info "Deploying Documentation"
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/DynamicSampling.jl.git",
        target = "build",
        push_preview = true,
        devbranch = "main",
    )
end
println("Finished boulding and deploying docs.")
