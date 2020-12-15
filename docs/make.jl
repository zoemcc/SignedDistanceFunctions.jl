using SignedDistanceFunctions
using Documenter

makedocs(;
    modules=[SignedDistanceFunctions],
    authors="Zoe McCarthy <zoemccarthy12@gmail.com> and contributors",
    repo="https://github.com/zoemcc/SignedDistanceFunctions.jl/blob/{commit}{path}#L{line}",
    sitename="SignedDistanceFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://zoemcc.github.io/SignedDistanceFunctions.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zoemcc/SignedDistanceFunctions.jl",
)
