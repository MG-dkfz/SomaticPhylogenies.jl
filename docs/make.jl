using SomaticPhylogenies
using Documenter

DocMeta.setdocmeta!(SomaticPhylogenies, :DocTestSetup, :(using SomaticPhylogenies); recursive=true)

makedocs(;
    modules=[SomaticPhylogenies],
    authors="Matthias Günther <ma.guenther@dkfz.de>",
    sitename="SomaticPhylogenies.jl",
    format=Documenter.HTML(;
        canonical="https://mg-dkfz.github.io/SomaticPhylogenies.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "References" => "ref.md"
    ],
)

deploydocs(;
    repo="github.com/mg-dkfz/SomaticPhylogenies.jl",
    devbranch="main",
)
