using SomaticPhylogenies
using Documenter

DocMeta.setdocmeta!(SomaticPhylogenies, :DocTestSetup, :(using SomaticPhylogenies); recursive=true)

makedocs(;
    modules=[SomaticPhylogenies],
    authors="Matthias Günther <ma.guenther@dkfz.de>",
    sitename="SomaticPhylogenies.jl",
    format = Documenter.HTML(
        prettyurls = false
    ),
    pages=[
        "Home" => "index.md",
        "References" => "ref.md"
    ],
)

deploydocs(;
    repo="github.com/MG-dkfz/SomaticPhylogenies.jl.git",
    devbranch="main",
)
