# see documentation at https://juliadocs.github.io/Documenter.jl/stable/
push!(LOAD_PATH,"../src/")

using Documenter, DocumenterCitations, OpticalFibers
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

ENV["GKSwstype"] = "100"

makedocs(
    #bib,
    modules = [OpticalFibers],
    format = Documenter.HTML(assets=String["assets/citations.css"],),
    authors = "Olivier Vanvincq",
    sitename = "OpticalFibers.jl",
    pages = ["Home"=>"index.md","Common functions"=>"Common.md","PhysicalData"=>"PhysicalData.md","ModeSolvers-Tutorial"=>"ModeSolvers-Tutorial.md","ModeSolvers-Solvers"=>"ModeSolvers-Solvers.md","ModeSolvers-Modes and Fields"=>"ModeSolvers.md","Bibliography"=>"Bibliography.md"],
    plugins=[bib],
    #remotes = nothing,
    #checkdocs=:none,
    warnonly=true,
)
println("Finished makedocs")
# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "github.com/ovanvincq/OpticalFibers.jl.git",
    push_preview = true,
    devbranch = "main",
)
