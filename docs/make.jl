# see documentation at https://juliadocs.github.io/Documenter.jl/stable/
push!(LOAD_PATH,"../src/")

using Documenter, DocumenterCitations, Plots, OpticalFibers
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

ENV["GKSwstype"] = "100"

makedocs(
    #bib,
    modules = [OpticalFibers],
    format = Documenter.HTML(assets=String["assets/citations.css"],),
    authors = "Olivier Vanvincq",
    sitename = "Opticalfibers.jl",
    pages = ["Home"=>"index.md","PhysicalData"=>"PhysicalData.md","ModeSolvers-Tutorial"=>"ModeSolvers-Tutorial.md","ModeSolvers-Solvers"=>"ModeSolvers-Solvers.md","ModeSolvers-Modes and Fields"=>"ModeSolvers.md","Bibliography"=>"Bibliography.md"],
    plugins=[bib],
    #remotes = nothing,
    #checkdocs=:none,
    warnonly=true,
)
println("Finished makedocs")
# Some setup is needed for documentation deployment, see “Hosting Documentation” and
# deploydocs() in the Documenter manual for more information.
deploydocs(
    repo = "https://github.com/ovanvincq/OpticalFibers.jl.git",
    #repo = Remotes.GitHub("OVanvincq","OpticalFibers.jl"),
    push_preview = true
)
