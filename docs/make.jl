using Documenter, ProjectedGPE

makedocs(
    modules = [ProjectedGPE],
    doctest=true
)

deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/AshtonSBradley/ProjectedGPE.jl.git",
    julia  = "0.5.2",
    osname = "macos")


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
