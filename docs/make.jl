using Documenter, Elfel

makedocs(
	modules = [Elfel, Elfel.RefShapes, Elfel.FElements, Elfel.FESpaces, Elfel.FEIterators, Elfel.QPIterators],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "Elfel.jl",
	pages = Any[
			"Home" => "index.md",
			"How to guide" => "guide/guide.md",
			"Reference" => Any[
				"man/types.md",
				"man/functions.md"],
			"Concepts" => "concepts/concepts.md"	
		],
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/Elfel.jl.git",
)
