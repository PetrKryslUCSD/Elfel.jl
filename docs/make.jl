using Documenter, Elfel

makedocs(
	modules = [Elfel],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "Elfel.jl",
	pages = Any[
			"Home" => "index.md",
			"How to guide" => "guide/guide.md",
			"Reference" => "man/reference.md",
			"Concepts" => "concepts/concepts.md"	
		],
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/Elfel.jl.git",
)
