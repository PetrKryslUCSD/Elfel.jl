using Literate

for t in readdir(".")
    if occursin(r"tut_.*.jl", t)
        println("\nTutorial $t in $(pwd())\n")
        Literate.markdown(t, "."; documenter=false);
    end
end
