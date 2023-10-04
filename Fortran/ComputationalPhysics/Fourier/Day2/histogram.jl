using Plots, DelimitedFiles
file = "serierayleigh"
Data = readdlm(file*".txt")[:,3]
histogram(Data,normalize=:pdf)
savefig(file*".png")
