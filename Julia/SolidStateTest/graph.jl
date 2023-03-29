using DelimitedFiles
using Plots
bands = readdlm("bandas.out")
p = plot((bands[:,1],[bands[:,2] bands[:,3] bands[:,4] bands[:,5]]),legend=false)
savefig(p,"bandas.png")
