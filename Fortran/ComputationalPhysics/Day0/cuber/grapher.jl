using Plots, DelimitedFiles
data = readdlm(cubes.dat)
plot(data[1],data[2],data[3])
