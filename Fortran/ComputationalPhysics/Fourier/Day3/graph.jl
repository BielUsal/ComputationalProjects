using Plots, DelimitedFiles
name = "heavi10_"
Time = readdlm(name*"time_series.txt")
Fourier = readdlm(name*"fourier.txt")
plot(Fourier[:,1], Fourier[:,2],label="Real")
plot!(Fourier[:,1], Fourier[:,3], label="Imaginary")
savefig(name*"fourier.png")
plot(Fourier[:,1], Fourier[:,4], label="Absolute")
savefig(name*"Freqs.png")
plot(Time[:,1], Time[:,2])
savefig(name*"timeseries.png")
