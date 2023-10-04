using Plots, DelimitedFiles
VelTime = readdlm("Veltime.txt")
v = VelTime[:,2]
t = VelTime[:,1]
t = t *1e12
avg = readdlm("average.txt")
Autocor=readdlm("autocorr.txt")
PSD = readdlm("PSD.txt")
PSD4 = readdlm("PSD_Vel4.txt")
PSD6 = readdlm("PSD_Vel6.txt")
PSDFFT = readdlm("PSD_VelFFT.txt")
Nsample = 500
plot(t[1:501],v[1:501],label="Velocity",xlabel="Time(ps)",ylabel="Velocity(m/s)")
hline!(avg, label="Average")
savefig("velol20.png")
plot(Autocor[:,1]*1e12,Autocor[:,2],label="",title="Autocorrelation of Velocity",xlabel="Time(ps)", ylabel="Autocorr(10^14 cm^2/s^2)")
savefig("autocor20.png")
plot(PSD[:,1],PSD[:,2]*1e4,yaxis=:log,title="Power Spectrum Density",label="From Autocorrelation",xlabel="Frequency",ylabel="PSD (cm^2/s)")
savefig("PSD.png")
plot!(PSD4[:,1],PSD4[:,2]*1e4,label="N=10⁴")
plot!(PSD6[:,1],PSD6[:,2]*1e4,label="N=10⁶")
savefig("PSD_Vels.png")

plot(PSD[:,1],PSD[:,2]*1e4,yaxis=:log,title="Power Spectrum Density",label="From Autocorrelation",xlabel="Frequency",ylabel="PSD (cm^2/s)")
plot!(PSDFFT[:,1],PSDFFT[:,2]*1e4,label="Improved-Fourier")
savefig("PSD_FFT.png")
