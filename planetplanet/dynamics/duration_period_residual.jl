using PyPlot
using CGS
nplanet = 7
t = [36.4,42.37,49.13,57.21,62.6,68.4,75.6]
#b = [0.126,0.161,0.17,0.12, 0.382,0.421,0.45]
sigt = [0.17,0.22,0.65,0.71,0.6,0.66,1.15]
rr_obs = [0.08524846380071895,0.08285759783610132, 0.06054351753969651, 0.07199543643724797, 0.08204346713215337, 0.08841765826354074, 0.059324407529117607]

p = [1.51087081,2.4218233,4.0949610,6.099615,9.20669,12.35294,18.7659]
b0 = 0.1
b = b0.*(p./p[1]).^(2./3.)
bobs = [0.126,0.161,0.17,0.12,0.382,0.421,0.45]
#t = t./sqrt(1.-b.^2)
t = t./sqrt((1.+rr_obs).^2-bobs.^2)

clf()
#scatter(log10(p),log10(t))
#scatter(p,t)
p0 = linspace(0,1.3,100)

rho_sun= MSUN/(4pi/3*RSUN^3)
fn = zeros(1,nplanet)
fn[1,:] = (3/pi^2/GRAV/rho_sun*p*24.*3600.).^(1./3.)/60. 
coeff,cov = regress(fn,t,sigt)
rho = rho_sun/coeff[1]^3
println("Density: ",1./coeff[1]^3)
dur=(3/pi^2/GRAV/rho*10.^p0*24.*3600.).^(1./3.)/60.
#plot(10.^p0,dur,linewidth=3,alpha=0.5)
#errorbar(p,t,yerr = sigt, fmt = "o",linewidth=3,alpha=0.5)
#xlabel("Period [d]",fontsize=20,fontweight="bold")
#ylabel("Rescaled transit duration [min]",fontsize=20,fontweight="bold")
#axis([0,20,30,80])
#tick_params("both",labelsize=20)


subplot2grid((4,1), (0,0),rowspan = 3)
subplots_adjust(hspace=0.0) # Set the vertical spacing between axes
plot(10.^p0,dur,linewidth=3,alpha=0.5)
errorbar(p,t,yerr = sigt, fmt = "o",linewidth=3,alpha=0.5)
ylabel("Rescaled transit duration [min]",fontsize=20,fontweight="bold")
axis([0,20,30,80])
#setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
#plot(rand(10))
subplot2grid((4,1), (3,0), rowspan=1)
dur_mod=(3/pi^2/GRAV/rho*p*24.*3600.).^(1./3.)/60.
errorbar(p,t-dur_mod,yerr = sigt, fmt = "o",linewidth=3,alpha=0.5,color="g")
plot([0,20],[0,0],linewidth=3,alpha=0.5)
ylabel("Residual [min]",fontsize=20,fontweight="bold")
xlabel("Period [d]",fontsize=20,fontweight="bold")
#fig[:canvas][:draw]() # Update the figure

#fig = figure("pyplot_subplot_touching",figsize=(10,10))
#subplots_adjust(vspace=0.0) # Set the vertical spacing between axes
#subplot(211) # Create the 1st axis of a 2x1 array of axes
#ax1 = gca()
##ax1[:set_xscale]("log") # Set the x axis to a logarithmic scale
#setp(ax1[:get_xticklabels](),visible=false) # Disable x tick labels
#grid("on")
#title("Title")
#yticks(0.1:0.2:0.9) # Set the y-tick range and step size, 0.1 to 0.9 in increments of 0.2
#ylim(0.0,1.0) # Set the y-limits from 0.0 to 1.0
#subplot(212,sharex=ax1) # Create the 2nd axis of a 2x1 array of axes
#ax2 = gca()
##ax2[:set_xscale]("log") # Set the x axis to a logarithmic scale
#setp(ax2[:get_xticklabels](),visible=false) # Disable x tick labels
#grid("on")
##ylabel("Log Scale")
#yticks(0.1:0.2:0.9)
#ylim(0.0,1.0)
#fig[:canvas][:draw]() # Update the figure
#suptitle("2x1 Shared Axis")

tight_layout()
savefig("duration_period.pdf", bbox_inches="tight")
