# Simulate two eccentric planets, and time offsets of planet-planet occultation.
#include("/Users/ericagol/Courses/PlanetaryDynamics/ExoJulia/ExoJulia/exojulia.jl")
using CGS
using PyPlot
include("regress.jl")

# All we need is the X-coordinate of the orbit:
# Reference direction is sky plane:

function calc_position(time,gmstar,elements)
# We'll use t0, time of transit; need to convert to time of periastron:
  period = elements[2]
  psec = period*24.*3600.
  semi = ((psec/2/pi)^2*gmstar)^(1./3.)  # semi-major axis in cm if gmstar is cgs
  t0 = elements[3] # Time of transit
  ecc = elements[4]^2+elements[5]^2  # sqrt(e)*cos(omega),sqrt(e)*sin(omega) (omega is omega of planet, not star)
  omega = atan2(elements[5],elements[4])
  theta =  -.5*pi-omega  # equation 7.53 from Marion & Thornton
  tp = t0 - period/(2.*pi)*(2.*atan(sqrt((1.-ecc)/(1.+ecc))*tan(.5*theta))-ecc*sqrt(1.-ecc^2)*sin(theta)/(1.+ecc*cos(theta)))
  m = (time-tp).*2pi./period
  ntime = length(time)
  eanom = zeros(ntime)
  for j=1:ntime
    eanom[j]=ExoJulia.Orbit.kepler_solve(m[j],ecc)
  end
  # equation 6.3.23 from Danby:
  X = semi.*(cos(eanom)-ecc)
  Y = semi.*sqrt(1.-ecc^2).*sin(eanom)
  R = sqrt(X.^2+Y.^2)
  # equation 6.3.24 from Danby:
  n = 2pi/psec
  dXdt = -n*semi^2./R.*sin(eanom)
  dYdt = -n*semi^2*sqrt(1.-ecc^2)./R.*cos(eanom)
  # Now, rotate to observer's frame of reference:
  Xobs = X.*cos(omega)-Y.*sin(omega) # when omega=pi/2, Y is in negative X direction
  Yobs = Y.*cos(omega)+X.*sin(omega) # when omega=pi/2, Y is in negative X direction
  dXdtobs = dXdt.*cos(omega)-dYdt.*sin(omega)
return Xobs,Yobs,dXdtobs
end

# For specificity, use Trappist-1b/c parameters:
#eb = abs(randn()*.02)
eb = 0.1
omb = 2pi*0.125
pb= 1.51087081
t0b=322.51736
elements_b = [1e-4,pb,t0b,sqrt(eb)*cos(omb),sqrt(eb)*sin(omb)]
p1sec = pb*24.*3600.
gmstar = 0.0802*MSUN*GRAV
semi1 = ((p1sec/2/pi)^2*gmstar)^(1./3.)  # semi-major axis in cm if gmstar is cgs
#ec = abs(randn()*.02)
ec = 0.1
omc = 2pi*0.375
pc= 2.4218233
t0c = 282.80728
elements_c = [1e-4,pc,t0c,sqrt(ec)*cos(omc),sqrt(ec)*sin(omc)]
p2sec = pc*24.*3600.
semi2 = ((p2sec/2/pi)^2*gmstar)^(1./3.)  # semi-major axis in cm if gmstar is cgs

t1 = 300.0
t2 = 400.0
time = collect(t1:1e-4:t2)
ntime = length(time)

xb,yb,vb = calc_position(time,gmstar,elements_b)
# Now, run with zero eccentricity for comparison:
elements_b0 = copy(elements_b)
elements_b0[4:5]=0.
xb0,yb0,vb0 = calc_position(time,gmstar,elements_b0)
xc,yc,vc = calc_position(time,gmstar,elements_c)
# Now, run with zero eccentricity for comparison:
elements_c0 = copy(elements_c)
elements_c0[4:5]=0.
xc0,yc0,vc0 = calc_position(time,gmstar,elements_c0)

#plot(time,xb)
#plot(time,xc)
#plot(time,xb0)
#plot(time,xc0)

# Check that transits occur at correct times:
for i=-100:100
  t0 = elements_b[3]+i*elements_b[2]
  if(t0 > t1 && t0 < t2)
#    plot([t0,t0],[-4,4])
  end
  t0 = elements_c[3]+i*elements_c[2]
  if(t0 > t1 && t0 < t2)
#    plot([t0,t0],[-4,4])
  end
end

# Now, find the timing offsets of planet-planet occultation events:
dx = xc-xb
tpp = zeros(0)
xoff = zeros(0)
for i=2:ntime
  if (dx[i-1] < 0 && dx[i] > 0) || (dx[i-1] > 0 && dx[i] < 0)
    tpp_0 = time[i-1]-dx[i-1]*(time[i]-time[i-1])/(dx[i]-dx[i-1])
    push!(tpp,tpp_0)
    xbe,ybe,vbe = calc_position(tpp_0,gmstar,elements_b)
    push!(xoff,xbe[1])
  end
end
# Run circular case:
dx = xc0-xb0
tpp0 = zeros(0)
xoff0 = zeros(0)
dvobs = zeros(0)
dxobs = zeros(0)
dx1obs = zeros(0)
dx2obs = zeros(0)
for i=2:ntime
  if (dx[i-1] < 0 && dx[i] > 0) || (dx[i-1] > 0 && dx[i] < 0)
    tpp_0 = time[i-1]-dx[i-1]*(time[i]-time[i-1])/(dx[i]-dx[i-1])
    push!(tpp0,tpp_0)
    xb00,yb00,vb00 = calc_position(tpp_0,gmstar,elements_b0)
    xc00,yc00,vc00 = calc_position(tpp_0,gmstar,elements_c0)
    push!(xoff0,xb00[1])
    push!(dvobs,(vc00[1]-vb00[1])*24.*3600.)
    xbe,ybe,vbe = calc_position(tpp_0,gmstar,elements_b)
    xce,yce,vce = calc_position(tpp_0,gmstar,elements_c)
    push!(dxobs,xce[1]-xbe[1])
    push!(dx1obs,xbe[1]-xb0[1])
    push!(dx2obs,xce[1]-xc0[1])
  end
end
#scatter(tpp,xoff)
#scatter(tpp0,xoff0)
npp = size(tpp)
println(size(tpp))
println(size(tpp0))
#read(STDIN,Char)

# Time to make a diagram:

tb = linspace(t0b,t0b+pb,1000)
tc = linspace(t0c,t0c+pc,1000)
xbd,ybd,vbd = calc_position(tb,gmstar,elements_b)
xbd0,ybd0,vbd0 = calc_position(tb,gmstar,elements_b0)
xcd,ycd,vcd = calc_position(tc,gmstar,elements_c)
xcd0,ycd0,vcd0 = calc_position(tc,gmstar,elements_c0)
# Now, pick a conjunction:
#for iconj=1:10
for iconj=9:9
clf()
#iconj = 1
xbpp,ybpp,vbpp = calc_position(tpp[iconj],gmstar,elements_b)
xbpp=xbpp[1]
ybpp=ybpp[1]
xbpp0,ybpp0,vbpp0 = calc_position(tpp[iconj],gmstar,elements_b0)
xbpp0=xbpp0[1]
ybpp0=ybpp0[1]
xbpp00,ybpp00,vbpp0 = calc_position(tpp0[iconj],gmstar,elements_b0)
xbpp00=xbpp00[1]
ybpp00=ybpp00[1]
xcpp,ycpp,vcpp = calc_position(tpp[iconj],gmstar,elements_c)
xcpp=xcpp[1]
ycpp=ycpp[1]
xcpp0,ycpp0,vcpp0 = calc_position(tpp[iconj],gmstar,elements_c0)
xcpp0=xcpp0[1]
ycpp0=ycpp0[1]
xcpp00,ycpp00,vcpp0 = calc_position(tpp0[iconj],gmstar,elements_c0)
xcpp00=xcpp00[1]
ycpp00=ycpp00[1]
# Now, do some plotting:

clf()
#colors = [ "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
colors = [ "firebrick", "coral", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", ]
linewidth=3
fontsize=10
plot(xbpp/AU+zeros(100),linspace(-semi2/AU*1.05,semi2/AU*1.05,100),ls="dashed",linewidth=linewidth,color="green")
plot(xbpp00/AU+zeros(100),linspace(-semi2/AU*1.05,semi2/AU*1.05,100),ls="dashed",linewidth=linewidth,color="blue")
aw=semi2*.02
x1 = minimum([xbpp,xbpp00])
x2 = maximum([xbpp,xbpp00])
plot([x1,x1+aw,x1,x1+aw,x1,x2,x2-aw,x2,x2-aw]/AU,semi2/AU+[0.,aw,0.,-aw,0.,0.,aw,0.,-aw]/AU+5e-4,linewidth=linewidth,color="black")
annotate(L"\Delta x",xy=[.5*(x1+x2)/AU;semi2/AU+aw/AU+1e-3],fontsize=15,ha="center")
plot(xbd/AU,ybd/AU,linewidth=linewidth,color=colors[1])
plot(xbd0/AU,ybd0/AU,linewidth=linewidth,alpha=0.5,colors[1])
plot(xcd/AU,ycd/AU,linewidth=linewidth,color=colors[2])
plot(xcd0/AU,ycd0/AU,linewidth=linewidth,alpha=0.5,colors[2])
scatter(xbpp/AU,ybpp/AU,color=colors[1],s=64)
scatter(xbpp00/AU,ybpp00/AU,color=colors[1],alpha=0.5,s=64)
#scatter(xbpp00,ybpp00,color=colors[1],s=64)
scatter(xcpp/AU,ycpp/AU,color=colors[2],s=64)
scatter(xcpp00/AU,ycpp00/AU,color=colors[2],alpha=0.5,s=64)
#scatter(xcpp00,ycpp00,color=colors[2],s=64)
scatter([0,0],[0,0],color=colors[8],alpha=0.5,s=256)
xlabel("x [AU]",fontsize=20,fontweight="bold")
ylabel("y [AU]",fontsize=20,fontweight="bold")
#axis([-.02,.02,-.02,.02])
axis("equal")

tight_layout()
savefig("timing_offset_diagram2.pdf",bbox_inches="tight")

tick_params("both",labelsize=15)
## Try to draw an epicycle:
#theta = linspace(0,2pi,100)
#xepb0 = -semi1*eb*cos(theta)
#yepb0 = 2*semi1*eb*sin(theta)
#xepb = xepb0*xbpp0/sqrt(xbpp0^2+ybpp0^2)-yepb0*ybpp0/sqrt(xbpp0^2+ybpp0^2)
#yepb = yepb0*xbpp0/sqrt(xbpp0^2+ybpp0^2)+xepb0*ybpp0/sqrt(xbpp0^2+ybpp0^2)
#plot(xbpp0+xepb,ybpp0+yepb,linewidth=linewidth,color=colors[1],ls="dashed")
## Try to draw an epicycle:
#xepc0 = -semi2*ec*cos(theta)
#yepc0 = 2*semi2*ec*sin(theta)
#xepc = xepc0*xcpp0/sqrt(xcpp0^2+ycpp0^2)-yepc0*ycpp0/sqrt(xcpp0^2+ycpp0^2)
#yepc = yepc0*xcpp0/sqrt(xcpp0^2+ycpp0^2)+xepc0*ycpp0/sqrt(xcpp0^2+ycpp0^2)
#plot(xcpp0+xepc,ycpp0+yepc,linewidth=linewidth,color=colors[2],ls="dashed")
println(iconj)
read(STDIN,Char)
end
