# Simulate two eccentric planets, and time offsets of planet-planet occultation.
include("/Users/ericagol/Courses/PlanetaryDynamics/ExoJulia/ExoJulia/exojulia.jl")
using CGS
using PyPlot
include("regress.jl")

# All we need is the X-coordinate of the orbit:
# Reference direction is sky plane:

function calc_xposition(time,gmstar,elements)
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
  dXdtobs = dXdt.*cos(omega)-dYdt.*sin(omega)
return Xobs,dXdtobs
end

# For specificity, use Trappist-1b/c parameters:
#eb = abs(randn()*.02)
eb = 0.025886097
#omb = 2pi*rand()
omb = -2.0880258
period = 1.51087081
elements_b = [1e-4,period,322.51736,sqrt(eb)*cos(omb),sqrt(eb)*sin(omb)]
p1sec = period*24.*3600.
gmstar = 0.0802*MSUN*GRAV
semi1 = ((p1sec/2/pi)^2*gmstar)^(1./3.)  # semi-major axis in cm if gmstar is cgs
#ec = abs(randn()*.02)
ec = 0.008543
#omc = 2pi*rand()
omc = 0.6692565
period = 2.4218233
elements_c = [1e-4,period,282.80728,sqrt(ec)*cos(omc),sqrt(ec)*sin(omc)]
p2sec = period*24.*3600.
semi2 = ((p2sec/2/pi)^2*gmstar)^(1./3.)  # semi-major axis in cm if gmstar is cgs

t1 = 300.0
t2 = 330.0
time = collect(t1:1e-3:t2)
ntime = length(time)

xb,vb = calc_xposition(time,gmstar,elements_b)
# Now, run with zero eccentricity for comparison:
elements_b0 = copy(elements_b)
elements_b0[4:5]=0.
xb0,vb0 = calc_xposition(time,gmstar,elements_b0)
xc,vc = calc_xposition(time,gmstar,elements_c)
# Now, run with zero eccentricity for comparison:
elements_c0 = copy(elements_c)
elements_c0[4:5]=0.
xc0,vc0 = calc_xposition(time,gmstar,elements_c0)


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
    xbe,vbe = calc_xposition(tpp_0,gmstar,elements_b)
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
    xb00,vb00 = calc_xposition(tpp_0,gmstar,elements_b0)
    xc00,vc00 = calc_xposition(tpp_0,gmstar,elements_c0)
    push!(xoff0,xb00[1])
    push!(dvobs,(vc00[1]-vb00[1])*24.*3600.)
    xbe,vbe = calc_xposition(tpp_0,gmstar,elements_b)
    xce,vce = calc_xposition(tpp_0,gmstar,elements_c)
    push!(dxobs,xce[1]-xbe[1])
    push!(dx1obs,xbe[1]-xb0[1])
    push!(dx2obs,xce[1]-xc0[1])
  end
end
scatter(tpp,xoff)
scatter(tpp0,xoff0)
npp = size(tpp)[1]
npp0 = size(tpp0)[1]
if npp != npp0
  tpp_new = zeros(0)
  tpp0_new = zeros(0)
  i = 0
  while i < npp
    i += 1
    j = 0
    not_found = true
    while j < npp0 && not_found
      j +=1
      if abs(tpp[i]-tpp0[j]) < 0.01*period
        push!(tpp_new,tpp[i])
        push!(tpp0_new,tpp0[j])
        not_found = false
      end
    end
  end
  tpp = tpp_new
  tpp0 = tpp0_new
end
npp = size(tpp)[1]
npp0 = size(tpp0)[1]
println(npp)
println(npp0)
#read(STDIN,Char)

#clf()
dtobs = tpp-tpp0
# Now, use the analytic formula (from 5/15/2017 notes):
# Planet b:
p1 = elements_b[2]
t01 = elements_b[3]
e1 = elements_b[4]^2+elements_b[5]^2
om1 = atan2(elements_b[5],elements_b[4])
# Planet c:
p2 = elements_c[2]
t02 = elements_c[3]
e2 = elements_c[4]^2+elements_c[5]^2
om2 = atan2(elements_c[5],elements_c[4])
alpha = (p1/p2)^(2./3.)
# longitude is M+\omega, but planet crosses sky planet pi/2 after transit:
theta =  -.5*pi-om1 # equation 7.53 from Marion & Thornton
tp = t01 - p1/(2.*pi)*(2.*atan(sqrt((1.-e1)/(1.+e1))*tan(.5*theta))-e1*sqrt(1.-e1^2)*sin(theta)/(1.+e1*cos(theta)))
#lam1 = 2pi/p1.*(tpp0-t01)+pi/2
lam1 = 2pi/p1.*(tpp-t01)-pi/2 #+2*e1*cos(om1)
#lam1 = 2pi/p1.*(tpp0-t01)+pi/2+2*e1*cos(om1)
#lam1 = 2pi/p1.*(tpp0-tp)+om1
theta =  -.5*pi-om1 # equation 7.53 from Marion & Thornton
tp = t02 - p2/(2.*pi)*(2.*atan(sqrt((1.-e2)/(1.+e2))*tan(.5*theta))-e2*sqrt(1.-e2^2)*sin(theta)/(1.+e2*cos(theta)))
#lam2 = 2pi/p2.*(tpp0-t02)+pi/2
lam2 = 2pi/p2.*(tpp-t02)-pi/2 #+2*e2*cos(om2)
#lam2 = 2pi/p2.*(tpp0-t02)+pi/2+2*e2*cos(om2)
#lam2 = 2pi/p2.*(tpp0-tp)+om2
dt1 = p2/(2pi).*(e2.*cos(lam2).*cos(lam2-om2)+2*e2.*sin(lam2).*sin(lam2-om2)-
                alpha*(e1.*cos(lam1).*cos(lam1-om1)+2*e1.*sin(lam1).*sin(lam1-om1)))./
                (sin(lam2)-sin(lam1)/sqrt(alpha))
#dx1 = -semi1*e1*(cos(lam1-om1).*cos(lam1)+2*sin(lam1-om1).*sin(lam1))
#dx2 = -semi2*e2*(cos(lam2-om2).*cos(lam2)+2*sin(lam2-om2).*sin(lam2))
#dx1 = -semi1*e1*(cos(lam1-om1).*cos(lam1)+2*sin(lam1-om1).*sin(lam1))-2*semi1*e1*sin(lam1)*cos(om1)
#dx2 = -semi2*e2*(cos(lam2-om2).*cos(lam2)+2*sin(lam2-om2).*sin(lam2))-2*semi2*e2*sin(lam2)*cos(om2)
# 5/17/2017 notes:
dx1 = semi1*e1*(0.5*cos(2*lam1-om1)-1.5*cos(om1)-2*sin(lam1)*cos(om1))
dx2 = semi2*e2*(0.5*cos(2*lam2-om2)-1.5*cos(om2)-2*sin(lam2)*cos(om2))
dv = 2pi*semi2/p2*(sin(lam2)-sin(lam1)/sqrt(alpha))
dt1 = (dx2-dx1)./dv
#dacc = -4pi^2/p2^2*(cos(lam2)-cos(lam1)/alpha^2)
dacc = -4pi^2*semi2/p2^2*(cos(lam2)-cos(lam1)/alpha^2)
#dt2p = (-dv+sqrt(dv.^2+2*dacc.*(dx2-dx1)))./dacc
#dt2m = (-dv-sqrt(dv.^2+2*dacc.*(dx2-dx1)))./dacc
dt2 = (-dv+dv.*sqrt(abs(1.+2.*dacc.*(dx2-dx1)./dv.^2)))./dacc
#dt2 = zeros(tpp)
#for i=1:npp[1]
#  if dv[i] < 0
#    dt2[i] = dt2m[i]
#  else
#    dt2[i] = dt2p[i]
#  end
#end
#scatter(tpp,dt2,color="r")
#scatter(tpp,dt2m,color="g")

#read(STDIN,Char)
#clf()
#scatter(dtobs,dt)
#plot(dtobs,dtobs)

# Let's fit for eccentricity vectors!

fn = zeros(4,npp[1])
dvnorm = (sin(lam2)-sin(lam1)/sqrt(alpha))*2pi/p2
#fn[1,:]=-alpha*(cos(2*lam1)-2*sin(lam1)-1.5)./dvnorm  # e1*cos(om1) term
#fn[2,:]=-.5*alpha*sin(2*lam1)./dvnorm  # e1*sin(om1) term
#fn[3,:]=(cos(2*lam2)-2*sin(lam2)-1.5)./dvnorm  # e2*cos(om2) term
#fn[4,:]=.5*(2*lam2)./dvnorm  # e2*sin(om2) term
fn[1,:]=-.5*alpha*(-3.+cos(2*lam1)-4.*sin(lam1))./dvnorm  # e1*cos(om1) term
fn[2,:]=-.5*alpha*sin(2*lam1)./dvnorm  # e1*sin(om1) term
fn[3,:]=.5*(-3.+cos(2*lam2)-4.*sin(lam2))./dvnorm  # e2*cos(om2) term
fn[4,:]=.5*sin(2*lam2)./dvnorm  # e2*sin(om2) term
sigma = 1./24./12.  # assume 5-minute uncertainty
noise = randn(npp)*sigma
evec, cov, logdetA = regress(fn,dtobs+noise,sigma+zeros(tpp))
#evec = [eb*cos(omb),eb*sin(omb),ec*cos(omc),ec*sin(omc)]
dtfit = zeros(npp[1])
for i=1:4
  dtfit += evec[i]*fn[i,:]
end
println("e1*cos(om1): ",eb*cos(omb)," ",evec[1]," ",sqrt(cov[1,1]))
println("e1*sin(om1): ",eb*sin(omb)," ",evec[2]," ",sqrt(cov[2,2]))
println("e2*cos(om1): ",ec*cos(omc)," ",evec[3]," ",sqrt(cov[3,3]))
println("e2*sin(om1): ",ec*sin(omc)," ",evec[4]," ",sqrt(cov[4,4]))
tpp -= 300
clf()
fig,axes = subplots(2,1)
subplots_adjust(hspace=0.001)
ax = axes[1]
#ax[:scatter](tpp,(dtobs+noise)*24.*60.,color="r",s=64)
ax[:errorbar](tpp,(dtobs+noise)*24.*60.,color="r",yerr=sigma*24.*60.,fmt="o")
ax[:plot](tpp,dtobs*24.*60.,linewidth=3,alpha=0.8,label="noise free")
ax[:plot](tpp,dtfit*24.*60.,linewidth=3,alpha=0.8,label="total",ls="dashed")
ax[:axis]([0,30,-40,25])
ax[:tick_params]("x",labelsize=0)
ax[:tick_params]("y",labelsize=15)
ax[:set_ylabel]("Timing offset [min]",fontsize=15,fontweight="bold")
ax = axes[2]
ax[:plot](tpp,fn[1,:]*evec[1]*24.*60.,linewidth=2,label=L"e_1 \cos{\omega_1}",ls="solid",alpha=0.8)
ax[:plot](tpp,fn[2,:]*evec[2]*24.*60.,linewidth=2,label=L"e_1 \sin{\omega_1}",ls="dashed",alpha=0.8)
ax[:plot](tpp,fn[3,:]*evec[3]*24.*60.,linewidth=2,label=L"e_2 \cos{\omega_2}",ls="dotted",alpha=0.8)
ax[:plot](tpp,fn[4,:]*evec[4]*24.*60.,linewidth=2,label=L"e_2 \sin{\omega_2}",ls="dashdot",alpha=0.8)
ax[:set_xlabel]("Time [days]",fontsize=15,fontweight="bold")
ax[:set_ylabel]("Timing offset [min]",fontsize=15,fontweight="bold")
ax[:tick_params]("both",labelsize=15)
#scatter(tpp,dt1)
ax[:axis]([0,30,-10,18])
ax[:legend](loc="upper center",fontsize=15,ncol=4)
tight_layout()
savefig("timing_offset_recovered_components_noisy.pdf",bbox_inches="tight")
