#!/usr/bin/env python3


from pylab import *
from magic import *

def plotSliceOutline(r,theta):

    plot(r[0]*sin(theta),r[0]*cos(theta),'k')
    plot(r[-1]*sin(theta),r[-1]*cos(theta),'k')
    plot([0,0], [r[0],r[-1]],'k')
    plot([0,0], [-r[0],-r[-1]],'k')

nr = 21
ntheta = 64

eta = 0.35

ri = eta/(1-eta)

ro = 1./(1-eta)

r = chebgrid(nr,ri,ro)
theta, gauss = legendre.legendre.gauleg(ntheta)

figure(figsize=(5.5,10))
plotSliceOutline(r,theta)
for k in range(nr):

    plot(r[k]*sin(theta),r[k]*cos(theta),'k:')
    plot(r[k]*sin(theta[ntheta//4]),r[k]*cos(theta[ntheta//4]),'o',color='#1f77b4')
    
axis('off')
tight_layout()
#show()
savefig('rLines.pdf',dpi=400)

figure(figsize=(5.5,10))
plotSliceOutline(r,theta)
for j in range(ntheta):

    plot(r*sin(theta[j]),r*cos(theta[j]),'#1f77b4')

axis('off')
tight_layout()
savefig('thLines.pdf',dpi=400)

figure(figsize=(5.5,10))
plotSliceOutline(r,theta)
for k in range(nr):
    plot(r[k]*sin(theta),r[k]*cos(theta),'#1f77b4')
for j in range(ntheta):
    plot(r*sin(theta[j]),r*cos(theta[j]),'#1f77b4')

axis('off')
tight_layout()
#savefig('GridLines.pdf',dpi=400)
show()
