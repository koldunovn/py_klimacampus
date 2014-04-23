## so you don't have to execute the command plt.show(), run from inside of ipython --pylab
## or ipython --matplotlib

import ferr
import mcmath
import mcp
import numpy as np
import scipy as sp
import matplotlib as mpl
from ferr import use
from matplotlib import cm
import matplotlib.pyplot as plt

d1 = use('sst.v3b.mnmean.nc')

sst = d1.gv('sst')


tax = sst.t
xax = sst.x
yax = sst.y

# text placements
titl_plc = (0.5,1.03)
rtxt_plc = (0.00,0.98)
lett_plc = (0.0,1.02)


# subplot placements
rect_ll = 0.025, 0.08, 0.470,0.26
rect_lr = 0.515, 0.08, 0.470,0.26
rect_ml = 0.025, 0.39, 0.470,0.26
rect_mr = 0.515, 0.39, 0.470,0.26
rect_ul = 0.025, 0.70, 0.470,0.26
rect_ur = 0.515, 0.70, 0.470,0.26

ex=0.025
xe=0.015
sx=0.015

cmapuse = cm.spectral

bnds3 = np.arange(-2,32)
norm3 = mpl.colors.BoundaryNorm(bnds3, 256)

fig = plt.figure(figsize=(6.7,6))

################ A
ax = fig.add_axes(rect_ul)
m1 = mcp.shade(xax,yax,sst[56],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'a.)',ha='left',va='bottom',size=10,transform=ax.transAxes)

################ B
ax = fig.add_axes(rect_ur)
m1 = mcp.shade(xax,yax,sst[159],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'b.)',ha='left',va='bottom',size=10,transform=ax.transAxes)

################ C
ax = fig.add_axes(rect_ml)
m1 = mcp.shade(xax,yax,sst[275],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'c.)',ha='left',va='bottom',size=10,transform=ax.transAxes)

################ D
ax = fig.add_axes(rect_mr)
m1 = mcp.shade(xax,yax,sst[595],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'d.)',ha='left',va='bottom',size=10,transform=ax.transAxes)

################ E
ax = fig.add_axes(rect_ll)
m1 = mcp.shade(xax,yax,sst[1275],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'e.)',ha='left',va='bottom',size=10,transform=ax.transAxes)

################ F
ax = fig.add_axes(rect_lr)
m1 = mcp.shade(xax,yax,sst[1875],lon0=35,ex=ex,xe=xe,sx=sx,norm=norm3,cmap=cmapuse,add=1,cbticksize=10,cllw=0.5,nocb=True)
plt.text(lett_plc[0],lett_plc[1],'f.)',ha='left',va='bottom',size=10,transform=ax.transAxes)


ax = fig.add_axes((0.2,0.03,0.6,0.02))

cb2 = mpl.colorbar.ColorbarBase(ax, cmap=cmapuse,
                                     norm=norm3,
                                     extend='both',
                                     #ticks=bnds3, # optional
                                     #spacing='proportional',
                                     orientation='horizontal')


plt.draw()
