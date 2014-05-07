
d1 = use('/home/zmaw/u241180/ferret/ocbasin-ar5.nc')
basin = d1.gv('basin1')

nls = mpl.colors.BoundaryNorm(np.linspace(0.,1.,11), 80)

m1,cs = mcp.shade(basin.x,basin.y,basin,lon0=35,cmap=cm.Blues,norm=nls,nocb=True)
m1.fillcontinents(color='coral')


x1 = arange(290.,340.)
y1 = arange(39.,51,0.241)

m1.scatter(x1,y1,marker='*',color='black',latlon=False)

