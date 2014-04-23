import numpy as np
#import numpy.ma
import scipy as sp
import scipy.signal as sig
import scipy.stats as stats
import datetime
from datetime import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import date2num as nd2n
from netCDF4 import num2date as nn2d

global R_earth
R_earth = 6371000.  ## in meters


class MCMathError(Exception):
	def __init__(self,valus):
		self.value = valus
	def __str__(self):
		return repr(self.value)

global plt

def acorr_mc(x,ret_intts=False):
	"""An \'mcmath\' module version of this function
	This function calculates the normed, positive lag autocorrelation function."""
	if isinstance(x, np.ma.MaskedArray):
		if x.count() <= 3:
			return None
	M = x.size
	R = np.zeros(M-1)
	x = x - x.mean()
	for i in range(M-1):
		R[i] = (x[:M-i]*x[i:]).sum()

	R = R/R[0]
	if bool(ret_intts):
		R1 = R.copy()
		if np.where(R1 < 0)[0].size == 0:
			intts = None
		else:
			atzero = np.where(R1 < 0)[0][0]
			xblah=sp.arange(atzero,R1.size)
			np.put(R1,xblah,0)
			intts = sum(R1)
		return [R,intts]
	else:
		return R

def acorr_mlab(x,endlag):
	"""% This function calculates the autocorrelation function for the 
	% data set 'origdata' for a lag timescale of 0 to 'endlag' and outputs 
	% the autocorrelation function in to 'a'."""

	N = len(x)
	a = []
	for lag in xrange(endlag+1):
		x1 = x[:N-lag]
		x1 = x1 - x1.mean()
		x2 = x[lag:N]
		x2 = x2 - x2.mean()
		a.append(sum(x1*x2)/np.sqrt(sum(x1**2) * sum(x2**2)))

	a = np.array(a)
	return a

def detrend_norm(x):
	from pylab import detrend_linear
	xdt = detrend_linear(x)
	xnrm = xdt / xdt.std()
	return xnrm


def acorr_mskd(x,ret_intts=False):
	return acorr_mc(x,ret_intts)


def blockave(x,blk,ax=0):
	"""Calculate a block average of ndarray 'x'.  I.e., average every 'blk' number of
	elements of 'x' together and return to a new array.  Output array is length x.size / blk
	(any remaining elements on end of x which are less than 'blk' in number are discarded)."""

	xshp = np.array(x.shape)
	newshp = xshp.copy()
	tshp = xshp[0] / blk
	newshp[ax] = tshp
	out = np.ma.masked_all(shape=newshp,dtype=x.dtype)
	for i in xrange(out.shape[ax]):
		out[i] = x[i*blk:i*blk + blk].mean(axis=ax)
	return out


def run_mean_win_mskd(x,win=sig.boxcar(5),cutoff=0):
	"""An \'mcmath\' module function
	Calculates the running mean of (masked) array 'x' using a window 'win', which is a passed array;
	uses a boxcar of width 5 if 'win' is not specified."""
	if isinstance(x, np.ma.MaskedArray):
		x = x.filled(np.NaN)
	N = len(win)
	if N % 2 == 0:
		raise MCMathError("Window must have odd width.")
	X=x.size
	new_arr_len = X - 2*(N / 2)
	new_ind = np.zeros(new_arr_len,dtype=np.int16)
	new_arr = np.zeros(new_arr_len)
	for i in range(new_arr_len):
		ngood = 0
		sum = 0.0
		index = i + N/2
		new_ind[i] = index
		for j in range(N):
			if np.isnan(x[i+j]): continue
			ngood = ngood + 1
			sum = sum + x[i+j]*win[j]
		if ngood == 0:
			new_arr[i] = np.nan
		elif ngood < cutoff:
			new_arr[i] = np.nan
		else:
			new_arr[i] = sum / ngood
	return [new_ind,new_arr]

def run_mean_win_mskd2(x,win=sig.boxcar(5)):
	
	"""An \'mcmath\' module function:
	Just calls 'np.convolve(x,win,mode='valid'), but also returns new index array for
	independent variable as well.
	uses a boxcar of width 5 if 'win' is not specified."""
	if isinstance(x, np.ma.MaskedArray):
		x = x.filled(np.NaN)
	N = len(win)
	if N % 2 == 0:
		raise MCMathError("Window must have odd width.")
	X=x.size
	new_arr_len = X - 2*(N / 2)
	new_ind = np.zeros(new_arr_len,dtype=np.int16)
	for i in  range(new_arr_len):
		index = i + N/2
		new_ind[i] = index
	new_arr = np.convolve(x,win,mode='valid')
	new_arr = np.ma.masked_invalid(new_arr)
	new_arr = new_arr / win.sum()
	return [new_ind,new_arr]

def rm_calc(x,win='b',npts=5,dtrnd=None):
	"""Wrapper for running mean mcmath function 'run_mean_win_mskd2'.
	Calculates running mean over just the first axis (which should
	almost always be time), so that a temporal running mean is calculated for each
	spatial location in 3-d and 4-d datasets passed to this function.

	Parameters: x = data to process
                win = a single letter for window type: 'b'-boxcar, 'h'-hanning,
				      'g'-gaussian, 'm'-hamming, 'p'-parzen (default is boxcar)
                npts = number of points in window (default = 5)
    """
	if not isinstance(x,np.ndarray):
		raise MCMathError("input array must be numpy array.")
	
	N = npts
	if (win == 'b'):
		wind = sig.boxcar(N)
	elif (win == 'h'):
		wind = sig.hanning(N)
	elif (win == 'g'):
		sd = np.int(N*0.16666666667)
		wind = sig.gaussian(N,sd)
	elif (win == 'p'):
		wind = sig.parzen(N)
	elif (win == 'm'):
		wind = sig.hamming(N)
	else:
		raise MCMathError("window type unknown!")

	if x.ndim > 1:
		xshp = np.array(x.shape)
		new_arr_len = xshp[0] - 2*(N / 2)
		rsdshp = xshp[1:]                        ## residual shape (minus axis 0)
		newshp = np.append(new_arr_len,rsdshp)   ## new shape for output array
		out = np.ma.masked_all(newshp,dtype=np.float32)
		out.data[:] = 1e20
		for i in np.ndindex(tuple(rsdshp)):
			rsdind = [slice(None)]
			for xx in xrange(rsdshp.size):
				rsdind.append(i[xx])
			if x[rsdind].mask.all() == True: continue
			if bool(dtrnd):
				ts1 = mpl.mlab.detrend_linear(x[rsdind])
			else:
				ts1 = x[rsdind]
			ind,out[rsdind] = run_mean_win_mskd2(ts1,wind)
			
	else:
		ind,out = run_mean_win_mskd2(x,wind)

	return [ind,out]

def extrema_find(x,kind='all',sd=None,t=None,dt=None):
	"""Returns the indices of array 'x' which have local extrema of kind 'kind', at least 'sd' standard
	deviations away from the mean.

	x: a 1-d ndarray

	kind: 'max','min','all' to return either maxima, minima, or indices for all local extrema. Default is 'all'.

	sd: float which excludes extrema indices from being included in output within 'sd' standard deviations
	of the mean.  Default is to return all extrema indices."""

	diff_2 = np.diff(np.sign(np.diff(x)))

	if ('min' == kind):
		extr_ind = np.ma.where( diff_2 > 0.5 )[0] + 1
	elif ('max' == kind):
		extr_ind = np.ma.where( diff_2 < -0.5 )[0] + 1
	elif ('all' == kind):
		q = np.ma.where(diff_2 > 0.5)[0] + 1
		w = np.ma.where(diff_2 < -0.5)[0] + 1
		extr_ind= np.sort(np.concatenate((q,w)))
	else:
		raise MCMathError("error for 'kind': not recognized, must be one of 'max', 'min', or 'all'.")

	# block for filtering out extrema "too close" to each other along indep. variable
	if (t is not None) and (dt is not None):
		if ('min' == kind):
			pass
		elif ('max' == kind):
			max_ind = extr_ind.copy()
			new_ext_ind = []
			while max_ind.size > 0:
				cur1 = np.where(x[max_ind] == x[max_ind].max())[0][0]
				new_ext_ind.append(max_ind[cur1])
				remain_ind = np.where((t[max_ind] > (t[max_ind][cur1]+dt)) | (t[max_ind] < (t[max_ind][cur1]-dt)))[0]
				max_ind = max_ind[remain_ind]
			extr_ind = np.array(new_ext_ind,dtype=np.int64)
		elif ('all' == kind):
			pass


	# end of block for "closeness" filtering


	# block for s.d. filtering
	if (sd is not None):
		da_mean = np.ma.mean(x)
		da_std = np.ma.std(x)
		new_ind = np.array([],dtype=np.int64)
		sd_lim = sd * da_std
		for i in extr_ind:
			if (abs(x[i] - da_mean) > sd_lim):
				new_ind = np.append(new_ind,i)
		if (new_ind.size > 0):
			extr_ind = new_ind
		else:
			print "error: no surviving extrema for s.d. filter of size %d\n" % sd
			return None
	
	# end of block for s.d. filtering

	return extr_ind

def area_cswt_mean(x,lats,ax=(1,2),units='deg'):
	"""Returns the cos(latitude)-scaled, area-weighted mean of a dataset x, where the dimensions
	being averaged over are evenly-spaced.

    Input: x - the data being area-weight averaged
        lats - the latitudes being averaged over (must correspond to
	         the length of one of the given axes)
        ax - the axes over which to average (usually x- and y- axes,
	         but not always the right after the t-axis)

    Returns: the mean, along other axes not specified"""

	x_shape = np.array(x.shape)
	if max(ax) >= len(x_shape):
		print "Axes given to area_wt_mean() not within data space\n"
		return None

	if ('deg' == units):
		cslats = np.array(np.abs(np.cos(np.deg2rad(lats))))
	elif ('rad' in units):
		cslats = np.array(np.abs(np.cos(lats)))
	else:
		print "function only takes 'deg' or 'rad' for latitude units\n"
		return None

	ax_arr = np.array(ax)
	ave_shape = x_shape[ax_arr]  # 'averaging' shape
	resid_ax = np.ones_like(x_shape)
	resid_ax = resid_ax.astype(bool)
	resid_ax[ax_arr] = False
	resid_shp = x_shape[resid_ax]
	
	## guarantee a reshaped 'cslats' has the proper
	## dimensions to multiply 'x'
	reshape_cslats = np.ones(len(ave_shape)) - 2
	rshp = np.where(ave_shape == len(cslats),reshape_cslats,1)
	oo = np.ones_like(rshp)
	if (oo == rshp).all():
		print "\'lats\' given to area_wt_mean() not of proper length for \'x\'\n"
		return None
	cslats = cslats.reshape(rshp)

	result_arr = np.zeros(resid_shp)
	result_arr = np.ma.MaskedArray(result_arr,mask=False)

	## Main Loop: loop through all dims not being averaged over,
	##   yielding an averaged slice per point in those dims
	for resind in np.ndindex(tuple(resid_shp)):
		x_ind_list = []
		ind_loop = 0
		for x_ind in range(len(resid_ax)):
			if resid_ax[x_ind]:
				x_ind_list.append(resind[ind_loop])
				ind_loop = ind_loop + 1
			else:
				x_ind_list.append(slice(None))

		# get the slice of the data over which to average
		x_slice = x[x_ind_list]

		# make the cosine-scaled latitudes array masked
		# so that a simple .sum() adds up only the pieces which went into scaling the data
		cslats_mskd = cslats*np.ones_like(x_slice)

		csmod_vals = x_slice*cslats_mskd
		ans_slice = csmod_vals.sum() / cslats_mskd.sum()
		if ans_slice.all():
			result_arr.mask[resind] = True
		result_arr[resind] = ans_slice

	return result_arr

def area_wt_mean(x,area,axes=(1,2),units='deg'):
	"""WORK IN PROGRESS!  This function calculates the area-weighted mean for the x,y part
	of a data set.  One must provide the data and the grid area data.  Axes for x,y must
	be specified at the moment, but should be automated to find the two x,y axes in the
	given data in the future."""

	if axes == (1,2):
		x_w = x * area.reshape(1,area.shape[0],area.shape[1])
		x_wm = x_w.sum(axis=2).sum(axis=1) / area.sum()
	elif axes == (0,1):
		x_w = x * area
		x_wm = x_w.sum() / area.sum()
	return x_wm

##
#### BEGIN corr2d

def corr2d(x,y):
	"""Given two 3-d arrays 'x' and 'y', which are both dimensionally (T,Y,X) or at least are 3-d where the
	correlation in the first dimension is desired for each Y,X 'spot', returns a 2-d array with (Y,X) dimensions
	with the T correlation value at each 2-d point."""

	if x.shape != y.shape:
		raise MCMathError("x and y arrays must be of the same shape")

	tsz,ysz,xsz = x.shape
	out = np.ma.masked_all((ysz,xsz))
	for (i,j) in np.ndindex(ysz,xsz):
		if isinstance(x, np.ma.MaskedArray):
			if (x[:,i,j].count() == 0):
				continue
		out[i,j] = np.corrcoef(x[:,i,j],y[:,i,j])[0,1]

	return out

#### END corr2d
##


##
#### BEGIN trnd_calc

def trnd_calc(t,x,bs=None,pval=0.025,intts=None):
	if isinstance(t[0],dt):
		t = d2n(t)
	if isinstance(x, np.ma.MaskedArray):
		if x.mask.all() == True:
			return None
		if x.mask.any() == True:
			x1 = x.data[~x.mask]
			t1 = t[~x.mask]
		else:
			x1 = x.data
			t1 = t
	else:
		x1 = x
		t1 = t
	if t1.size < 3: return None

	if (bs is not None) & (intts is None):
		tacorr1 = acorr_mskd(x)
		if np.where(tacorr1 < 0)[0].size == 0:
			return None
		atzero = np.where(tacorr1 < 0)[0][0]
		xblah=sp.arange(atzero,tacorr1.size)
		np.put(tacorr1,xblah,0)
		intts = sum(tacorr1)

	out1 = stats.linregress(t1,x1)
	intercep1 = out1[1]
	slope1 = out1[0]
	r1 = out1[2]
	P1 = out1[3]
	SE_est1 = out1[4]
	yhat = slope1*t1 + intercep1
	ssy = stats.ss(yhat - x1)
	ssx = stats.ss(t1 - t1.mean())

	if bs is not None:
		dof = t1.size / intts
		if dof < 3:
			dof = 3
		SE_est1 = SE_est1*(sp.sqrt((t1.size-2.)/(dof-2.)))
		P1 = stats.t.sf(abs(slope1) / SE_est1, dof-2) *2

		SE1_slope = SE_est1*stats.t.isf(pval,dof-2)
		#SE1_slope = sp.sqrt(ssy/(x1.size - 2)) / sp.sqrt(ssx)  # same as stats.linregress' SE_est1
		#P1_slope = stats.t.sf(abs(slope1) / SE1_slope, dof-2) * 2  #WRONG!!

	delT1 = slope1*(t1[-1] - t1[0])

	if bs is None:
		return (slope1,intercep1,r1,SE_est1,P1,delT1)
	else:
		return (slope1,intercep1,r1,SE_est1,P1,SE1_slope,delT1)

#### END trnd_calc
##

##
#### BEGIN trnd_map_calc

def trnd_map_calc(t,x,ind_slice=slice(None,None,None),fullstat=None,intts=None):
	"""Calculates trends across a 3-d (t,y,x) map to yield
	a 2-d map (y,x) of trends across the t-dimension.  Units will
	be property-units per time unit (e.g., degC per day) if t is
	datetime instances, or in whatever time units 't' is in if
	provided directly from the source netCDF data file.

	If fullstat is not None, then the return array has all the returned values
	from trnd_calc(), otherwise, just the slope.

	If intts is not None, it implies fullstat is not None;
	intts has to be passed an array the same size as the input array
	without the time axis (usually 2-d, like y,x), with the integral time scale
	being the value at each grid point.
	
	(ind_slice still not implemented)"""
	
	if isinstance(t[0],dt):
		t = d2n(t)

	xshp = x.shape
	if len(xshp) != 3:
		raise MCMathError('size of input array must be 3-dim (t,y,x)')

	if fullstat is None:
		trnd = np.ma.masked_all((xshp[1:]),dtype=x.dtype)
	else:
		outshp = list(xshp)
		outshp[0] = 7
		outshp = tuple(outshp)
		trnd = np.ma.masked_all(outshp,dtype=x.dtype)

	for i in np.ndindex(xshp[1:]):
		if isinstance(x, np.ma.MaskedArray):
			if x[:,i[0],i[1]].mask.all() == True:
				continue
			if x[:,i[0],i[1]].mask.any() == True:
				x1 = x[:,i[0],i[1]].data[~x[:,i[0],i[1]].mask]
				t1 = t[~x[:,i[0],i[1]].mask]
			else:
				x1 = x[:,i[0],i[1]].data
				t1 = t
		else:
			x1 = x[:,i[0],i[1]]
			t1 = t

		if t1.size < 3: continue

		if fullstat is None:
			blah = trnd_calc(t1,x1)
			if blah is not None:
				trnd[i[0],i[1]] = blah[0]
		else:
			if intts is None:
				blah = trnd_calc(t1,x1,bs=1)
				if blah is not None:
					trnd[:,i[0],i[1]] = np.array(blah)
			else:
				if intts.shape != xshp[1:]:
					raise MCMathError('size of intts must be the same as the input array minus the first dimension (time axis)')
				blah = trnd_calc(t1,x1,bs=1,intts=intts[i[0],i[1]])
				if blah is not None:
					trnd[:,i[0],i[1]] = np.array(blah)


	return trnd

#### END trnd_map_calc
##

##
#### BEGIN run_trnd_map

def run_trnd_map(t,x,tlen):
	"""Calculates a series of running trend magnitudes over a time interval
	of length 'tlen', and returns an output file with the magnitudes along the
	first axis.  [Previously, to avoid highly repetitious overlap, the intervals
	only overlapped by one quarter of their overall length (e.g., a 100-time point
	running trend map calculates the 100-timepoint trend every 25 time points).  Now,
	this is being actively changed, is currently 5, but should make this a passable arg.]"""

	#intvl = tlen / 4
	intvl = 5
	numtrnd = ((x.shape[0] - tlen) / intvl) + 1

	outarr = np.ma.masked_all((numtrnd,x.shape[1],x.shape[2]))

	for i in xrange(numtrnd):
		outarr[i] = trnd_map_calc(t[i*intvl:i*intvl+tlen],x[i*intvl:i*intvl+tlen])

	return t[:-tlen+1:intvl],outarr


#### END run_trnd_map
##

##
#### BEGIN trnd_cmp_test

def trnd_cmp_test(t1,x1,t2,x2,pval=0.025,bs=None,eqs=None):

	if isinstance(x1, np.ma.MaskedArray):
		if x1.mask.any() == True:
			x11 = x1.data[~x1.mask]
			t1 = t1[~x1.mask]
		else:
			x11 = x1.data
	else:
		x11 = x1

	if isinstance(x2, np.ma.MaskedArray):
		if x2.mask.any() == True:
			x22 = x2.data[~x2.mask]
			t2 = t2[~x2.mask]
		else:
			x22 = x2.data
	else:
		x22 = x2

	## Don't bother "fitting" a line through 1 or 2 points
	if t1.size < 3: return None
	if t2.size < 3: return None

	tacorr1 = acorr_mskd(x11)
	if np.where(tacorr1 < 0)[0].size == 0:
		return None
	atzero = np.where(tacorr1 < 0)[0][0]
	xblah=sp.arange(atzero,tacorr1.size)
	np.put(tacorr1,xblah,0)
	intts1 = sum(tacorr1)

	tacorr2 = acorr_mskd(x22)
	if np.where(tacorr2 < 0)[0].size == 0:
		return None
	atzero = np.where(tacorr2 < 0)[0][0]
	xblah2=sp.arange(atzero,tacorr2.size)
	np.put(tacorr2,xblah2,0)
	intts2 = sum(tacorr2)

	out1 = stats.linregress(t1,x11)
	intercep1 = out1[1]
	slope1 = out1[0]
	r1 = out1[2]
	P1 = out1[3]
	SE_est1 = out1[4]
	yhat = slope1*t1 + intercep1
	ssy = stats.ss(yhat - x11)
	ssx = stats.ss(t1 - t1.mean())

	out2 = stats.linregress(t2,x22)
	intercep2 = out2[1]
	slope2 = out2[0]
	r2 = out2[2]
	P2 = out2[3]
	SE_est2 = out2[4]
	yhat2 = slope2*t2 + intercep2
	ssy2 = stats.ss(yhat - x22)
	ssx2 = stats.ss(t2 - t2.mean())

	if bs is not None:
		if eqs is not None:
			alph1 = alph_est(x11,'MPK',t1.size)
			dof1 = eqn(t1.size,alph1)
			alph2 = alph_est(x22,'MPK',t2.size)
			dof2 = eqn(t2.size,alph2)
		else:
			dof1 = t1.size / intts1
			dof2 = t2.size / intts2

		if (dof1 == 2):
			dof1 = 3
		if (dof2 == 2):
			dof2 = 3

		SE_est1 = SE_est1*(sp.sqrt((t1.size-2.)/(dof1-2.)))
		SE_est2 = SE_est2*(sp.sqrt((t2.size-2.)/(dof2-2.)))

##		print "original d.o.f.1: %d" % dof1
##		print "original d.o.f.2: %d" % dof2

	else:
		dof1 = t1.size
		dof2 = t2.size

	P1 = stats.t.sf(abs(slope1) / SE_est1, dof1 - 2) * 2
	P2 = stats.t.sf(abs(slope2) / SE_est2, dof2 - 2) * 2
	SE12_est = sp.sqrt(SE_est1**2 + SE_est2**2)
	t12 = abs(slope1 - slope2) / SE12_est
	P12 = stats.t.sf(t12,dof1 - 2) * 2
	err1_slope = SE_est1*stats.t.isf(pval,dof1 - 2)
	err2_slope = SE_est2*stats.t.isf(pval,dof2 - 2)
	err12 = SE12_est*stats.t.isf(pval,dof1 - 2)

##	if P12 > 0.05:
##		Ptest = P12
##		new_dof = dof1
##		while Ptest > 0.05:
##			new_dof = new_dof + 1
##			Ptest = stats.t.sf(t12,new_dof - 2) * 2
##		print "original d.o.f.: %d" % dof1
##		print "# d.o.f needed to get p-value of 0.05: %d" % (new_dof)

	return [P12,slope1,err1_slope,slope2,err2_slope,t12,SE12_est,dof1,dof2,err12]

#### END trnd_cmp_test
##

##
#### BEGIN my_dtrnd

def my_dtrnd(x,t=None):
	"""My version of a linear detrend of 3-d data (will add better code
	for 4-d data later), input as (T,Y,X) data.  This version, unlike the
	Numpy / Scipy code, can detrend a time series with masked holes properly
	(i.e., no slopes on the order of 1e13 or all NaNs), by the same method
	as trnd_calc(), which provides (t,x) data pairs for only those unmasked
	data points.

	Input:  x - data to detrend
	           t - time axis, if data.t does not contain it

	Output: dtrnd (same size as original array)
	"""
	if t is None:
		if hasattr(x,'t'):
			t = x.t
		else:
			raise MCMathError("no time axis specified for trend calculation")
	
	if isinstance(t[0],dt):
		t = d2n(t)

	xshp = x.shape
	if len(xshp) != 3:
		raise MCMathError('size of input array must be 3-dim (t,y,x)')

	dtrnd = np.ma.masked_all((xshp),dtype=x.dtype)

	for i in np.ndindex(xshp[1:]):
		if isinstance(x, np.ma.MaskedArray):
			if x[:,i[0],i[1]].mask.all() == True:
				continue
			if x[:,i[0],i[1]].mask.any() == True:
				x1 = x[:,i[0],i[1]].data[~x[:,i[0],i[1]].mask]
				t1 = t[~x[:,i[0],i[1]].mask]
			else:
				x1 = x[:,i[0],i[1]].data
				t1 = t
		else:
			x1 = x[:,i[0],i[1]]
			t1 = t

		if t1.size < 3: continue

		blah = trnd_calc(t1,x1)
		if blah is not None:
			trendline = blah[0]*t + blah[1]
			dtrnd[:,i[0],i[1]] = x[:,i[0],i[1]] - trendline

	return dtrnd

#### END my_dtrnd
##


##
#### BEGIN mean_diff_test

def mean_diff_test(x1,x2,pval=0.025,bs=None,eqs=None,alph=None,intts=None):

	if isinstance(x1, np.ma.MaskedArray):
		if x1.mask.any() == True:
			x11 = x1.data[~x1.mask]
		else:
			x11 = x1.data
	else:
		x11 = x1

	if isinstance(x2, np.ma.MaskedArray):
		if x2.mask.any() == True:
			x22 = x2.data[~x2.mask]
		else:
			x22 = x2.data
	else:
		x22 = x2

	if (intts is None) & (eqs is None):
		tacorr1,intts1 = acorr_mc(x11,True)
		tacorr2,intts2 = acorr_mc(x22,True)
	else:
		intts1 = intts
		intts2 = intts

	if bs is not None:
		if eqs is not None:
			if alph is None:
				alph1 = alph_est(x11,'MPK',x11.size)
				dof1 = eqn(x11.size,alph1)
				alph2 = alph_est(x22,'MPK',x22.size)
				dof2 = eqn(x22.size,alph2)
			else:
				dof1 = eqn(x11.size,alph)
				dof2 = eqn(x22.size,alph)
		else:
			dof1 = x11.size / intts1
			dof2 = x22.size / intts2

	else:
		dof1 = x11.size
		dof2 = x22.size

	if (dof1 <= 2):
		dof1 = 2
	if (dof2 <= 2):
		dof2 = 2

	mindof = min(dof1,dof2)
	mean1 = x11.mean()
	mean2 = x22.mean()
	meandiff = mean2 - mean1
	SE1_mean = x11.var(ddof=1) / dof1
	SE2_mean = x22.var(ddof=1) / dof2
	SE12_mean = sp.sqrt(SE1_mean + SE2_mean)

	# just for fun
	nu = (SE1_mean + SE2_mean)**2 / ((SE1_mean**2/(dof1-1)) + (SE2_mean**2 / (dof2-1)))

	if nu < 1:
		nu = 1
	
	t12 = abs(mean2 - mean1) / SE12_mean
	P12 = stats.t.sf(t12,nu) * 2
	err12 = SE12_mean*stats.t.isf(pval,nu)
	Ptest = stats.t.sf(t12,nu+1) * 2
	Ptest2 = stats.t.sf(t12,nu+2) * 2
	Ptest3 = stats.t.sf(t12,nu+3) * 2

	return np.array([meandiff,P12,err12,dof1,dof2,SE1_mean,SE2_mean,SE12_mean,nu,Ptest,Ptest2,Ptest3])

#### END mean_diff_test
##


##
#### BEGIN bmapinterp

def bmapinterp(data,xin1,yin1,xin2,yin2,**kwargs):
	from mpl_toolkits.basemap import interp as bminterp

	if (xin1 > 360).any():
		xin1_0_360 = np.where(xin1 > 360, xin1 - 360, xin1)
		roll_val = xin1_0_360.size -  np.where(xin1_0_360 - np.roll(xin1_0_360,1) < 0)[0]
		xin1 = np.roll(xin1_0_360,roll_val)
		data = np.roll(data,roll_val,axis=1)
		
	if (xin2 > 360).any():
		xin2_0_360 = np.where(xin2 > 360, xin2 - 360, xin2)
		roll_val = xin2_0_360.size -  np.where(xin2_0_360 - np.roll(xin2_0_360,1) < 0)[0]
		xin2 = np.roll(xin2_0_360,roll_val)

	(x2,y2) = np.meshgrid(xin2,yin2)
	new_data_1 = bminterp(data,xin1,yin1,x2,y2,**kwargs)
	new_data_0 = bminterp(data,xin1,yin1,x2,y2,order=0,**kwargs)
	new_data = np.ma.where(new_data_1.mask == True,new_data_0,new_data_1)

	##return (new_data_0,new_data_1)
	return (new_data,xin2,yin2)

#### END bmapinterp
##

##
#### BEGIN edges

def edges(xin,mod=None):
	if xin.size <= 1:
		raise MCMathError('size of input array has to be larger than 1')

	xedg = np.ndarray(shape=xin.size+1)
	for i in range(xedg.size):
		if i == 0:
			if mod is not None:
				xedg[i] = ((xin[-1] - mod) + xin[0]) / 2
			else:
				xedg[i] = xin[0] - ((xin[1] - xin[0]) / 2)
			continue
		if i < (xedg.size - 1):
			xedg[i] = xin[i-1] + (xin[i] - xin[i-1]) / 2
		else:
			if mod is not None:
				xedg[i] = xedg[0] + mod
			else:
				xedg[i] = xin[i-1] + ((xin[i-1] - xin[i-2]) / 2)

	return xedg.copy()

#### END edges
##

##
#### BEGIN intvls

def intvls(xin,edgs=None):
	if xin.size <= 1:
		raise MCMathError('size of input array has to be larger than 1')

	xedg = edges(xin)
	xdiff = xedg[1:] - xedg[0:-1]
	return xdiff

#### END intvls
##

##
#### BEGIN basmask

def basmask(bas,grid='edr20',returnxy=True):
	import ferr
	dtmp = ferr.use('/pf/u/u241180/py/basin.nc',silent=1)
	if grid == 'edr20':
		v = 'basin_edr_20'
	elif grid == 'edr0':
		v = 'basin_edr_0'
	elif grid == '1dr':
		v = 'basin_1_0'
	else:
		print "grid '%s' not found\n" % grid
		return None

	basin = dtmp.v[v][:]
	baslon = dtmp.getax(v,'x')
	baslat = dtmp.getax(v,'y')
	dtmp.f.close()
	
	if bas == 'atl':
		
		atl = np.ma.where(basin != 0.,np.ma.masked, 1.)
		atl.grid = grid
		if returnxy is True:
			atl.x = baslon
			atl.y = baslat
		return atl
	
	elif bas == 'pac':
		
		pac = np.ma.where(basin != 1.,np.ma.masked, 1.)
		pac.grid = grid
		if returnxy is True:
			pac.x = baslon
			pac.y = baslat
		return pac
	
	elif bas == 'ind':

		ind = np.ma.where(basin != 4.,np.ma.masked, 1.)
		ind.grid = grid
		if returnxy is True:
			ind.x = baslon
			ind.y = baslat
		return ind

	else:
		basin.grid = grid
		if returnxy is True:
			basin.x = baslon
			basin.y = baslat
		return basin

#### END basmask
##

def gapmax(data):
	if data.ndim > 1:
		print 'data needs to be a 1-d array\n'
		return None
	gaptmp = 0
	gapbig = 0
	for i in xrange(data.size):
		if data.mask[i] == True:
			gaptmp = gaptmp + 1
		else:
			gaptmp = 0
		if gaptmp > gapbig:
			gapbig = gaptmp

	return gapbig


def make_clim(data,ax=0,stride=12):
	shp = np.array(data.shape)
	shp[ax] = stride
	clim = np.ma.masked_all(shp,dtype=data.dtype)
	
	for j in xrange(stride):
		slicer = []
		for i in xrange(len(shp)):
			spot = j
			if i == ax:
				slicer.append(slice(spot,None,stride))
			else:
				slicer.append(slice(None))

		binned = data[slicer]
		slicer[ax] = j
		clim[slicer] = binned.mean(axis=ax)

	return clim

def partial_sum_clim(data,ax=0,stride=12):
	shp = np.array(data.shape)
	shp[ax] = stride
	clim = np.ma.masked_all(shp)
	
	for j in xrange(stride):
		slicer = []
		for i in xrange(len(shp)):
			spot = j
			if i == ax:
				slicer.append(slice(spot,None,stride))
			else:
				slicer.append(slice(None))

		binned = data[slicer]
		slicer[ax] = j
		clim[slicer] = binned.sum(axis=ax)

	return clim

def anom_from_clim(data,clim):
	tlength = clim.shape[0]
	dlength = data.shape[0]
	tax_clim = dlength / tlength + 1

	if data.ndim == 1:
		newclim = np.tile(clim,tax_clim)
	if data.ndim == 3:
		newclim = np.tile(clim,((tax_clim),1,1))
	if data.ndim == 4:
		newclim = np.tile(clim,((tax_clim),1,1,1))
	newclim = newclim[:dlength]
	anom = data - newclim
	return anom


def ts_median(d,t,length):
	medt = np.zeros(t.size - length, dtype=t.dtype)
	dmed = np.zeros(d.size - length, dtype=d.dtype)
	for i in xrange(t.size - length):
		end = i+length
		medt[i] = np.ma.median(t[i:end])
		dmed[i] = d[i+int(length/2)]

	return (dmed,medt)



def eof(x):
	import time
	xshp = x.shape
	if len(xshp) != 3:
		print "error: x must be a (t,y,x) data array)\n"
		return None
	print "mcmath.eof() started: %s" % (time.ctime())
##	out1 = np.zeros(shape=(xshp[0],1))
##	loc1 = np.array([(0,0)])
##	firstwrite = 0
##	for i in np.ndindex(xshp[1:3]):
##		ts = x[:,i[0],i[1]]
##		##print i,ts
##		if ts.mask.any() == True:
##			#print "all mask tripped"
##			continue
##		if firstwrite == 0:
##			##print "firstwrite tripped"
##			out1[:,0] = ts
##			loc1[0] = i
##			firstwrite = 1
##			continue
##		i1 = np.reshape(np.array(i),(1,2))
##		ts = np.reshape(ts,(ts.size,1))
##		##print "in main body of for loop\n"
##		out1 = np.append(out1,ts,axis=1)
##		loc1 = np.append(loc1,i1,axis=0)

	out1,loc1 = tyx2tm_arr(x)
	
	print "out1 size is: %d,%d" % (out1.shape)
	print "Done collecting matrix; now performing SVD..."
	print "%s\n" % (time.ctime())
	[U,s,V] = sp.linalg.svd(out1,full_matrices=0)
	## The rows of V (this 'svd' returns V' actually, so eigenvectors are rows) are the
	## spatial eigenvectors (the way I have it constructed here).  The columns of U are
	## the PC time series corresponding to the spatial eigenvector rows of V.

	print "Done calculating SVD; now parsing results back into x,y and t arrays"
	print "%s" % (time.ctime())
	neofs = min(s.size,15)

	space = np.ma.masked_all(shape=(neofs,xshp[1],xshp[2]),dtype=np.float32)
	tim1 = np.ma.masked_all(shape=(xshp[0],neofs),dtype=np.float32)
	
	tim1[:,0:neofs] = U[:,0:neofs].copy()
	for eof in xrange(neofs):
		ctr = 0
		for i,j in loc1:
			space[eof,i,j] = V[eof,ctr].copy()
			ctr = ctr + 1

	space.data[space.mask] = space.get_fill_value()
	print "Completed... ",
	print "%s" % (time.ctime())
	
	return space,s,tim1

def tyx2tm_arr(x):
	"""Takes (t,y,x) array and returns (t,M) array where M are the columns
	representing the spatial locations of the time series concatenated into one
	row per time point.

	Returns: out1 (the t x M array) and loc1 (an M x 2 array carrying the
             original indices of the y,x locations of the columns)"""
	
	xshp = x.shape
	if len(xshp) != 3:
		print "error: x must be a (t,y,x) data array)\n"
		return None

	##out1 = np.zeros(shape=(xshp[0],1),dtype=np.float32)
	loc1 = np.array([(0,0)],dtype=np.int16)
	firstwrite = 0
	for i in np.ndindex(xshp[1:3]):
		ts = x[:,i[0],i[1]]
		##print i,ts
		if ts.mask.any() == True:
			#print "all mask tripped"
			continue
		if firstwrite == 0:
			##print "firstwrite tripped"
			##out1[:,0] = ts
			loc1[0] = i
			firstwrite = 1
			continue
		i1 = np.reshape(np.array(i),(1,2))
		##ts = np.reshape(ts,(ts.size,1))
		##print "in main body of for loop\n"
		##out1 = np.append(out1,ts,axis=1)
		loc1 = np.append(loc1,i1,axis=0)

	out1 = x[:,loc1[:,0],loc1[:,1]]

	return out1,loc1

def tm2tyx_arr(x,loc1,outshp):
	"""reconstructs (t,y,x) array from (t,M) array when given an (M,2) array of location indices.
	Also needs 'outshp' shape to construct the proper output array's (y,x) dimensions. """
	if (outshp[0] * outshp[1]) < loc1.shape[0]:
		raise MCMathError("Output shape not large enough to contain the number of locations in loc1 index array.")
	if x.shape[-1] != loc1.shape[0]:
		raise MCMathError("Input array does not have the same number of spatial entries as loc1 index array.")
	if len(x.shape) == 1:
		firstdim = 1
		x = x.reshape(1,-1)
	else:
		firstdim = x.shape[0]
	
	out1 = np.ma.masked_all(shape=(firstdim,outshp[0],outshp[1]),dtype=np.float32)
	ctr = 0
	for i,j in loc1:
		out1[:,i,j] = x[:,ctr]
		ctr = ctr + 1

	out1 = np.squeeze(out1)
	return out1

def pt_corrmap(pt1,x):
	"""Input array 'x', with dims ('t,y,x'), is correlated with a single time series 'pt1' at all valid spatial locations.
	Requires (at present) that the input array doss not have any missing values in the time series;
	it can have missing values for an entire spatial location (i.e., no time values there at all). """
	
	newarr,loc1 = tyx2tm_arr(x)
	pt_arr = np.insert(newarr,0,pt1,axis=1)
	cc1 = np.corrcoef(pt_arr,rowvar=0)[0]
	out1 = tm2tyx_arr(cc1[1:],loc1,x.shape[1:])

	return out1

	
def area_vol_calc(depth=750,grid='edr20'):
	from ferr import use

	mperdeg = np.deg2rad(1.)*R_earth   ## meters per degree latitude
	sqdeg_box = 2.5 ** 2  ## square degrees per lon/lat grid box

	all = basmask('all',returnxy=True)
	abmarg = np.ma.masked_values(all,3.0)
	allmaj = np.ma.masked_values(abmarg,2.0)
	cslats = np.abs(np.cos(np.deg2rad(all.y)))
	cslats = cslats.reshape((all.y.size,1))
	all_ones = np.ones_like(allmaj)
	sq_deg_xy = all_ones*sqdeg_box
	area = (sq_deg_xy*cslats*all_ones)*(mperdeg ** 2)

	drose = use('/pf/u/u241180/py/etopo_edr.nc',silent=1)
	rose = drose.v['rose_edr']
	rosedepth = np.where(rose[:] < -depth, -depth, rose)
	rosedeptha = np.ma.masked_greater(rosedepth,0)
	rosedeptha = np.abs(rosedeptha)

	vol = area*rosedeptha

	area.x = all.x
	area.y = all.y
	vol.x = all.x
	vol.y = all.y
	
	return (area,vol)

def area_calc(x,y):
	"""Calculates area of a REGULAR grid from the input list of longitudes (x) and latitudes (y).
	Returns a 2-d grid of square meters with dimensions derived the input grid: size(y) x size(x)."""

	mperdeg = np.deg2rad(1.)*R_earth   ## meters per degree latitude
	
	xedgs = edges(x)
	xsiz = xedgs[1:] - xedgs[:-1]

	yedgs = edges(y)
	ysiz = yedgs[1:] - yedgs[:-1]

	# this should do nothing but "fix" y-coordinate-box sizes derived
	# from decreasing y-values instead of increasing y-values
	ysiz = abs(ysiz)
	
	x1,y1 = np.meshgrid(xsiz,ysiz)
	xblah,ylatgrid = np.meshgrid(x,y)
	cosscal = sp.special.cosdg(ylatgrid)
	xscald = x1*cosscal

	area1 = xscald*y1

	area_m2 = area1 * (mperdeg)**2

	return area_m2


def regs(trop=20,sp=60):
	atl = basmask('atl',returnxy=True)
	pac = basmask('pac')
	ind = basmask('ind')
	lon = atl.x
	lat = atl.y
	n_in = np.where(lat >= trop)[0]
	nsp_in = np.where(lat >= sp)[0]
	nml_in = np.where((lat >= trop) & (lat < sp))[0]
	t_in = np.where((lat >= -trop) & (lat <= trop))[0]
	s_in = np.where(lat <= -trop)[0]
	sml_in = np.where((lat <= -trop) & (lat > -sp))[0]
	ssp_in = np.where(lat <= -sp)[0]

	n_bas = np.ma.masked_all(atl.shape)
	n_bas[n_in,:] = 1.

	nsp_bas = np.ma.masked_all(atl.shape)
	nsp_bas[nsp_in,:] = 1.

	nml_bas = np.ma.masked_all(atl.shape)
	nml_bas[nml_in,:] = 1.

	t_bas = np.ma.masked_all(atl.shape)
	t_bas[t_in,:] = 1.

	s_bas = np.ma.masked_all(atl.shape)
	s_bas[s_in,:] = 1.

	ssp_bas = np.ma.masked_all(atl.shape)
	ssp_bas[ssp_in,:] = 1.

	sml_bas = np.ma.masked_all(atl.shape)
	sml_bas[sml_in,:] = 1.


	natl = n_bas*atl
	nmlatl = nml_bas*atl
	nspatl = nsp_bas*atl
	tatl = t_bas*atl
	satl = s_bas*atl
	smlatl = sml_bas*atl
	sspatl = ssp_bas*atl

	npac = n_bas*pac
	nmlpac = nml_bas*pac
	nsppac = nsp_bas*pac
	tpac = t_bas*pac
	spac = s_bas*pac
	smlpac = sml_bas*pac
	ssppac = ssp_bas*pac

	nind = n_bas*ind
	tind = t_bas*ind
	sind = s_bas*ind
	smlind = sml_bas*ind
	sspind = ssp_bas*ind

	q = EmbArray()
	q.lon = atl.x
	q.lat = atl.y
	q.atl = atl
	q.pac = pac
	q.ind = ind
	q.natl = natl
	q.nmlatl = nmlatl
	q.nspatl = nspatl
	q.tatl = tatl
	q.satl = satl
	q.smlatl = smlatl
	q.sspatl = sspatl

	q.npac = npac
	q.nmlpac = nmlpac
	q.nsppac = nsppac
	q.tpac = tpac
	q.spac = spac
	q.smlpac = smlpac
	q.ssppac = ssppac

	q.nind = nind
	q.tind = tind
	q.sind = sind
	q.smlind = smlind
	q.sspind = sspind
	
	q.n_in = n_in
	q.nsp_in = nsp_in
	q.nml_in = nml_in
	q.t_in = t_in
	q.s_in = s_in
	q.sml_in = sml_in
	q.ssp_in = ssp_in


	return q




class EmbArray(object):
	def __init__(self):
		self.data = np.ndarray(shape=(0,),dtype=np.float64)
		return None
	
	def totalsize(self):
		tot = 0
		for i in self.__dict__.keys():
			iatt = getattr(self,i)
			if hasattr(iatt,'size'):
				tot += iatt.size
			elif hasattr(iatt,'__len__'):
				tot += len(iatt)
		
		return tot
	
	def totalbytes(self):
		tot = 0
		for i in self.__dict__.keys():
			iatt = getattr(self,i)
			if hasattr(iatt,'nbytes'):
				tot += iatt.nbytes
			elif hasattr(iatt,'__len__'):
				tot += len(iatt)
		
		return tot
	
	def __repr__(self):
		s = "Contents:\n"
		for i in np.sort(self.__dict__.keys()):
			iatt = getattr(self,i)
			if isinstance(iatt,np.ma.MaskedArray):
				s += "%s:\tMaskedArray, shape: %s\n" % (i,iatt.shape)
			elif isinstance(iatt,np.ndarray):
				s += "%s:\tArray     , shape: %s\n" % (i,iatt.shape)

		s += "\nTotal size (in bytes): %d\n" % self.totalbytes()
		return s
	
##	def __len__(self):
##		return 10


def edr0_20(data,lon,ax=2):
	"""Function which rolls data and lon to edr grid starting at 20E
	when original lon starts at 0E
	Assumes data is in (t,y,x) dimensions.  If not, set ax equal to axis
	for the x-dimension"""
	
	data20 = np.roll(data,-8,axis=2)
	lon20 = np.roll(lon,-8)
	lon20 = np.where(lon20 < 20, lon20 + 360., lon20)
	return data20,lon20

def tgrid(start_time,end_time=None,npts=None,intvl='mth',edges=True,retedges=True):
	if (end_time is None) & (npts is None):
		raise MCMathError("must define either 'end_time' or 'npts'")

	ntim = None
	if npts is not None:
		ntim = npts
	elif end_time is not None:
		if 'day' not in intvl:
			num_yrs = end_time.year - (start_time.year + 1)
			if start_time.month > end_time.month:
				num_mths = (12 - start_time.month) + (end_time.month - 1)
			else:
				num_mths = 11 + (end_time.month - start_time.month)
			if end_time.day >= start_time.day:
				num_mths = num_mths + 1
			if 'mth' in intvl:
				ntim = num_yrs*12 + num_mths
			elif 'seas' in intvl:
				ntim = num_yrs*4 + (num_mths / 3)
			elif ('yr' in intvl) | ('year' in intvl):
				ntim = num_yrs + num_mths / 12
		else:
			ntim = (end_time - start_time).days

	tim = np.empty(ntim,dtype=object)
	if edges is True:
		timedges = np.empty((ntim,2),dtype=object)

	for i in xrange(ntim):
		strt_yr = start_time.year
		strt_mth = start_time.month
		strt_day = start_time.day
		if 'mth' in intvl:
			mth_add = i + (strt_mth - 1)
			yr = mth_add / 12
			mth = (mth_add % 12) + 1
			if mth == 12:
				timedges[i,0] = dt(strt_yr+yr,mth,strt_day)
				timedges[i,1] = dt(strt_yr+yr+1,1,strt_day)
			else:
				timedges[i,0] = dt(strt_yr+yr,mth,strt_day)
				timedges[i,1] = dt(strt_yr+yr,mth+1,strt_day)
			tim[i] = (timedges[i,1] - timedges[i,0])/2 + timedges[i,0]

		if 'seas' in intvl:
			mth_add = 3*i + (strt_mth - 1)
			yr = mth_add / 12
			mth = (mth_add % 12) + 1
			if mth == 10:
				timedges[i,0] = dt(strt_yr+yr,mth,strt_day)
				timedges[i,1] = dt(strt_yr+yr+1,1,strt_day)
			else:
				timedges[i,0] = dt(strt_yr+yr,mth,strt_day)
				timedges[i,1] = dt(strt_yr+yr,mth+3,strt_day)
			tim[i] = (timedges[i,1] - timedges[i,0])/2 + timedges[i,0]

		if 'yr' in intvl:
			yr = i
			timedges[i,0] = dt(strt_yr+yr,1,1)
			timedges[i,1] = dt(strt_yr+yr+1,1,1)
			tim[i] = timedges[i,0] + (timedges[i,1] - timedges[i,0])/2



	return [tim,timedges]


def gcdist2(p1,p2):
	"""Calculates great circle distance between two points, defined by lat and lon values
	p1=(x1,y1) and p2=(x2,y2), on the surface of a sphere having the same volume as
	the Earth."""

	from scipy.special import cosdg,sindg
	
	x1,y1 = p1
	x2,y2 = p2
	d = np.arccos(cosdg(y1)*cosdg(y2)*cosdg(x1-x2) + sindg(y1)*sindg(y2))
	gcdist = d*R_earth

	return gcdist


def gcdist(p1,p2):
	"""Calculates great circle distance between two points, defined by lat and lon values
	p1=(x1,y1) and p2=(x2,y2), on the surface of a sphere having the same volume as
	the Earth.  Uses an approach based on the haversine formula."""

	from scipy.special import cosdg,sindg
	
	x1,y1 = p1
	x2,y2 = p2
	d = 2. * np.arcsin(min(1.,np.sqrt(sindg((y2-y1)/2.)**2 + cosdg(y1)*cosdg(y2)*sindg((x1-x2)/2.)**2)))
	gcdist = d*R_earth

	return gcdist

def OLSAR1(y):
	Nobs = y.size

	ave1 = y[1:].mean()
	ave2 = y[0:-1].mean()

	sum_nom = ((y[1:] - ave1)*(y[0:-1] - ave2)).sum()
	sum_denom = ((y[0:-1] - ave2)**2).sum()

	if sum_denom > 0.:
		OLSAR1 = sum_nom / sum_denom
	else:
		OLSAR1 = 0.

	return OLSAR1

def alph_est(x,RNest,Nsub):
	"""From the VBA scripts by S.N. Rodionov
	Calculates the bias-corrected estimate of alpha for AR(1)
	process like the given time series"""

	if (Nsub < 5) & (RNest == 'MPK'):
		Nsub = 5
		print "subsample size too small for MPK; changed to 5"

	if (Nsub < 3) & (RNest is 'IPN4'):
		Nsub = 3
		print "subsample size too small for IPN4; changed to 3"


	if Nsub > x.size:
		Nsub = x.size

	# subsampling
	N_ct = x.size - Nsub + 1
	ss = np.ma.masked_all((N_ct,))
	for i in xrange(N_ct):
		ss[i] = OLSAR1(x[i:i+Nsub])

	est1 = np.median(ss)
	if RNest is 'MPK':
		alpha = MPK(est1,Nsub)
	elif RNest is 'IPN4':
		alpha = ipn4(est1,Nsub)
	else:
		alpha = est1

	if alpha < 0.0:
		alpha = 0.0
	if alpha > 0.99:
		alpha = 0.99

	return alpha

def ipn4(est,Nsub):
	est = np.float64(est)
	Nsub = np.float64(Nsub)
	ipn4 = est + 1 / Nsub
	for i in xrange(3):
		ipn4 = ipn4 + np.abs(ipn4) / Nsub

	return ipn4

def MPK(est,Nsub):
	est = np.float64(est)
	Nsub = np.float64(Nsub)
	if Nsub > 4:
		mpk = ((Nsub - 1.)*est + 1.) / (Nsub - 4.)
	else:
		# should not be here due to check in alph_est
		mpk = est

	return mpk

def eqn(ns,alpha):
	
	my_sum = 0.
	for i in xrange(1,ns):
		i = np.float64(i)
		ns = np.float64(ns)
		my_sum = my_sum + (1 - (i / ns))*(alpha**i)

	eqn = ns / (1.0 + my_sum)

	if eqn <= 2.: eqn = 2.
	if eqn > ns: eqn = ns

	return eqn


##
#### BEGIN corrtest

def corrtest(x1,x2,mcintts=None,npintts1=None,pval=0.05):
	"""Tests the correlation of two (time) series of equal length.
	Input: x1,x2: the two time series to test correlation
	       mcintts: None or a value of the integral time scale estimate from mcmath.acorr_mskd()
		   npintts: None or a value of the integral time scale estimate from pylab.acorr()
		            if None is given for these, they are calculated within the function
           pval: probability value to test against (default: 0.05)

    Output: R, Rsquared,
	        (CLtests - these give the value of Rsquared needed to pass at the pval sig level)
			IPN,
			OLS,
			MCIntTS,
			NPIntTS
    """
	
	if x1.size != x2.size:
		raise MCMathError("size of inputs must be equal and 1-d")

	x1 = x1 - x1.mean()
	x2 = x2 - x2.mean()
	R = np.corrcoef(x1,x2)[0,1]
	R2 = R**2
	szm1 = x1.size - 1
	if mcintts is None:
		mcintts = max(acorr_mskd(x1,ret_intts=True)[1],acorr_mskd(x2,ret_intts=True)[1])
	if npintts1 is None:
		npintts1 = plt.acorr(x1,normed=True,maxlags=szm1)[1][szm1:][0: np.where(plt.acorr(x1,normed=True,maxlags=szm1)[1][szm1:] < 0)[0][0]].sum()
		npintts2 = plt.acorr(x2,normed=True,maxlags=szm1)[1][szm1:][0: np.where(plt.acorr(x2,normed=True,maxlags=szm1)[1][szm1:] < 0)[0][0]].sum()
		npintts = max(npintts1,npintts2)

	ipn1 = alph_est(x1,'IPN4',x1.size)
	ols1 = alph_est(x1,None,x1.size)
	ipn2 = alph_est(x2,'IPN4',x2.size)
	ols2 = alph_est(x2,None,x2.size)
	ipnm = np.max([ipn1,ipn2])
	olsm = np.max([ols1,ols2])

	ipn_dof = eqn(x1.size,ipnm)
	ols_dof = eqn(x1.size,olsm)
	mcits_dof = x1.size / mcintts
	npits_dof = x1.size / npintts
	#print "ipn2: %.4f :: ols2: %.4f :: ipn: %.4f :: ols: %.4f" % (ipn2,ols2,ipnm,olsm)
	#print "ipn_dof: %.3f :: ols_dof: %.3f" % (ipn_dof,ols_dof)

	cltest_ipn = sp.stats.t.isf(pval,ipn_dof)**2 / (sp.stats.t.isf(pval,ipn_dof)**2 + (ipn_dof - 2.))
	cltest_ols = sp.stats.t.isf(pval,ols_dof)**2 / (sp.stats.t.isf(pval,ols_dof)**2 + (ols_dof - 2.))
	cltest_mcits = sp.stats.t.isf(pval,mcits_dof)**2 / (sp.stats.t.isf(pval,mcits_dof)**2 + (mcits_dof - 2.))
	cltest_npits = sp.stats.t.isf(pval,npits_dof)**2 / (sp.stats.t.isf(pval,npits_dof)**2 + (npits_dof - 2.))

	return [R,R2,cltest_ipn,cltest_ols,cltest_mcits,cltest_npits]

#### END corrtest
##

##
#### BEGIN lagcorr

def lagcorr(x1,x2):
	"""Calculates the lagged correlation between two (time) series, via the
	sp.signal.correlate() function (which 'slides' along the time series
	similar to 'convolve')."""
	x1 = x1 - x1.mean()
	x2 = x2 - x2.mean()
	corr1 = sp.signal.correlate(x1,x2)
	corr1 = (corr1 / (x1.std() * x2.std())) / x1.size

	return corr1

#### END lagcorr
##



def alph_chk(x):
	#print "Ns:\t MPK    \t IPN4   \t OLS"
	N = x.size
	mpks = np.zeros(N)
	ipns = np.zeros(N)
	olss = np.zeros(N)
	for i in xrange(5,x.size):
		#print "%d:\t %.2f   \t %.2f   \t %.2f" % (i,alph_est(x,RNest,i),alph_est(x,'IPN4',i),alph_est(x,1,i))
		mpks[i] = alph_est(x,'MPK',i)
		ipns[i] = alph_est(x,'IPN4',i)
		olss[i] = alph_est(x,1,i)
	mpk_mode = sp.stats.mode(mpks[5:])
	ipn_mode = sp.stats.mode(ipns[5:])
	ols_mode = sp.stats.mode(olss[5:])

	return [mpk_mode,ipn_mode,ols_mode]

def mth_ann_mean(x,dtin=None):
	"""monthly means to annual means, blindly with no accounting for leap years.

	dtin - not operational yet; will allow for true weighting depending on given time axis."""
	N = x.shape
	N12 = N[0] / 12
	outshp = np.array(N)
	outshp[0] = N12
	xout = np.ma.zeros(outshp)

	mths = np.array([31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])
	mthshp = np.ones(outshp.shape)
	mthshp[0] = 12
	mths = mths.reshape(mthshp)

	for i in xrange(N12):
		yrx = x[i*12:(i+1)*12]
		xout[i] = (yrx*(mths/mths.sum())).sum(axis=0)

	return xout

def d2n(date_in,units=None,cal='proleptic_gregorian'):
	if units is None:
		numout = nd2n(date_in,'days since 0001-01-01 00:00:00',cal)
	else:
		numout = nd2n(date_in,units,cal)

	return numout

def n2d(num_in,units=None,cal='proleptic_gregorian'):
	if units is None:
		dateout = nn2d(num_in,'days since 0001-01-01 00:00:00',cal)
	else:
		dateout = nn2d(num_in,units,cal)

	if cal != 'proleptic_gregorian':
		if hasattr(dateout,'size'):
			datenew = np.ndarray(shape=dateout.shape,dtype=object)
			for i in xrange(datenew.size):
				datenew[i] = dt(dateout[i].year,dateout[i].month,dateout[i].day,dateout[i].hour,
								dateout[i].minute,dateout[i].second)
		else:
			datenew = dt(dateout.year,dateout.month,dateout.day,dateout.hour,
								dateout.minute,dateout.second)
		dateout = datenew
		
	return dateout

def redspec(rho,h,M):
	R = (1 - rho**2) / (1 - 2*rho*np.cos(h*np.pi/M) + rho**2)
	return R

def eff(lats,mult=None):
	"""Function for calculating f-values (Coriolis parameter values)
	at each of the latitudes given in 'lats'."""
	omega = 2*np.pi/86400
	if mult is not None:
		y2 = np.tile(lats.reshape(-1,1),(1,mult))
	else:
		y2 = lats
	
	f1 = sp.special.sindg(y2) * 2 * omega

	return f1

#############
##
##  Code for Millenium project files
##
#############

def tstmp_dt(txx):
	"""timestamp in the form of yyyymmdd.hh turned into datetime instances of months, centered
       (instead of default month's end dates)."""
	N = txx.size
	dtxx = np.empty(N,dtype=object)
	dtbnds = np.empty((N,2),dtype=object)
	for i in xrange(len(txx)):
		yr = int(str(int(txx[i]))[:-4])
		mth = int(str(int(txx[i]))[-4:-2])
		day = int(str(int(txx[i]))[-2:])
		d1 = dt(yr,mth,1)
		if mth == 12:
			yr = yr + 1
			d2 = dt(yr,1,1)
		else:
			d2 = dt(yr,mth+1,1)

		dtxx[i] = d1 + (d2 - d1)/2
		dtbnds[i,0] = d1
		dtbnds[i,1] = d2

	return dtxx,dtbnds


def grid_crnr(xcrnr):
	"""Takes the MPI grid corner coord_vars and returns a (n+1)x(m+1) array
	that pcolor likes for plotting 2-d grids"""
	ysize,xsize,blah = xcrnr.shape
	newcrnr = np.zeros((ysize+1,xsize+1))
	for i in xrange(ysize):
		for j in xrange(xsize):
			newcrnr[i,j] = xcrnr[i,j,0]
			if i == (ysize-1):
				newcrnr[i+1,j] = xcrnr[i,j,1]
				if j == (xsize-1):
					newcrnr[i+1,j+1] = xcrnr[i,j,2]
					newcrnr[i,j+1] = xcrnr[i,j,3]
			elif j == (xsize-1):
				newcrnr[i,j+1] = xcrnr[i,j,3]

	return newcrnr

def mpi_gridsize(basin=None):
	from ferr import use
	dgdist = use('/scratch/local2/u241180/data/mill/MPI_grid-distance-in-meter-pressure-points_ocean-grid.nc',silent=1)
	distx = dgdist.v['grid_x_distance_at_pressure_point'][:]
	disty = dgdist.v['grid_y_distance_at_pressure_point'][:]
	gdist = distx*disty

	dbas = use('/scratch/local2/u241180/data/mill/MPI_basin-index_ocean-grid.nc',silent=1)
	bas = dbas.v['basinindex'][0]
	if basin is None:
		bas1 = np.ma.where(bas > 0, 1,np.ma.masked)
	else:
		bas1 = np.ma.where(bas == basin, 1,np.ma.masked)

	gdist = bas1*gdist

	return gdist


##
#### BEGIN depth_sum

def depth_sum(fset,var1,k=None,zind=1):
	"""Function sums up values to a certain depth index 'k', defaults to whole water column.
	For variables in files 'ready' to be integrated by simple addition (e.g., thermosteric grid box
	contributions)
	   Parameters:  fset = globbable string location of files to process
	                var1 = variable name in fset files to process
                    k    = index of maximum depth (index, not depth in meters)
					zind = index of depth dimension of array"""

	from ferr import use
	from glob import glob
	blah = glob(fset)
	list1 = np.sort(blah)
	sizes = []

	if k is None:
		k = slice(None)
	else:
		k = slice(None,k+1)
	for f1 in xrange(list1.size):
		d1 = use(list1[f1],silent=1)
		sizes.append(d1.v[var1].shape[0])
		xshp = d1.v[var1][:].shape
		d1.f.close()
		del d1

	sizes = np.array(sizes,dtype=np.int16)
	xshp = np.array(xshp,dtype=np.int16)
	xshp = np.delete(xshp,zind)
	totsiz = sizes.sum()

	xshp[0] = totsiz
	vout = np.ma.masked_all(shape=xshp,dtype=np.float32)
	tout = np.zeros(totsiz,dtype=np.float64)
	start = 0

	for f1 in xrange(list1.size):
		d1 = use(list1[f1],silent=1)
		t1 = d1.getax(var1,'t')
		v1 = d1.v[var1][:,k,:,:]
		v1_fill = d1.v[var1][:,k,:,:].get_fill_value()

		tout[start:start+sizes[f1]] = t1.copy()
		vout[start:start+sizes[f1]] = v1.sum(axis=zind)

		vout = np.ma.masked_values(vout,v1_fill)

		start = start + sizes[f1]
		del v1,t1
		d1.f.close()
		del d1

	del sizes,xshp,totsiz,start


	return (tout,vout)

#### END depth_sum
##


#############
##
##  Code for other projects
##
#############


##
#### BEGIN psmsl_interp

def psmsl_interp(file1,gaps=1):
	"""Gets PSMSL data from yearly PSMSL files,
	filling data gaps no longer than 'gaps' in length via linear interpolation.
	(Note: 'gaps' kw not functional yet)"""
	
	import mcread

	yr1,msl = mcread.psmsl(file1)

    ##msl = file1

	holes = np.arange(msl.size)[msl.mask]
	cut = np.empty(0,dtype=np.int16)
	for i in xrange(holes.size):
		if (i < holes.size - 1):
			if ((holes[i]+1) == holes[i+1]):
				cut = np.append(cut,i)
		if i > 0:
			if ((holes[i] - 1) == holes[i-1]):
				cut = np.append(cut,i)

	hol1 = np.delete(holes,np.unique(cut))
	ind1 = np.arange(msl.size)
	ind2 = np.delete(ind1,hol1)

	msl_interp = np.interp(ind1,ind2,msl[ind2])
	msl_interp = np.ma.masked_values(msl_interp,-99999.)

	##return holes,cut,hol1,ind1,ind2,msl_interp
	return yr1,msl_interp

#### END psmsl_interp
##


##
#### BEGIN std_analy

def std_analy(data,t=None,fac=None):
	"""Takes monthly input data and calculates the standard deviation,
	the anomaly standard deviation, the 10-year running mean anomaly standard
	deviation, and the trends.  Input data must be monthly and in (T,Y,X)
	axis-order.

	Input: data - monthly input data set
               t - time axis
               fac - factor to multiply 'data' by to obtain the units desired

	Output: [ostd,anomstd,10ystd,trends]
	"""
	if t is None:
		if hasattr(data,'t'):
			t = data.t
		else:
			raise MCMathError("no time axis specified for trend calculation")

	if fac is not None:
		data = fac * data
	data1 = my_dtrnd(data,t)

	hclim = make_clim(data1)
	hanom = anom_from_clim(data1,hclim)

	ann1 = mth_ann_mean(hanom)
	dt1 = t[6::12]
	ann2 = mth_ann_mean(data)

	ostd = data1.std(axis=0,ddof=1)
	anomstd = ann1.std(axis=0,ddof=1)
	h10 = rm_calc(ann1,'h',11)
	ystd10 = h10[1].std(axis=0)
	trends = trnd_map_calc(dt1,ann2)

	return [ostd,anomstd,ystd10,trends]

#### END std_analy
##

