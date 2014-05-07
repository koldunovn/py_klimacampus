## script is meant to be run line-by-line,
## and not using execfile() or the ipython run magic function

## a short example

m = ones((2,5))

p = arange(5.)

q = arange(2.)

# works (last dim of m is same as dim of p)
m*p

# doesn't work (last dim of m is different than dim of q)
m*q

# works (q is reshaped to have the same number of dims as m, with the dim the same
# for the axis over which the multiplication is to take place)
m*q.reshape(2,1)


d1 = use('/data/icdc/ocean/ishii/DATA/temp.2004.nc')
temp = d1.gv('WTMP_GDS0_DBSL')


zintvl = mcmath.intvls(temp.z)

zintvl



mths = np.array([31.,29.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.])


temp * mths

t12 = temp * mths.reshape(12,1,1,1)

tann = t12.sum(axis=0) / mths.sum()

tdep1 = temp * zintvl.reshape(1,-1,1,1)

tdep = temp.sum(axis=1) / zintvl.sum()

## when DOESN'T this work?
## when one must apply a function to a portion of the array,
## and the function doesn't support multi-dimensional arrays


# let's do a lo-pass filter on the entire data set for the time series at each x,y point

d2 = use('/data/icdc/ocean/ishii/DATA/temp.2005.nc')
temp2 = d2.gv('WTMP_GDS0_DBSL')


temp1_0 = temp[:,0,...]  # surface slice, z=0
temp2_0 = temp2[:,0,...]

t0 = r_[temp.t,temp2.t]
temp0 = r_[temp1_0, temp2_0]

temp0 = np.ma.masked_greater(temp0,1e19)


from scipy.signal import butter, freqz, filtfilt

b,a=butter(6,0.25);
h,w=freqz(b,a,128);

lp1out = np.ma.masked_all(temp0.shape,dtype=np.float32)

for i in np.ndindex(temp0.shape[1:]):
        if temp0[:,i[0],i[1]].mask.any(): continue
        lp1out[:,i[0],i[1]] = filtfilt(b,a,temp0[:,i[0],i[1]])


plot(t0,temp0[:,90,170])

plot(t0,lp1out[:,90,170],'r--',lw=2)

