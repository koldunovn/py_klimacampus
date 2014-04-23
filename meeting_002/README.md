Meeting #2
====
26.03.2014

Presentations:

* IPython tips and tricks and mcocean module

Mark's modules
==============

In order to get `sst.v3b.mnmean.nc` use [this link](https://drive.google.com/file/d/0B0WRA0YFFl_tb3lrZWhZWk44RGs/edit?usp=sharing).

Note: all of my modules use `Numpy`, `Scipy`, `Python NetCDF4`, and `Matplotlib` modules.

* **ferr.py**:  gives a simple way to open and access data in NetCDF files formatted using CF conventions.  It has some issues still with curvilinear coordinates and time axes for climatologies, but is under revision.  Also has a function for saving simple NetCDF files.

* **mcmath.py**: Math functions of various levels of usefulness.  My most used functions within it include: 
  * **trnd_map_calc** (for calculating the trends of time series over an entire map of time series; i.e., t,y,x data),
  * **area_cswt_mean** (for calculating the planview area of a dataset; allows for global spatial means to be calculated leaving other dimensions intact, like Ferret can),
  * **my_dtrnd** (can detrend time series with masked values in it), 
  * **mth_ann_mean** [can calculate the annual mean from monthly data with (mostly) correct weighting],
  * **num2date (n2d)**, **date2num (d2n)**: the ones from Python NetCDF4 module are MUCH better than the ones from matplotlib, especially if interacting with files back and forth between matplotlib and Ferret.

* **mcp.py**: simplifies some of the most-common everyday plots
  * **dateplt1**, **dateplt_m**, **dateplt_m2** – allows well-formatted plots with year dates along the bottom (I'm a long-time-scale kinda scientist), the _m / _m2 variants allow for adding time series to the same plot, with the plot scale expanding to cover any range outside the old plot range
  * **shade** – an attempt to replace Ferret's 'shade' function, so that making planview maps are easy. Works pretty well, but not on curvilinear data
  * **cvshade** – works on curvilinear data, but only if the grid is well defined in the NetCDF files, and you have to pass this function the corners of the grid cells
  
Examples with Mark's modules
============================

* [Use of `ferr` and `shade`](http://nbviewer.ipython.org/github/koldunovn/py_klimacampus/blob/master/meeting_002/Tutorial.ipynb?create=1)

* [Couple of 1-d plot examples](http://nbviewer.ipython.org/github/koldunovn/py_klimacampus/blob/master/meeting_002/mcp_1d_plots.ipynb?create=1)





  
