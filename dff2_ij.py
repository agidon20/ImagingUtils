# Jython imports
from ij import IJ, ImagePlus, WindowManager
from ij.plugin.frame import RoiManager
from ij.plugin import Duplicator
import os
from java.util import ArrayList
from java.util import Arrays
from java.lang import Math
from glob import glob
from ij.gui import Plot

# Import Apache Commons Math classes
from org.apache.commons.math3.fitting import CurveFitter
from org.apache.commons.math3.fitting import WeightedObservedPoints
from org.apache.commons.math3.analysis import ParametricUnivariateFunction
from org.apache.commons.math3.optim.nonlinear.vector.jacobian import LevenbergMarquardtOptimizer



global roiman, ij


def initialize_imagej():
    global roiman, ij
    ij = IJ.getInstance()
    print("ImageJ2 version: ", ij.getInfo())
    roiman = RoiManager.getRoiManager()

def duplicate_imageplus(fname):
	print fname
	tmp = IJ.openImage(fname)
	imp = Duplicator().run(tmp)
	return imp

def remove_slices(imp, rng):
    stack = imp.getStack()
    imp.setT(stack.size())  # get fluorescence of the last slice.
    for r in sorted(rng, reverse=True):
        stack.deleteSlice(r)  # remove from the end to not change indexes
        imp.setStack(stack)  # without this line deleteSlice will not work properly

def imp_open(rng):
	fname = glob('*.tif')[0]
	imp = duplicate_imageplus(os.getcwd() + "\\" + fname)  # do not work on the original file. too risky
	remove_slices(imp, rng)  # remove all the slices that are noise (mostly before and after recordings)
	return imp

def mean(imp, label):
    roi = roiman.getROIs()[label]
    imp.setRoi(roi)
    stats = imp.getStatistics()
    return stats.mean



class ExpFunc(ParametricUnivariateFunction):
    def value(self, x, params):
        a, b, c, k, l = params
        return a * Math.exp(-x / b) + k * Math.exp(-x / l) + c

    def gradient(self, x, params):
        a, b, c, k, l = params
        da = Math.exp(-x / b)
        db = a * x / (b * b) * Math.exp(-x / b)
        dc = 1.0
        dk = Math.exp(-x / l)
        dl = k * x / (l * l) * Math.exp(-x / l)
        return [da, db, dc, dk, dl]

def fit_f0(stim_start, t, f):
    xfit = ArrayList()
    for i in range(stim_start):
        xfit.add(t[i])
    for i in range(len(t)-3, len(t)-1):
        xfit.add(t[i])
    
    yfit = ArrayList()
    for i in range(stim_start):
        yfit.add(f[i])
    for i in range(len(f)-3, len(f)-1):
        yfit.add(f[i])

    optimizer = LevenbergMarquardtOptimizer()
    fitter = CurveFitter(optimizer)
    for i in range(len(xfit)):
        fitter.addObservedPoint(xfit[i], yfit[i])

    initial_guess = [1.0, 1.0, 1.0, 1.0, 1.0]
    params = fitter.fit(ExpFunc(), initial_guess)

    print("fit: a=%5.3f, b=%5.3f, c=%5.3f, k=%5.3f, l=%5.3f" % tuple(params))
    return [ExpFunc().value(x, params) for x in t]

    
    
def stack_mean_roi(imp):
    t = [i + 1 for i in range(imp.getStack().size())]
    f, bg = ArrayList(), ArrayList()
    for slice in t:
        imp.setSlice(slice)
        f.add(mean(imp, "cell"))
        bg.add(mean(imp, "bg"))
    return t, f, bg

def dff(exp_dir, cell, repetitions, slices_to_remove, stim_start, folder, folderRoi):
    t, f, bg = ArrayList(), ArrayList(), ArrayList()
    for id in repetitions:
        os.chdir(exp_dir + cell + "_" + str(id))
        if folderRoi == ".":
            roiman.runCommand("Open", "./RoiSet.zip")
        else:
            roiman.runCommand("Open", exp_dir + folderRoi + "/RoiSet.zip")
        imp = imp_open(slices_to_remove)
        t, f1, bg1 = stack_mean_roi(imp)
        if f.isEmpty():
            f = f1
        else:
            for i in range(len(f)):
                f.set(i, f.get(i) + f1.get(i) / len(repetitions))
        if bg.isEmpty():
            bg = bg1
        else:
            for i in range(len(bg)):
                bg.set(i, bg.get(i) + bg1.get(i) / len(repetitions))
        imp.close()
        roiman.runCommand("Delete")
    f0 = fit_f0(stim_start, t, f)  # f0 is the fitted version to remove bleaching
    return t, [(f.get(i) - f0[i]) / (f0[i] - bg.get(i)) for i in range(len(f))], f, bg, f0

def save_data(dir1, cell, t,dff1):
	os.chdir(dir1)
	filename = os.path.join(dir1, cell + ".er.txt")
	with open(filename, "w") as f:
		for i in range(len(t)):
			f.write( str(t[i]) + " " + str(dff1[i]) + "\n")
	print("Data saved to "+ filename)

            
def plot_data(t, dff1, f, bg, f0):
	
	t_array = [t[i] for i in range(len(t))]
	f_array = [f[i] for i in range(len(f))]
	dff_array = [dff1[i] for i in range(len(dff1))]    
	bg_array = [bg[i] for i in range(len(bg))]
	f0_array = [f0[i] for i in range(len(f0))]

	plot = Plot("Fluorescence Data (parameters)", "Time", "Intensity")
	plot.setColor("black")
	plot.add("line", t_array, f_array)
	plot.setColor("blue")	
	plot.add("line", t_array, bg_array)
	plot.setColor("red")
	plot.add("line", t_array, f0_array)
	plotdff = Plot("Fluorescence Data (dff)", "Time", "Intensity")	
	plotdff.setColor("orange")	
	plotdff.add("line", t_array, dff_array)    
	plot.show()
	plotdff.show()	



ephys_dir = "C:/port/hudrive/OneDrive - hu-berlin.de/lab/sfv/ephys/"
exp_dirs = [
    ephys_dir + "sfv190321/Sina_SFV_Ccamp/",
    ephys_dir + "sfv190516/imaging/"
]
imagej_dir = "C:/port/hudrive/OneDrive - hu-berlin.de/lab/sfv/imaging/analisys/Fiji.app"
plugins_dir = "C:/port/hudrive/OneDrive - hu-berlin.de/lab/sfv/imaging/analisys/Fiji.app/plugins"
save_data_dir = "C:/port/hudrive/OneDrive - hu-berlin.de/lab/sfv/imaging/withephys"

# [cellid, repetitions, slices to remove, stimulus first time, folder, folderroi]
cell = (ephys_dir + "20231122/imaging/", "C4S14", [1, 2, 3, 4, 5], list([1,2]), 15, "231122", "C4S12_1")


initialize_imagej()
t, dff1, f, bg, f0 = dff(*cell)
plot_data(t,dff1, f, bg, f0)
save_data(".", cell[5] + cell[1], t, dff1)
