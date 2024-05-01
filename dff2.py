# https://imagej.nih.gov/ij/developer/macro/functions.html#getMeasurementValue
# https://github.com/imagej/pyimagej/blob/main/doc/Initialization.md
# https://pyimagej.readthedocs.io/en/latest/index.html
#  set the auto contrast option on - just for visibility
#  ij.IJ.run("Appearance...", "auto menu=15 gui=1 16-bit=Automatic")


import imagej
import scyjava
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as mpl
from scipy.optimize import curve_fit
from glob import glob
import sqlite3
import pandas as pd

#matplotlib.use(backend='TkAgg', force=True)
matplotlib.use(backend='Qt5Agg', force=True)

global roiman, ij
sfvdir = "C:/port/hudrive/OneDrive - hu-berlin.de/lab/sfv"
ephys_dir = f"{sfvdir}/ephys/"
imagej_dir = f"{sfvdir}/imaging/analisys/Fiji.app"
plugins_dir = f"{sfvdir}/imaging/analisys/Fiji.app/plugins"
save_data_dir = f"{sfvdir}/imaging/withephys"

def initialize_imagej():
    global roiman, ij
    scyjava.config.add_option('-Xmx6g')
    scyjava.config.add_option(f'-Dplugins.dir={plugins_dir}')
    ij = imagej.init(imagej_dir, mode=imagej.Mode.INTERACTIVE)
    print(f"ImageJ2 version: {ij.getVersion()}")
    roiman = ij.RoiManager.getRoiManager()


def duplicate_imageplus(fname):
    tmp = ij.IJ.openImage(fname)
    imp = tmp.duplicate()  # do not work on the original file. too risky
    tmp.close()
    return imp


def remove_slices(imp, slices_to_remove_left,slices_to_remove_right):
    stack = imp.getStack()
    imp.setT(stack.size())  # get fluorescence of the last slice.
    rng = list(range(1,slices_to_remove_left+1))
    if(slices_to_remove_right > -99):
         rng = rng + list(range(slices_to_remove_right,stack.size()+1))
    for r in np.sort(rng)[::-1]:
        stack.deleteSlice(r)  # remove from the end to not change indexes
        imp.setStack(stack)  # without this line deleteSlice will not work properly

    # rng = list(range(1,slices_to_remove_left))
    # if(slices_to_remove_right > -99):
    #      rng = rng + list(range(slices_to_remove_right,130+1))
    # for r in np.sort(rng)[::-1]:
    #     print(r)


def imp_open(slices_to_remove_left,slices_to_remove_right):
    fname = glob('*.tif')[0]
    imp = duplicate_imageplus(fname)  # do not work on teh original file. too risky
    remove_slices(imp, slices_to_remove_left,slices_to_remove_right)  # remove all the slices that are noise (mostly before and afer recordings)
    return imp


def mean(imp, label):
    roi = roiman.getROIs()[label]
    imp.setRoi(roi)
    return float(imp.getRawStatistics().mean)


def fit_f0(stim_start, t, f):
    def expfunc(x, a, b, c, k, l):
        return a * np.exp(-x / b) + k * np.exp(-x / l) + c

    xfit = np.append(t[0:stim_start], t[-3:-1])
    yfit = np.append(f[0:stim_start], f[-3:-1])

    # xfit = t[1:stim_start]
    # yfit = f[1:stim_start]

    popt, pcov = curve_fit(expfunc, xfit, yfit,
                           bounds=(0.0, [1.2 * max(yfit), 10 * max(t), 1.2 * max(yfit), 1.2 * max(yfit), 10 * max(t)]))
    print('fit: a=%5.3f, b=%5.3f, c=%5.3f, k=%5.3f, l=%5.3f' % tuple(popt))
    return expfunc(t, *popt)


def stack_mean_roi(imp):
    t = np.arange(1, imp.getStack().size() + 1)
    f, bg = np.array([]), np.array([])
    for slice in t:
        imp.setSlice(int(slice))
        f = np.append(f, mean(imp, "cell"))
        bg = np.append(bg, mean(imp, "bg"))
    return t, f, bg


# a session is something that can be averaged somethings like
#  once cell with the same number of repetitions and stimulus etc'
def dff(exp_dir, cell, repetitions, slices_to_remove_left, slices_to_remove_right, stim_start,folderRoi):
    t, f, bg = 0, 0, 0
    for id in repetitions:
        os.chdir(f"{exp_dir}{cell}_{id}")
        if(folderRoi == "."):
            roiman.runCommand("Open", f"./RoiSet.zip")
        else:
            roiman.runCommand("Open", f"{exp_dir}{folderRoi}/RoiSet.zip")
        imp = imp_open(slices_to_remove_left,slices_to_remove_right)
        t, f1, bg1 = stack_mean_roi(imp)
        f = f1 if f is None else f + f1 / len(repetitions)
        bg = bg1 if bg is None else bg + bg1 / len(repetitions)
        imp.close()
        roiman.runCommand("Delete")
    f0 = fit_f0(stim_start, t, f)  # f0 is the fitted version to remove bleaching
    return t, (f - f0) / (f0 - bg), f, bg, f0


def save_data(dir, cell, xy_data):
    os.chdir(dir)
    with open(cell + '.npy', 'wb') as f:
        np.save(f, np.array(xy_data))


#####################################################################################
## Here starts implementation, it is rather simple depending on the way that you store your data.
## I stored the data in sqlite db file. But can be done any other way.
## as long as you have these two lines:
##        arg = (path + file, exclude_left,exclude_right, stimulus_start,ROI_folder)
##        t, dff1, f, bg, f0 = dff(*arg)

## the data is specified here:
## path + file: path + file name for to the image stacks
## exclude_left: the number of images to exclude from the begining of the image stack
## exclude_right: the number of images to exclude at the end of the image stack
## stimulus_start:  what image the stimulus starts
## ROI_folder: folder with the RIO file (name rioset.zip). The RIO file should include three ROIs, cell.roi, bg.roi and all.roi.


sql = "SELECT * FROM DFF WHERE avoid is null and virus = 'AAV' AND cellstate regexp 'C[0-9]S[0-9]*' AND inputrate = 1;"
# Read sqlite query results into a pandas DataFrame
con = sqlite3.connect("C:\\port\\hudrive\\OneDrive - hu-berlin.de\\lab\\sfv\\database\\data\\_sfv_.db")
con.enable_load_extension(True)
con.load_extension("C:\\port\\hudrive\\OneDrive - hu-berlin.de\\lab\\sfv\\database\\tools\\SQLiteStudio\\extensions\\regexp")
df = pd.read_sql_query(sql, con)
# Verify that result of SQL query is stored in the dataframe
con.close()
df = df.reset_index() 

initialize_imagej()
veusz_cmd = ""

## example of a row from the database:
## sfv190826/Sina_AAV_20190826	| C2S9	| 5	| 2	| 81 | 190826 | C2S8_1 | 15
df['path'] = df['path'].astype('string')
df['cellstate'] = df['cellstate'].astype('string')
df['repetitions'] = df['repetitions'].astype('Int64')
df['exclude_left'] = df['exclude_left'].astype('Int64') 
df['exclude_right'] = df['exclude_right'].astype('Int64') 
df['stimulus_start'] = df['stimulus_start'].astype('Int64')
df['ROI_folder'] = df['ROI_folder'].astype('string')

for col in df:
    df[col] = df[col].fillna(-99)
maxdff = np.array([])
for index, cell in df.iterrows():
    meandff = []
    for rep in range(1,cell['repetitions']+1):
        #t, dff1, f, bg, f0 = dff(*cell[0:7])
        arg = (ephys_dir + cell['path'] + '/',cell['cellstate'],f"{rep:d}",cell['exclude_left'],cell['exclude_right'],cell['stimulus_start'],cell['ROI_folder'])
        print (arg)
        t, dff1, f, bg, f0 = dff(*arg)
        if(meandff == []):
            meandff = dff1
        else:
            meandff = meandff + dff1
    #     save_data(save_data_dir, cell[5] + cell[1], [t, dff1])
        # mpl.plot(t, f,'--b')
        # mpl.plot(t, f0,'-g')
        # mpl.plot(t, bg,'.r')
         #mpl.plot(t, dff1)   
    meandff = meandff/cell['repetitions']
    mpl.plot(t, meandff)
    maxdff = np.append(maxdff, np.max(meandff))
#     veusz_cmd += f"ImportFilePlugin(u'Numpy NPY import', u'{cell[5]}{cell[1]}.npy', linked=False, errorsin2d=False, name=u'{cell[5]}{cell[1]}')\n"
print(maxdff)
mpl.show()
print(veusz_cmd)
