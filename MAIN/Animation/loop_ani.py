"""
This will produce an animation of the percentile bootstrap test for education purposes
Input data must be in csv format and two dimensional (rows = trials, columns = time frames)
So, two sheets like this are needed (one for each condition)
Also, a third csv file is needed which has the time vector in milliseconds that corresponds to the time frames

The main function is loop_ani.run which calls most of the other functions and sets their parameters. Feel free to adjust the parameters.
Also, here you can change the name of the time vector file. You will also need to change file names in loop_ani.calc
Note that the function loop_ani.getOpts controls much of the details in the vizualization (labels, ticks, offsets, scaling, etc...)

Copyright (C) 2017 Allan Campopiano

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
from scipy import stats
from matplotlib.patches import ConnectionPatch
import matplotlib
matplotlib.rcParams['toolbar'] = 'None'

def calc(ntrials, nboot, baseline):

    # baseline example baseline = [0, 102]

    # data
    cond1 = pd.read_csv('/home/allan/HCDSB/Research/python_code/condition1.csv', header=None).T.values[0:ntrials,0:600]
    cond2 = pd.read_csv('/home/allan/HCDSB/Research/python_code/condition2.csv', header=None).T.values[0:ntrials,0:600]

    # baseline
    m=cond1[:,baseline[0]:baseline[-1]].mean(axis=1)
    m = m[:, np.newaxis]
    cond1=cond1-m

    m=cond2[:,baseline[0]:baseline[-1]].mean(axis=1)
    m = m[:, np.newaxis]
    cond2=cond2-m

    # bootstrap cond 1
    boot_cond1= np.array([]).reshape(0, cond1.shape[1])
    trials=cond1.shape[0]
    boot_inds1=np.random.randint(trials, size=(nboot, trials)).astype(int)
    for i in range(nboot):
        samp=cond1[boot_inds1[i,:],:]
        tr = stats.trim_mean(samp, .2, 0)
        rem=tr[baseline[0]:baseline[-1]].mean()
        tr=tr-rem
        boot_cond1 = np.vstack([boot_cond1, tr])

    # bootstrap cond 2
    boot_cond2 = np.array([]).reshape(0, cond2.shape[1])
    trials=cond1.shape[0]
    boot_inds2=np.random.randint(trials, size=(nboot, trials)).astype(int)
    for i in range(nboot):
        samp=cond2[boot_inds2[i,:],:]
        tr = stats.trim_mean(samp, .2, 0)
        rem=tr[baseline[0]:baseline[-1]].mean()
        tr=tr-rem
        boot_cond2 = np.vstack([boot_cond2, tr])

    #diff
    diff_array=boot_cond1-boot_cond2


    return boot_inds1, boot_inds2, boot_cond1, boot_cond2, cond1, cond2, diff_array

def init(ax1, cond1, tv, opts):
    # init(ntrials, cond1, fig, ax1, ax2, ax3, tv):

    # colors
    cmap = cm.get_cmap('viridis')
    c = cmap(np.linspace(0, 1, opts['plot1']['nshow']))

    #numSamples, numRows = len(tv), ntrials
    #data = cond1[0:ntrials, :]

    #numSamples, numRows = len(tv), ntrials
    #data = cond1[0:ntrials, :]

    for ii in range(opts['plot1']['nshow']):
        #ii in range(numRows):
        ax1.plot(tv, cond1[ii, :] + opts['plot1']['offsets'][ii], '-', markersize=1.5, lw=1.5, c=c[ii])
        #ax1.plot(tv, data[ii, :] + opts['plot1']['offsets'][ii], '-', markersize=1.5, lw=1.5, c=c[ii])

    plt.pause(.3)

    return c

def resamp(ax1, ax2, opts, ntrials, inds, cond1, tv, c):

    off = np.flipud(opts['plot2']['offsets'])
    boot_data=cond1[inds, :]
    aro_list=[]

    for i in range(opts['plot2']['nshow']):
        #i in range(0, ntrials):

        x = boot_data[i, :] + off[i]
        y = tv
        ax2.plot(y, x, '-', lw=1.5, markersize=1.5, c=c[inds[i]])

        # Create annotation
        xy1 = (y[0], boot_data[i, 0] + off[i])
        xy2 = (y[-1], boot_data[i, 0] + opts['plot2']['offsets'][inds[i]])
        coordsA = "data"
        coordsB = "data"
        ann = ConnectionPatch(xyA=xy1, xyB=xy2, coordsA=coordsA, coordsB=coordsB,
                              axesA=ax2, axesB=ax1,
                              arrowstyle="->", shrinkB=0)
        aro_list.append(ann)
        ax2.add_artist(ann)

        plt.pause(.001)

        for aro in aro_list:
            aro.remove()
        aro_list[:] = []

    return ax1, ax2, ntrials, tv, c

def trim(ax2, ax3, opts, row, boot_data, tv):

    rec_list=[]
    off = np.flipud(opts['plot3']['offsets'])
    tr_data = boot_data

    lin=[]
    step=20
    for i in range(0, len(tv), step):

        for li in lin:
            li.pop(0).remove()
        lin=[]

        l=ax3.plot(tv[0:i+step], tr_data[0:i+step]+off[row], '-', lw=1.5, markersize=1.5, c='dimgray')
        lin.append(l)

        rect = ax2.axvline(tv[i], color='indianred', zorder=10)

        rec_list.append(rect)

        plt.pause(.1)

        for rec in rec_list:
            rec.remove()
        rec_list[:] = []

def resamp_trim(ax1, ax2, ax3, opts, ntrials, cond1, boot_cond1, tv, c):
    #resamp_trim(ax1, ax2, ax3, opts, ntrials, boot_inds1, cond1, boot_cond1, tv, c):

    for i in range(0, 4):

        # only here when using a subset of the data for sake of visualization
        inds = np.random.randint(opts['plot1']['nshow'], size=(1, opts['plot1']['nshow'])).astype(int)[0]

        #ax1, ax2, ntrials, tv, c = resamp(ax1, ax2, opts, ntrials, boot_inds1[i,:], cond1, tv, c)
        ax1, ax2, ntrials, tv, c = resamp(ax1, ax2, opts, ntrials, inds, cond1, tv, c)
        trim(ax2, ax3, opts, i, boot_cond1[i,:], tv)

        # clear data from plot
        for li in ax2.lines + ax2.collections:
            li.remove()

def show_msg(msg_ind, fig, axs, wait, last=False):

    align_opts=['left', 'left', 'center', 'center', 'center']
    pos_opts=[.05, .05, .5, .5, .5]

    msgs=['• The following animation demonstrates how to carry out the percentile bootstrap test using trimmed means\n\n'
          '• This variation of the bootstrap test yields confidence intervals and p-values for a single-subject\n',

          '• Begin by randomly resampling with replacement N trials from the original set of N trials\n\n'
          '• Using this new set of trials, compute the trimmed ERP (which is based on the trimmed mean)\n\n'
          '• Repeat the resampling and averaging B times where B is usually in the thousands\n\n'
          '• For the sake of demonstration, this is only repeated a few times in the following animation',

          'The bootstrapping procedure just shown is repeated again for a second condition (not shown),\n\n'
          'and the differences between conditions are taken as follows',

          'The final animation sorts the differences and then uses them to derive confidence intervals\n',

          'More information about robust statistics as well as Python code for this animation\n'
          'can be found at https://github.com/Alcampopiano/STATSLAB\n'

          ]

    for a in axs:
        a.cla()
        a.axis('off')

    big=fig.add_subplot(111)
    big.axis('off')

    te=plt.text(pos_opts[msg_ind],.5, msgs[msg_ind], ha=align_opts[msg_ind], va='center', visible=True, fontsize=18, color='dimgray')
    plt.pause(wait)
    te.remove()

    if last:
        for a in axs:
            a.axis('off')
    else:
        for a in axs:
            a.axis('on')

def diff(ax1, ax2, ax3, boot_cond1, boot_cond2, tv, opts):

    # b2 = boot_cond2[0:nboot, :]
    # b1 = boot_cond1[0:nboot, :]
    b2 = boot_cond2[range(opts['plot5']['nshow']), :]
    b1 = boot_cond1[range(opts['plot4']['nshow']), :]

    diff_array = b1 - b2

    for ii in range(opts['plot4']['nshow']):
    #for ii in range(nboot):
        ax1.plot(tv, b1[ii, :] + opts['plot4']['offsets'][ii], '-', markersize=1.5, lw=1.5, c='dimgray')
        ax2.plot(tv, b2[ii, :] + opts['plot5']['offsets'][ii], '-', markersize=1.5, lw=1.5, c='dimgray')

    # viz
    step=20
    rec_list1=[]
    rec_list2=[]
    lin=[]
    for i in range(0, len(tv), step):

        for j in range(opts['plot6']['nshow']):
        #for j in range(0, nboot):
            l=ax3.plot(tv[0:i+step], diff_array[j, 0:i+step]+opts['plot6']['offsets'][j], c='dimgray', lw=1.5)
            lin.append(l)

        rect1 = ax1.axvline(tv[i], color='indianred', zorder=10)
        rect2= ax2.axvline(tv[i], color='indianred', zorder=10)
        rec_list1.append(rect1)
        rec_list2.append(rect2)

        plt.pause(.001)

        for li in lin:
            li.pop(0).remove()
        lin=[]

        for rec in zip(rec_list1, rec_list2):
            rec[0].remove()
            rec[1].remove()
        rec_list1[:] = []
        rec_list2[:] = []


    #return diff_array

def set_axes(axs, opts, plot_keys):

    for a, pl in zip(axs, plot_keys):
        #print(opts[pl]['ytick_loc'])
        #print(opts[pl]['ytick_lab'])
        a.set_xlabel(opts[pl]['xlab'])
        a.set_ylabel(opts[pl]['ylab'])
        a.set_title(opts[pl]['title'])
        a.set_xlim(opts[pl]['xlim'])
        a.set_yticks(opts[pl]['ytick_loc'])
        a.set_yticklabels(opts[pl]['ytick_lab'])
        #a.set_xticks(opts[pl]['xticks'])
        #a.axes.get_xaxis().set_ticks([])
        a.set_xticks(opts[pl]['xtick_loc'])
        a.set_xticklabels(opts[pl]['xtick_lab'])
        #print(opts[pl]['ylim'])
        a.set_ylim(opts[pl]['ylim'])
        a.title.set_size(opts[pl]['title_fs'])
        a.xaxis.label.set_size(opts[pl]['xlab_fs'])
        a.yaxis.label.set_size(opts[pl]['ylab_fs'])

        for tick in a.yaxis.get_major_ticks():
            tick.label.set_fontsize(opts[pl]['ytick_fs'])

        for tick in a.xaxis.get_major_ticks():
            tick.label.set_fontsize(opts[pl]['xtick_fs'])

def run(ntrials, nboot):
    #plt.pause(5)
    boot_inds1, boot_inds2, boot_cond1, boot_cond2, cond1, cond2, diff_array=calc(ntrials,nboot, [0, 102])

    tv = pd.read_csv('/home/allan/HCDSB/Research/python_code/time_vector.csv', header=None).values[0][:600]

    opts=getOpts(boot_cond1, boot_cond2, cond1, tv)

    #parent figure to start
    size = [ 18.69,   9.86]
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=size)
    fig.subplots_adjust(right=.85)

    show_msg(0, fig, (ax1, ax2, ax3), 11)
    show_msg(1, fig, (ax1, ax2, ax3), 12)

    # legend
    leg_axs=makeCustomLegend(fig, which_leg=1)

    #1 init
    set_axes((ax1, ax2, ax3), opts, ['plot1', 'plot2', 'plot3'])
    c=init(ax1, cond1, tv, opts)

    #2 resample and trim
    resamp_trim(ax1, ax2, ax3, opts, ntrials, cond1, boot_cond1, tv, c)

    # legend
    for la in leg_axs:
        fig.delaxes(la)

    show_msg(2, fig, (ax1, ax2, ax3), 7)

    # legend
    leg_axs=makeCustomLegend(fig, which_leg=1)

    set_axes((ax1, ax2, ax3), opts, ['plot4', 'plot5', 'plot6'])
    diff(ax1, ax2, ax3, boot_cond1, boot_cond2, tv, opts)

    # legend
    for la in leg_axs:
        fig.delaxes(la)

    show_msg(3, fig, (ax1, ax2, ax3), 7)

    # legend
    leg_axs=makeCustomLegend(fig, which_leg=2)

    set_axes((ax1, ax2, ax3), opts, ['plot7', 'plot8', 'plot9'])
    sort_and_CI(ax1, ax2, ax3, diff_array, nboot, tv, opts)

    plt.pause(3)
    for la in leg_axs:
        fig.delaxes(la)

    show_msg(4, fig, (ax1, ax2, ax3), 7, last=True)

    return boot_cond1, boot_cond2, diff_array, cond1, cond2

def bubblesort(ax, data, nshow, rects):

    step=30
    b = list(range(0, nshow, step))
    b.append(len(data) - 1)

    for i in range(len(data)):
        for k in range(len(data) - 1, i, -1):
            if data[k] < data[k - 1]:
                tmp = data[k]
                data[k] = data[k-1]
                data[k-1] = tmp

        if i in b:
            for rect, h in zip(rects, data):
                rect.set_height(h)
                lims = (min(data), max(data))
                ax.set_ylim(lims[0], lims[1])

            plt.pause(.1)

def sort_and_CI(ax1, ax2, ax3, diff_array, nboot, tv, opts):

    rec_list= []

    for i in range(opts['plot7']['nshow']):
        ax1.plot(tv, diff_array[i,:]+opts['plot6']['offsets'][i], c='dimgray', lw=1.5)

    rects = ax2.bar(range(opts['plot7']['nshow']), diff_array[range(opts['plot7']['nshow']),0],.8, color='darkgray')
    #rects = ax2.bar(range(0, nboot), diff_array[:, 0], .8, color='darkgray')

    # centering becasue there was no trial cap
    # m=diff_array.mean(axis=1)
    # diff_array=diff_array-m[0]

    srt_dif=diff_array.copy()
    srt_dif.sort(axis=0)
    l_real = round(.05 * nboot / 2)
    u_real = nboot-l_real

    l_fake = round(.05 * opts['plot7']['nshow'] / 2)
    u_fake = opts['plot7']['nshow']-l_fake

    step=20
    objs = []
    for i in range(0, len(tv), step):

        for rec in rec_list:
            rec.remove()
        rec_list[:] = []

        for zobj in objs:
            zobj.pop(0).remove()
        objs=[]

        tmp = diff_array[:,i]

        # add rectangle
        rect = ax1.axvline(tv[i], color='indianred', zorder=10)

        # set ax2 label for current time
        ax2.set_xlabel('values at ' + str(int(round(tv[i]))) + ' ms')

        rec_list.append(rect)

        # draw CI
        # low and upper bound
        low_bnd=srt_dif[l_real, 0:i+step]
        up_bnd=srt_dif[u_real, 0:i+step]
        #mid=(low_bnd+up_bnd)/2
        obj1=ax3.plot(tv[0:i+step], low_bnd, '-', c='dimgray')
        obj2=ax3.plot(tv[0:i+step], up_bnd, '-', c='dimgray')
        #obj3 = ax3.plot(tv[0:i + step], mid, '-', c='white')
        fill=ax3.fill_between(tv[0:i+step], low_bnd, up_bnd, color='whitesmoke', alpha=1)
        sig_line = ax3.plot(tv[0:i+step], np.zeros(i+step), color='gold', zorder=-100)

        objs.append(obj1)
        objs.append(obj2)
        #objs.append(obj3)
        objs.append(sig_line)

        # show bubble
        bubblesort(ax2, tmp, opts['plot7']['nshow'], rects)

        rect_low = ax2.axvline(l_fake, color='royalblue', zorder=10)
        rect_up = ax2.axvline(u_fake, color='royalblue', zorder=10)

        # show circles tv[bgn + step], low_bnd
        c1 = ax3.scatter(tv[i + step-1], low_bnd[-1], s=100, c='royalblue', alpha=1, zorder=10)
        c2 = ax3.scatter(tv[i + step-1], up_bnd[-1], s=100, c='royalblue', alpha=1, zorder=10)

        plt.pause(.3)

        rec_list.append(rect_low)
        rec_list.append(rect_up)
        rec_list.append(c1)
        rec_list.append(c2)
        rec_list.append(fill)

def makeCustomLegend(fig, which_leg=1):

    leg_axs=[]

    if which_leg==1:
        ax4 = fig.add_axes((.86, .71, .05, .05))
        rect = patches.Rectangle((.5, 0), .01, 1, edgecolor='indianred', facecolor='indianred', alpha=1, zorder=10)
        ax4.add_patch(rect)
        ax4.text(0.7, 0.4, 'current time', color='k', size=15)
        ax4.axis('off')
        leg_axs.append(ax4)

    elif which_leg==2:
        ax4 = fig.add_axes((.86, .71, .05, .05))
        rect = patches.Rectangle((.5, 0), .01, 1, edgecolor='indianred', facecolor='indianred', alpha=1, zorder=10)
        ax4.add_patch(rect)
        ax4.text(0.7, 0.4, 'current time', color='k', size=15)

        ax5 = fig.add_axes((.86, .57, .05, .05))
        rect = patches.Rectangle((.45, 0), .01, 1, edgecolor='royalblue', facecolor='royalblue', alpha=1, zorder=10)
        ax5.add_patch(rect)
        rect = patches.Rectangle((.55, 0), .01, 1, edgecolor='royalblue', facecolor='royalblue', alpha=1, zorder=10)
        ax5.add_patch(rect)
        ax5.zorder = 100
        ax5.text(0.7, -.1, 'current CI', color='k', size=15, clip_on=False, zorder=100)

        ax6 = fig.add_axes((.86, .51, .05, .05))
        ax6.scatter([.5, .5], [.4, .5], 50, facecolor='royalblue', color='royalblue')

        ax7 = fig.add_axes((.8721, .39, .05, .05))
        rect = patches.Rectangle((0, .5), .5, .01, edgecolor='gold', facecolor='gold', alpha=1, zorder=1000)
        ax7.add_patch(rect)
        ax7.text(0.7, 0.4, 'p < α', color='k', size=15)

        ax4.axis('off')
        ax5.axis('off')
        ax6.axis('off')
        ax7.axis('off')

        leg_axs.append(ax4)
        leg_axs.append(ax5)
        leg_axs.append(ax6)
        leg_axs.append(ax7)



    return leg_axs

def getOpts(boot_cond1, boot_cond2, cond1, tv):

    diff_boot = boot_cond1 - boot_cond2
    scale_factors=[1.5, 1.5, 3.5, 3.5, 3.5, 3.5, 3.5]
    data_src=[cond1, cond1, boot_cond1, boot_cond1, boot_cond2, diff_boot, diff_boot]
    data_show = [25, 25, 100, 100, 100, 100, 100, 100]
    ks=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7', 'plot8', 'plot9']
    title_list=['original trials', 'resampled trials', 'trimmed ERPs', 'trimmed ERPs for condition 1', 'trimmed ERPs for condition 2', 'difference waves', 'difference waves']
    xlabel = ['time (ms)', 'time (ms)', 'time (ms)', 'time (ms)', 'time (ms)', 'time (ms)', 'time (ms)']
    yticklab=[['N', '1'], ['N', '1'], ['Β', '1'], ['Β', '1'], ['B', '1'], ['B', '1'], ['B', '1']]
    #ytick= [[], [], [], [], [], [], []]
    opts={k: {} for k in ks}

    # offsets
    for sf, ds, k , tit, xl, yt, dshow in zip(scale_factors, data_src, ks[0:7], title_list, xlabel, yticklab, data_show[0:-1]):

        offsets = []
        dmin = ds[range(dshow), :].min()
        dmax = ds[range(dshow), :].max()
        dr = (dmax - dmin) / sf

        for ii in range(dshow):
            offsets.append(ii * dr)

        if k=='plot2':
            miy = offsets[0] + ds[range(dshow), :].min()
            may = offsets[-1] + ds[range(dshow), :].max()

        else:
            # ylims
            miy = offsets[0] + min(ds[0, :])
            may = offsets[-1] + max(ds[dshow-1, :])
            #print(miy, may)

        opts[k]['offsets']=offsets
        opts[k]['ylim']=[miy, may]
        opts[k]['xlim'] = [tv[0], tv[-1]]
        opts[k]['scale']=sf
        opts[k]['title'] = tit
        opts[k]['xlab'] = xl
        opts[k]['xtick_loc']=[-200, 0, 200, 400, 600, 800]
        opts[k]['xtick_lab'] = [-200, 0, 200, 400, 600, 800]
        opts[k]['xlab_fs']=15
        opts[k]['title_fs']=15
        opts[k]['ylab'] = ''
        opts[k]['ylab_fs']=15
        opts[k]['nshow']=dshow

        #calc the tick positions and set their labels
        yloc_up=offsets[-1] + ds[-1, 0]
        yloc_low=offsets[0] + ds[0, 0]
        opts[k]['ytick_loc']=[yloc_low, yloc_up]
        opts[k]['ytick_lab']=yt
        opts[k]['ytick_fs']=15
        opts[k]['xtick_fs'] = 12

    # additional plots
    ks = ['plot8', 'plot9' ]
    title_list=['difference scores', 'confidence intervals']
    xlabel = ['values at ', 'time (ms)']
    #ytick= [[], []]

    for k, tit, xl in zip(ks, title_list, xlabel):

        if k=='plot8':
            #opts[k]['xlim']=[-1, boot_cond1.shape[0]]
            opts[k]['xlim'] = [-1, data_show[-1]]
            #opts[k]['xtick_loc']=[0, boot_cond1.shape[0]-1]
            opts[k]['xtick_loc'] = [0, data_show[-1] - 1]
            opts[k]['xtick_lab']=['1', 'B']
            opts[k]['xtick_fs'] = 15
            opts[k]['ytick_fs'] = 15
            opts[k]['ylab'] = ''
            opts[k]['ytick_loc']=[]
            opts[k]['ytick_lab'] = []

        else:
            opts[k]['xlim'] = [tv[0], tv[-1]]
            opts[k]['xtick_loc'] = [-200, 0, 200, 400, 600, 800]
            opts[k]['xtick_lab'] = [-200, 0, 200, 400, 600, 800]
            opts[k]['xtick_fs'] = 12
            opts[k]['ytick_fs'] = 15
            opts[k]['xticks']=tv
            opts[k]['ylab'] = 'microvolts'
            opts[k]['ytick_loc'] = [-10, 0, 10]
            opts[k]['ytick_lab'] = [-10, 0, 10]


        opts[k]['ylim']=[-10, 10]
        opts[k]['title'] = tit
        opts[k]['xlab'] = xl
        #opts[k]['ytick_lab'] = yt
        #opts[k]['ytick_loc'] = yt
        opts[k]['xlab_fs']=16
        opts[k]['title_fs']=16
        opts[k]['ylab_fs'] = 15




    return opts


