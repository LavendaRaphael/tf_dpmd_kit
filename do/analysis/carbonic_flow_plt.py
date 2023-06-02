import matplotlib.pyplot as plt
from tf_dpmd_kit import analysis
from tf_dpmd_kit import plot
import matplotlib as mpl
import numpy as np
import pandas as pd
import json
import matplotlib.colors as mcolors

def run(ax):

    fq = {}
    df = pd.read_csv('carbonic_statistic.csv', index_col='state')['rate(M/s)']
    fq['cc'] = df['CC']
    fq['ct'] = df['CT']
    fq['tt'] = df['TT']
    fq['xx'] = df['HCO3']

    df = pd.read_csv('carbonic_flow.csv', index_col=['from','to'])['rate(M/s)']
    dict_convert = {
        'cc_cc': ('CC', 'CT'),
        'cc_cc': ('CC', 'CT'),
        'cc_ct': ('CC', 'CT'),
        'cc_xx': ('CC', 'HCO3'),
        'cc_tt': ('CC', 'TT'),
        'ct_cc': ('CT', 'CC'),
        'ct_ct': ('CT', 'CT'),
        'ct_xx': ('CT', 'HCO3'),
        'ct_tt': ('CT', 'TT'),
        'xx_cc': ('HCO3', 'CC'),
        'xx_ct': ('HCO3', 'CT'),
        'xx_tt': ('HCO3', 'TT'),
        'tt_ct': ('TT', 'CT'),
        'tt_xx': ('TT', 'HCO3'),
    }
    for key, values in dict_convert.items():
        if values in df:
            fq[key] = df[values]
        else:
            fq[key] = 0

    rate = {}
    for key, value in fq.items():
        value /= 1e8
        rate[key] = f'{value:.3f}'

    plot.add_text(
        ax,
        dict_text = {
            #(0.9, 0.95): r'Frequency (ns$^{-1}$)',
            (0.9, 0.95): r'Rate ($\times 10^{8}$ M/s)',
        },
        va = 'top',
        ha = 'right',
    )

    ax.axis('off')

    # horizon fig
    w, h = ax.bbox.width, ax.bbox.height
    dx = 0.2
    dy = 1200/1600*w/h*dx
    # vertical fig
    dx1 = dy*h/w
    dy1 = dx*w/h

    # right, top
    r = np.array([ dx/2,    0])
    t = np.array([    0,   dy/2])
    r1 = np.array([ dx1/2,    0])
    t1 = np.array([    0,   dy1/2])
    # shift
    sx = np.array([dx/5,   0])
    sy = np.array([   0,dy/5])

    # pos
    p_tt = np.array([0.3, 0.13])
    p_ct = np.array([0.3, 0.5])
    p_cc = np.array([0.3, 0.87])
    p_xx = np.array([0.8, 0.5])

    # middle
    p_tt_ct = (p_tt+p_ct)/2
    p_ct_cc = (p_ct+p_cc)/2
    p_xx_cc = (p_xx-r1+t1+p_cc+r)/2
    p_xx_ct = (p_xx-r1+p_ct+r)/2
    p_xx_tt = (p_xx-r1-t1+p_tt+r)/2
    # color
    alpha = 0.7
    c = mcolors.to_rgb('tab:blue')
    c_cc = (c[0], c[1], c[2], alpha)
    c = mcolors.to_rgb('tab:orange')
    c_ct = (c[0], c[1], c[2], alpha)
    c = mcolors.to_rgb('tab:green')
    c_tt = (c[0], c[1], c[2], alpha)
    c = mcolors.to_rgb('tab:purple')
    c_xx = (c[0], c[1], c[2], alpha)

    # img

    dir_cp  = '/home/faye/research_d/202203_MDCarbonicAcid/server/01.init/H2CO3_CC_H2O_126/plm/'
    dir_tt  = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_nvt_velocity/330K/TT/plm/'
    dir_cc = '/home/faye/research_d/202203_MDCarbonicAcid/server/04.md_npt/330K/CC/snap/'
    [axin3] = plot.inset_img(
        ax,
        dict_img = {
            dir_cc + '0.360003.png': (p_xx[0]-0.5*dx1, p_xx[1]-dy1/2, dx1, dy1),
        },
        dict_spinecolor = {
            dir_cc + '0.360003.png': c_xx,
        },
        bool_rot90 = True,
    )
    [axin0, axin1, axin2] = plot.inset_img(
        ax,
        dict_img = {
            dir_cc +'3.000000.png': (p_cc[0]-0.5*dx, p_cc[1]-0.5*dy, dx, dy),
            dir_cc +'3.067520.png': (p_ct[0]-0.5*dx, p_ct[1]-0.5*dy, dx, dy),
            dir_cc +'3.170878.png': (p_tt[0]-0.5*dx, p_tt[1]-0.5*dy, dx, dy),
        },
        dict_spinecolor = {
            dir_cc +'3.000000.png': c_cc,
            dir_cc +'3.067520.png': c_ct,
            dir_cc +'3.170878.png': c_tt,
        }
    )

    headstyle = 'head_length=4, head_width=4'
    # CC
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_cc-t-sx, p_ct+t-sx],
            [p_cc+r+sy, p_xx+t1],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = c_cc,
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_cc-r, p_tt-r],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.5',
        color = c_cc,
    )
    plot.add_text(
        axin0,
        dict_text = {
            (0.03, 0.95): rate['cc'],
        },
        transform = axin0.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc=c_cc, lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_ct_cc-sx+0.3*sy): rate['cc_ct'],
            tuple(p_xx_cc+sy): rate['cc_xx'],
            #tuple(p_ct_cc-1.7*r): f['cc_tt'],
            #tuple(p_cc-1.6*r+sy): f['cc_cen']
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec=c_cc, fc='white')
    )
    # CT
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct+t+sx, p_cc-t+sx],
            [p_ct-t-sx, p_tt+t-sx],
            [p_ct+r+sy, p_xx-r1+sy],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = c_ct,
    )
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_ct-r+[0,0.1], p_ct-r-[0,0.1]],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        connectionstyle = 'arc3, rad=0.4',
        color = c_ct,
    )
    plot.add_text(
        axin1,
        dict_text = {
            (0.03, 0.95): rate['ct'],
        },
        transform = axin1.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc=c_ct, lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_tt_ct-sx+0.3*sy): rate['ct_tt'],
            tuple(p_ct_cc+sx-0.3*sy): rate['ct_cc'],
            tuple(p_xx_ct+sy): rate['ct_xx'],
            tuple(p_ct-r-sx): rate['ct_ct'],
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec=c_ct, fc='white')
    )

    # TT
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_tt+t+sx, p_ct-t+sx],
            [p_tt+r+sy, p_xx-r1-sy*3],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = c_tt,
    )
    plot.add_text(
        axin2,
        dict_text = {
            (0.03, 0.95): rate['tt'],
        },
        transform = axin2.transAxes,
        va = 'top',
        ha = 'left',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc=c_tt, lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_tt_ct+sx-0.3*sy): rate['tt_ct'],
            tuple(p_xx_tt+sy): rate['tt_xx'],
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec=c_tt, fc='white')
    )
    # HCO3
    plot.add_arrow(
        ax,
        list_arrow = [
            [p_xx-r1+sy*3, p_cc+r-sy],
            [p_xx-r1-sy,   p_ct+r-sy],
            [p_xx-t1, p_tt+r-sy],
        ],
        arrowstyle = f'simple, {headstyle}, tail_width=0.2',
        color = c_xx,
    )
    plot.add_text(
        axin3,
        dict_text = {
            (0.95, 0.97): rate['xx'],
        },
        transform = axin3.transAxes,
        va = 'top',
        ha = 'right',
        color = 'white',
        fontweight = 'bold',
        bbox = dict(boxstyle='round', fc=c_xx, lw=0)
    )
    plot.add_text(
        ax,
        dict_text = {
            tuple(p_xx_ct-sy): rate['xx_ct'],
            tuple(p_xx_cc-sy): rate['xx_cc'],
            tuple(p_xx_tt-sy): rate['xx_tt'],
        },
        va = 'center',
        ha = 'center',
        bbox = dict(boxstyle='round', ec=c_xx, fc='white')
    )

    #
    #plot.add_arrow(
    #    ax,
    #    list_arrow = [
    #        [p_cc-r+sy, p_cc-2.2*r+sy],
    #    ],
    #    arrowstyle = f'-',
    #    color = 'tab:blue',
    #    shrinkA=5, shrinkB=5,
    #    linestyle = ':'
    #),
    #plot.add_arrow(
    #    ax,
    #    list_arrow = [
    #        [p_cc-r+sy, p_cc-2.2*r+sy],
    #    ],
    #    arrowstyle = f'<|-|>, head_length=4, head_width=2',
    #    color = 'tab:blue',
    #    linewidth = 0,
    #)
    #plot.add_text(
    #    ax,
    #    dict_text = {
    #        tuple(p_cc-1.6*r+2*sy): 'Censored',
    #    },
    #    va = 'center',
    #    ha = 'center',
    #)

def main():

    plot.set_rcparam()
    cm = 1/2.54
    mpl.rcParams['figure.dpi'] = 300
    fig, ax = plt.subplots(figsize=(8.6*cm, 6*cm))

    run(ax)

    plot.save(
        fig,
        file_save = 'carbonic_flow',
        list_type = ['pdf', 'svg']
    )

    plt.show()

main()

