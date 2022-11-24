from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

dict_label = {
    'dist_vp_o_1': r'R(V$_P$O$_1$)',
    'dist_vp_o_2': r'R(V$_P$O$_2$)',
    'dist_vp_o_1_2': r'R(V$_P$O$_C$)',
    'dist_vp_c': r'R(V$_P$C)',

    'dist_o_0_h': r'R(O$_0$H$_0$)',
    'dist_o_2_h': r'R(O$_2$H$_1$)',

    'cn_o_h': r'CN(O$_C$H)',
    'cn_o_0_h': r'CN(O$_0$H$_w$)',

    'dhx_o_0_h': 'dhx(O$_0$H$_0$)',
    'dhx_o_2_h': 'dhx(O$_2$H$_1$)',

    'metadbias': 'V(s,t)',
    'metadrbias': 'M(s,t)', 
    'metadrct': 'c(t)', 
}

def gen_dict_header(
    list_header: list,
    dict_label: dict,
) -> dict:
    dict_header = {}
    for str_header in list_header:
        if str_header in dict_label:
            dict_header[str_header] = dict_label[str_header]
        else:
            dict_header[str_header] = str_header
    return dict_header


list_header = ['cn_o_h']

'''
fig, axs = plm.colvar_plt(
    dict_header = gen_dict_header(list_header, dict_label),
    #str_save = 'cn_o_h.time.png',
    tup_xlim = (0,None),
    tup_ylim = (-1,4.0),
    #float_timescale = 1/0.000484
    #float_timescale = 1/0.005
)

ax = axs[0]

plm.insert_img(
    fig = fig,
    ax = ax,
    dict_img = {
        '../img/5000.png': (0,0,0.2,0.3),
        '../img/11000.png': (0.2,0,0.2,0.3),
        '../img/100000.png': (0.4,0,0.2,0.3),
        '../img/130000.png': (0.6,0,0.2,0.3),
        '../img/170000.png': (0.8,0,0.2,0.3),
    },
    dict_arrow = {
        (5000*0.005,1.5): (0.03,0.3),
        (11000*0.005,1.0): (0.3,0.3),
        (100000*0.005,1.5): (0.5,0.3),
        (130000*0.005,1.0): (0.7,0.3),
        (170000*0.005,1.5): (0.9,0.3),
    },
    str_save = 'cn_o_h.time.img.png'
)
#'''

'''
plm.insert_img(
    fig = fig,
    ax = ax,
    dict_img = {
        '004000.png': (  0,  0,0.2,0.3),
        '010000.png': (0.0,0.7,0.2,0.3),
        '025000.png': (0.2,  0,0.2,0.3),
        '040000.png': (0.2,0.7,0.2,0.3),
        '050000.png': (0.4,  0,0.2,0.3),
    },
    dict_arrow = {
        (  4000*0.005,1.5): (0.03,0.3),
        ( 10000*0.005,2.0): (0.1,0.7),
        ( 25000*0.005,1.0): (0.3,0.3),
        ( 40000*0.005,2.0): (0.3,0.7),
        ( 50000*0.005,1.5): (0.5,0.3),
    },
    str_save = 'cn_o_h.time.img.png'
)
#'''

'''
fig, axs = plm.colvar_plt(
    dict_header = gen_dict_header(list_header, dict_label),
    tup_xlim = (21.9,22.3),
    tup_ylim = (1,2),
)
ax = axs[0]
plm.insert_img(
    fig,ax,
    dict_img = {
        '004390.png': (  0,0,0.2,0.3),
        '004422.png': (0.3,0,0.2,0.3),
        '004440.png': (0.6,0,0.2,0.3),
    },
    dict_arrow = {
        (4390*0.005,1.7): (0.1,0.3),
        (4422*0.005,1.45): (0.4,0.3),
        (4440*0.005,1.75): (0.7,0.3),
    },
    str_save = 'cn_o_h.time_21_22.img.png'
)
#'''

#'''
fig, axs = plm.colvar_plt(
    dict_header = gen_dict_header(list_header, dict_label),
    tup_xlim = (219.8,220.2),
    tup_ylim = (1,2),
)
ax = axs[0]
plm.insert_img(
    fig,ax,
    dict_img = {
        '043980.png': (  0,0,0.2,0.3),
        '043993.png': (0.3,0,0.2,0.3),
        '044005.png': (0.6,0,0.2,0.3),
    },
    dict_arrow = {
        (43980*0.005,1.7): (0.1,0.3),
        (43993*0.005,1.75): (0.4,0.3),
        (44005*0.005,1.75): (0.7,0.3),
    },
    str_save = 'cn_o_h.time_219_220.img.png'
)
#'''

plt.show()
