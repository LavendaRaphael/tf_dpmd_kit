from tf_dpmd_kit import plm
import matplotlib.pyplot as plt

dict_label = {
    'dist_vp_o_1': r'R(V$_p$O$_1$)',
    'dist_vp_o_2': r'R(V$_p$O$_2$)',
    'dist_vp_o_1_2': r'R(V$_p$O$_C$)',
    'dist_vp_c': r'R(V$_p$C)',
    'cost_o_h': 'cos dh(OH)',
    'dist_o_0_h': r'R(O$_0$H)',
    'metadbias': 'V(s,t)',
    'metadrbias': 'M(s,t)', 
    'metadrct': 'c(t)', 
}

plm.colvar_plt(
    list_header = [
        #['dist_vp_o_1'],
        #['dist_vp_o_2'],
        #['dist_vp_o_1_2'],

        #['dist_vp_c'],
        #['dist_vp_c', (0.5, 1.5)],
        #['lw_vp_cbias'],

        #['cn_o_1_2_h'],
        #['cn_o_h'],

        #['cn_c_h'],

        #['cn_c_o'],

        ['dhx_o_0_h'],
        ['dhx_o_2_h'],

        ['metadbias'],
        #['metadrbias'],
        #['metadrct'],

        ['dist_o_0_h'],
        ['dist_o_2_h'],
        #['dist_o_0_h', (1.1, 1.2)],
        #['uw_dist_o_0_hbias'],

        #['cn_o_0_h'],
        #['uw_cn_o_0_hbias'],
    ],
    float_timescale = 1/1000,
    #float_timescale = 1/0.0005/50,
    #float_timescale = 1/0.000048378/10
    #tup_xlim = (360800, 363300),
    #str_save = 'cv.png',
    dict_label = dict_label
)

plt.show()
