from matplotlib import pyplot as plt
from tf_dpmd_kit import analysis,plm

analysis.rdf_plt_compare(
    dict_data = {
        '280K': '280K/proton/msd.csv',
        '290K': '290K/proton/msd.csv',
        '300K': '300K/proton/msd.csv',
        '310K': '310K/proton/msd.csv',
        '320K': '320K/proton/msd.csv',
        '330K': '330K/proton/msd.csv',
    },
    #tup_xlim = ,
    #tup_ylim = dict_ylim[str_key],
    str_save = 'msd_temperature.pdf',
    #str_title = ,
    str_xlabel = 'Tims (ps)',
    str_ylabel = r'MSD (Å$^2$)',
    dict_temperature = {
        '280K': 280,
        '290K': 290,
        '300K': 300,
        '310K': 310,
        '320K': 320,
        '330K': 330,
    },
    tup_colormap = (280,330)
)

plt.show() 
