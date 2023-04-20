from tf_dpmd_kit import analysis

list_pair = [
    'cc_h_oh.o_w',
    'cc_o_oh.h_w',
    'cc_o_c.h_w',
    'o_w.o_w',
]

for str_pair in list_pair:
    analysis.error(
        list_file = [
            f'rdf.{str_pair}.0800000_1200000.csv',
            f'rdf.{str_pair}.1200000_1600000.csv',
            f'rdf.{str_pair}.1600000_2000000.csv',
        ],
        str_save = f'rdf.{str_pair}.ave.csv'
    )
