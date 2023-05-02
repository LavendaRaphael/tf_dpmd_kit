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
            f'rdf.{str_pair}.0002497_0004978.csv',
            f'rdf.{str_pair}.0004978_0007458.csv',
            f'rdf.{str_pair}.0007458_0009939.csv',
            f'rdf.{str_pair}.0009939_0012419.csv',
        ],
        str_save = f'rdf.{str_pair}.ave.csv'
    )
