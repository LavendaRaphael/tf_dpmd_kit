from tf_dpmd_kit import analysis

list_pair = [
    'h_0.o_w',
    'h_1.o_w',
]

for str_pair in list_pair:
    analysis.hbonds_ave(
        list_file = [
            f'hbonds.{str_pair}.status.csv',
        ],
        str_save = f'hbonds.{str_pair}.ave.csv'
    )
