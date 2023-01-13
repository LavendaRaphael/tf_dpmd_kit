from tf_dpmd_kit import analysis

list_range = [
    (100000,200000),
    (200000,300000),
    (300000,400000),
    (400000,500000),
    (500000,600000),
    (600000,700000),
    (700000,800000),
]

list_pair = [
    'h_0.o_w',
    'h_1.o_w'
]

for str_pair in list_pair:
    np_save = analysis.hbonds_status(
        str_file = f'hbonds.{str_pair}.csv',
        list_range = list_range,
        str_save = f'hbonds.{str_pair}.status.csv'
    )
    print(np_save)

