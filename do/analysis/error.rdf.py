from tf_dpmd_kit import analysis

list_pair = [
    'h_0_1.o_w',
    'o_1.h_w',
    'o_w.o_w',
]
str_dir_tmp = '../../../H2CO3_TT_H2O_126/330K/rdf/'

for str_pair in list_pair:
    analysis.error(
        list_file = [
            f'rdf.{str_pair}.0100000_0200000.csv',
            f'rdf.{str_pair}.0200000_0300000.csv',
            #f'rdf.{str_pair}.0300000_0400000.csv',
            str_dir_tmp+f'rdf.{str_pair}.0100000_0200000.csv',
            str_dir_tmp+f'rdf.{str_pair}.0200000_0300000.csv',
            str_dir_tmp+f'rdf.{str_pair}.0300000_0400000.csv',
        ],
        str_save = f'rdf.{str_pair}.ave.csv'
    )
