import pandas as pd

def run(
    list_data: list,
    str_save: str = None,
):

    df_final = None
    for str_file in list_data:
        print(str_file)
        df_data = pd.read_csv(str_file)
        df_sum = df_data.sum().to_frame().drop('frame', axis=0)
        df_sum.columns = ['nframe']
        float_sum = len(df_data)
        df_sum.loc['sum'] = [float_sum]
        df_sum['proportion'] = df_sum['nframe']/float_sum
        print(df_sum)
        if df_final is None:
            df_final = df_sum['nframe'].to_frame()
        else:
            df_final['nframe'] += df_sum['nframe']

    df_final['proportion'] = df_final['nframe']/df_final['nframe']['sum']
    print(str_save)
    print(df_final)
    if str_save:
        df_final.to_csv(str_save)

run(
    list_data = [
        'carbonic_state.product.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic_state.product.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic_state.product.csv',
    ],
    str_save = 'carbonic_status.csv'
)
