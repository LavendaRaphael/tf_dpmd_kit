import pandas as pd

def run(
    list_data: list,
    list_header: list = None,
):

    df_data = pd.read_csv(list_data[0])
    if list_header is None:
        list_header = df_data.columns[1:]

    df_life = pd.DataFrame()
    for str_file in list_data:
        print(str_file)
        df_data = pd.read_csv(str_file)
        df_tmp = pd.DataFrame()
        for header in list_header:
            ser_data = df_data[header]
            list_life = []
            life = 0
            for val in ser_data:
                if val == 1.0:
                    life += 1
                elif life != 0:
                    list_life.append(life)
                    life =0
            if life != 0:
                list_life.append(life)

            ser_life = pd.Series(list_life, name=header, dtype='int64')
            df_tmp = pd.concat([df_tmp, ser_life], axis=1)
        print(df_tmp)
        df_life = pd.concat([df_life, df_tmp], ignore_index=True)
    print(df_life)
    df_life.to_csv('carbonic_lifetime.csv')
    print(df_life.sum())

run(
    list_data = [
        'carbonic_state.product.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic_state.product.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic_state.product.csv',
    ],
)
