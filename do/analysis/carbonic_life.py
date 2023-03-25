import pandas as pd

def run(
    list_data: list,
    timestep: float,
    list_header: list = None,
    intermit_frame: int = 0,
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
            intermit = 0
            for val in ser_data:
                if val == 1.0:
                    if intermit > intermit_frame:
                        if life > 0:
                            list_life.append(life)
                            life = 0
                    intermit = 0
                    life += 1
                else:
                    intermit += 1
                    if life != 0 and life < intermit_frame:
                        life = 0
            if life != 0:
                list_life.append(life)

            ser_life = pd.Series(list_life, name=header, dtype='int64')
            df_tmp = pd.concat([df_tmp, ser_life], axis=1)
        print(df_tmp)
        df_life = pd.concat([df_life, df_tmp], ignore_index=True)

    file_lifeframe = 'carbonic_life_frame.csv'
    print(file_lifeframe)
    df_life.to_csv(file_lifeframe, index=False)
    print(df_life)

    file_lifetime = 'carbonic_life_time.csv'
    df_life = df_life*timestep
    df_life.to_csv(file_lifetime, index=False)
    print(file_lifetime)
    print(df_life)
'''
run(
    list_data = [
        'carbonic_state.csv',
    ],
    timestep = 0.0004837769,
    intermit_frame = 1000,
)

#'''
#'''
run(
    list_data = [
        'carbonic_state.product.csv',
        '../../../H2CO3_CT_H2O_126/330K/carbonic/carbonic_state.product.csv',
        '../../../H2CO3_CC_H2O_126/330K/carbonic/carbonic_state.product.csv',
    ],
    timestep = 0.005,
    intermit_frame = 100,
)
#'''
