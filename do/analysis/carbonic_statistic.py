import pandas as pd

def run(
    str_data: str,
    time_tot: float, # ps
    str_save: str = None,
):

    print(str_data)
    df_data = pd.read_csv(str_data)
    print(df_data)

    df_save = pd.DataFrame()

    ser_sum = df_data.sum()
    ser_sum['sum'] = time_tot
    ser_sum.name = 'timesum(ps)'
    df_save = pd.concat([df_save, ser_sum], axis=1)

    ser_prop = ser_sum/time_tot
    ser_prop.name = 'prop'
    df_save = pd.concat([df_save, ser_prop], axis=1)

    ser_mean = df_data.mean()
    ser_mean.name = 'timemean(ps)'
    df_save = pd.concat([df_save, ser_mean], axis=1)

    ser_count = df_data.count()
    ser_count.name = 'count'
    df_save = pd.concat([df_save, ser_count], axis=1)

    print(str_save)
    print(df_save)
    if str_save:
        df_save.to_csv(str_save)

run(
    time_tot = 72000,
    str_data = 'carbonic_life_time.csv',
    str_save = 'carbonic_statistic.csv',
)
