import pandas as pd

df_data = pd.read_csv('carbonic_state.csv')

df_new = df_data[df_data['CC']==1.0]
print(df_new)

df_new.to_csv('carbonic_filter.csv', index=False)
