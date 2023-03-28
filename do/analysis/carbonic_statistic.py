import pandas as pd
import json

with open('timelong.json', 'r') as fp:
    dict_timelong = json.load(fp)
    print(dict_timelong)

run(
    time_tot = dict_timelong['timelong(ps)'],
)
