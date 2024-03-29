import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from tf_dpmd_kit import plm
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import MDAnalysis.analysis.msd
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_dihedrals
import pandas as pd
import time
import json
from lifelines import KaplanMeierFitter

class CarbonicRDF(AnalysisBase):
    """
    Modify from Mdanalysis 2.4.0 InterRDF
    """
    def __init__(self,
                 water_o,
                 all_h,
                 carbonic,
                 nbins=200,
                 range=(0.0, 15.0),
                 **kwargs):

        super(CarbonicRDF, self).__init__(water_o.universe.trajectory, **kwargs)
        self.water_o = water_o
        self.all_h = all_h
        self.all_atoms = water_o.universe.atoms
        self.carbonic = carbonic

        self.rdf_settings = {'bins': nbins,
                             'range': range}

    def _prepare(self):
        # Empty histogram to store the RDF
        count, edges = np.histogram([-1], **self.rdf_settings)
        count = count.astype(np.float64)
        count *= 0.0
        self.results.edges = edges
        bins = 0.5 * (edges[:-1] + edges[1:])

        self.results.count = pd.DataFrame(
            data={
                'cc.o_nyl.h_w' : count,
                'cc.o_oh.h_w'  : count,
                'cc.h_oh.o_w'  : count,
                'ct.o_nyl.h_w' : count,
                'ct.o_oh_c.h_w': count,
                'ct.o_oh_t.h_w': count,
                'ct.h_oh_c.o_w': count,
                'ct.h_oh_t.o_w': count,
                'tt.o_nyl.h_w' : count,
                'tt.o_oh.h_w'  : count,
                'tt.h_oh.o_w'  : count,
                'o_w.o_w'  : count,
            },
            index = bins
        )
        self.results.count.index.name = 'r(ang)'
        self.results.nframes = pd.Series(data={'cc': 0, 'ct': 0, 'tt': 0}, name='nframe')

        # Cumulative volume for rdf normalization
        self.volume_cum = 0
        # Set the max range to filter the search radius
        self._maxrange = self.rdf_settings['range'][1]
        self._minrange = self.rdf_settings['range'][0]

    def _single_frame(self):
        self.volume_cum += self._ts.volume

        ser_carbonic= self.carbonic.loc[self._ts.frame]
        ncarbonyl = int(ser_carbonic['ncarbonyl'])
        if ncarbonyl != 1:
            # HCO3, CO3, H3CO3
            return

        dh0  = ser_carbonic['dh0(rad)']
        dh1  = ser_carbonic['dh1(rad)']
        pio2 = np.pi/2
        bool_dh0 = ((dh0 > -pio2) & (dh0 < pio2))
        bool_dh1 = ((dh1 > -pio2) & (dh1 < pio2))

        dho0 = int(ser_carbonic['dho0'])
        dh0o = int(ser_carbonic['dh0o'])
        dh1o = int(ser_carbonic['dh1o'])
        dh0h = int(ser_carbonic['dh0h'])
        dh1h = int(ser_carbonic['dh1h'])
        o_nyl = self.all_atoms[[dho0]]
        h_w = self.all_h - self.all_atoms[[dh0h, dh1h]]
        if bool_dh0:
            if bool_dh1:
                # CC
                o_oh = self.all_atoms[[dh0o, dh1o]]
                h_oh = self.all_atoms[[dh0h, dh1h]]
                self.results.count['cc.o_nyl.h_w'] += self._cal_count(o_nyl, h_w)
                self.results.count['cc.o_oh.h_w' ] += self._cal_count(o_oh , h_w)
                self.results.count['cc.h_oh.o_w' ] += self._cal_count(h_oh , self.water_o)
                self.results.nframes['cc'] += 1
            else:
                # CT
                o_oh_c = self.all_atoms[[dh0o]]
                o_oh_t = self.all_atoms[[dh1o]]
                h_oh_c = self.all_atoms[[dh0h]]
                h_oh_t = self.all_atoms[[dh1h]]
                self.results.count['ct.o_nyl.h_w' ] += self._cal_count(o_nyl , h_w)
                self.results.count['ct.o_oh_c.h_w'] += self._cal_count(o_oh_c, h_w)
                self.results.count['ct.o_oh_t.h_w'] += self._cal_count(o_oh_t, h_w)
                self.results.count['ct.h_oh_c.o_w'] += self._cal_count(h_oh_c, self.water_o)
                self.results.count['ct.h_oh_t.o_w'] += self._cal_count(h_oh_t, self.water_o)
                self.results.nframes['ct'] += 1
        else:
            if bool_dh1:
                # TC
                o_oh_t = self.all_atoms[[dh0o]]
                o_oh_c = self.all_atoms[[dh1o]]
                h_oh_t = self.all_atoms[[dh0h]]
                h_oh_c = self.all_atoms[[dh1h]]
                self.results.count['ct.o_nyl.h_w' ] += self._cal_count(o_nyl , h_w)
                self.results.count['ct.o_oh_c.h_w'] += self._cal_count(o_oh_c, h_w)
                self.results.count['ct.o_oh_t.h_w'] += self._cal_count(o_oh_t, h_w)
                self.results.count['ct.h_oh_c.o_w'] += self._cal_count(h_oh_c, self.water_o)
                self.results.count['ct.h_oh_t.o_w'] += self._cal_count(h_oh_t, self.water_o)
                self.results.nframes['ct'] += 1
            else:
                # TT
                o_oh = self.all_atoms[[dh0o, dh1o]]
                h_oh = self.all_atoms[[dh0h, dh1h]]
                self.results.count['tt.o_nyl.h_w'] += self._cal_count(o_nyl, h_w)
                self.results.count['tt.o_oh.h_w' ] += self._cal_count(o_oh , h_w)
                self.results.count['tt.h_oh.o_w' ] += self._cal_count(h_oh , self.water_o)
                self.results.nframes['tt'] += 1
        self.results.count['o_w.o_w'] += self._cal_count(self.water_o, self.water_o)

    def _cal_count(self, g1, g2):

        pairs, dist = capped_distance(
                    g1,
                    g2,
                    self._maxrange,
                    self._minrange,
                    box=self._ts.dimensions)
        count, _ = np.histogram(dist, **self.rdf_settings)
        return count

    def _conclude(self):
        # Volume in each radial shell
        vols = np.power(self.results.edges, 3)
        norm = 4/3 * np.pi * np.diff(vols)

        # Average number density
        self.results.box_vol = self.volume_cum / self.n_frames
        norm /= self.results.box_vol

        # Number of each selection
        n_h_w = len(self.all_h)-2
        n_o_w = len(self.water_o)

        nframe_cc = self.results.nframes['cc']
        nframe_ct = self.results.nframes['ct']
        nframe_tt = self.results.nframes['tt']

        self.results.rdf = self.results.count
        self.results.rdf['cc.o_nyl.h_w' ] = self.results.count['cc.o_nyl.h_w' ] / (norm *nframe_cc *n_h_w)
        self.results.rdf['cc.o_oh.h_w'  ] = self.results.count['cc.o_oh.h_w'  ] / (norm *nframe_cc *n_h_w *2)
        self.results.rdf['cc.h_oh.o_w'  ] = self.results.count['cc.h_oh.o_w'  ] / (norm *nframe_cc *n_o_w *2)
        self.results.rdf['ct.o_nyl.h_w' ] = self.results.count['ct.o_nyl.h_w' ] / (norm *nframe_ct *n_h_w)
        self.results.rdf['ct.o_oh_c.h_w'] = self.results.count['ct.o_oh_c.h_w'] / (norm *nframe_ct *n_h_w)
        self.results.rdf['ct.o_oh_t.h_w'] = self.results.count['ct.o_oh_t.h_w'] / (norm *nframe_ct *n_h_w)
        self.results.rdf['ct.h_oh_c.o_w'] = self.results.count['ct.h_oh_c.o_w'] / (norm *nframe_ct *n_o_w)
        self.results.rdf['ct.h_oh_t.o_w'] = self.results.count['ct.h_oh_t.o_w'] / (norm *nframe_ct *n_o_w)
        self.results.rdf['ct.o_oh.h_w'  ] = (self.results.count['ct.o_oh_c.h_w'] + self.results.count['ct.o_oh_t.h_w'])/2
        self.results.rdf['ct.h_oh.o_w'  ] = (self.results.count['ct.h_oh_c.o_w'] + self.results.count['ct.h_oh_t.o_w'])/2
        self.results.rdf['tt.o_nyl.h_w' ] = self.results.count['tt.o_nyl.h_w' ] / (norm *nframe_tt *n_h_w)
        self.results.rdf['tt.o_oh.h_w'  ] = self.results.count['tt.o_oh.h_w'  ] / (norm *nframe_tt *n_h_w *2)
        self.results.rdf['tt.h_oh.o_w'  ] = self.results.count['tt.h_oh.o_w'  ] / (norm *nframe_tt *n_o_w *2)
        self.results.rdf['o_w.o_w'      ] = self.results.count['o_w.o_w'      ] / (norm *(nframe_cc+nframe_ct+nframe_tt) *n_o_w *n_o_w)
    
def read_multidata(
    list_file: list,
):

    df = pd.DataFrame()
    for file in list_file:
        df_tmp = pd.read_csv(file)
        df = pd.concat([df, df_tmp])

    return df

def carbonic_lifetime(
    df,
    file_lifetime: str = 'carbonic_lifetime.csv',
):
    '''
    Lifetime = integral[s(t)*t]/integral[s(t)]
    '''

    dfgb = df.groupby(level='state')
    df_lifetime = pd.DataFrame()
    for state, gp in dfgb:

        timeline = gp.index.get_level_values('timeline(ps)')
        survival = gp['survival']
        lower = gp['lower']
        upper = gp['upper']

        lifetime = np.trapz(y=survival*timeline, x=timeline) / np.trapz(y=survival, x=timeline)
        error_lower = lifetime - np.trapz(y=lower*timeline, x=timeline) / np.trapz(y=lower, x=timeline)
        error_upper = np.trapz(y=upper*timeline, x=timeline) / np.trapz(y=upper, x=timeline) - lifetime
        df_tmp = pd.DataFrame(data={'lifetime(ps)': lifetime, 'lower': error_lower, 'upper': error_upper}, index=[state])
        df_tmp.index.name = 'state'
        df_lifetime = pd.concat([df_lifetime, df_tmp])

    print(file_lifetime)
    print(df_lifetime)
    df_lifetime.to_csv(file_lifetime)

def carbonic_survival_plt(
    ax,
    df,
    list_state: list,
    dict_color: dict = None,
    dict_label: dict = None,
):

    dfgb = df.groupby(level='state')
    for state in list_state:

        color = dict_color[state]
        label = state
        if state in dict_label:
            label = dict_label[state]

        gp = dfgb.get_group(state)
        timeline = gp.index.get_level_values('timeline(ps)')
        survival = gp['survival']
        lower = gp['lower']
        upper = gp['upper']

        ax.plot(timeline, survival, label=label, color=color, lw=1)
        ax.fill_between(timeline, lower, upper, alpha=0.5, color=color, lw=0)

    ax.legend(frameon=False, labelspacing=0.3, handlelength=1)
    ax.set_xscale('log')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Survial Probability')
    ax.set_ylim(0,1)

def carbonic_survival(
    list_file: list,
    func_step: float = 0.5,
):

    '''
    From the life data to calculate survival probability using Kaplan-Meier estimitor.
    '''

    print(list_file)
    df_data = read_multidata(list_file).loc[:, ['state','time(ps)','event']].dropna()

    kmf = KaplanMeierFitter()
    grouped = df_data.groupby('state')
    list_df = []
    list_state = []
    for state, gp in grouped:
        t_max = max(gp['time(ps)'])
        kmf.fit(gp['time(ps)'], gp['event'], timeline=np.linspace(0, t_max, num=int(t_max/func_step)))
        df_tmp = pd.concat([kmf.survival_function_, kmf.confidence_interval_], axis=1)
        df_tmp.rename(
            columns = {
                'KM_estimate': 'survival',
                'KM_estimate_lower_0.95': 'lower',
                'KM_estimate_upper_0.95': 'upper',
            }, 
            inplace=True
        )
        df_tmp.index.name = 'timeline(ps)'
        list_df.append( df_tmp )
        list_state.append(state)
        print(f'state: {state}, median survial time: {kmf.median_survival_time_}')
    df = pd.concat(list_df, keys=list_state, names=['state'])
    print(df)

    return df

def carbonic_statistic_mean(
    list_file: list,
    file_save = 'carbonic_statistic.csv',
):

    '''
    From the statictic data to calculate the mean and sem.
    '''

    list_df = []
    list_key = []
    for file in list_file:
        df_tmp = pd.read_csv(file, index_col=0)
        list_df.append(df_tmp)
        list_key.append(file)
    df_data = pd.concat(list_df, keys=list_key)
    df_data = df_data.fillna(0)
    print(df_data)

    df_save = df_data.groupby(level=1).mean()

    df_sem = df_data.groupby(level=1).sem()
    dict_column = {}
    for column in df_sem.columns:
        dict_column[column] = column+'_sem'
    df_sem.rename(columns=dict_column, inplace=True)

    df_save = pd.concat([df_save, df_sem], axis=1)
    print(file_save)
    print(df_save)
    df_save.to_csv(file_save)

def carbonic_statistic(
    time_tot: float, # ps
    volume: float, # ang3
    file_data: str = 'carbonic_lifedata.csv',
    file_state: str = 'carbonic_state.product.csv',
    file_save: str = 'carbonic_statistic.csv',
):

    '''
    From cabonic_state to calcualte state fraction & from lifedata to calulate frequency.
    '''

    print(file_state)
    df_data = pd.read_csv(file_state, index_col=0)
    ser_count = df_data.count()/len(df_data)
    df_save = ser_count.to_frame().rename(columns={0: 'frac'})

    print(file_data)
    df_data = pd.read_csv(file_data).loc[:, ['state', 'time(ps)']]
    print(df_data)
    gpby = df_data.groupby('state')

    #df_tmp = gpby.sum()/time_tot
    #df_tmp.rename(columns={'time(ps)': 'frac'}, inplace=True)
    #df_save = df_tmp

    df_tmp = gpby.count()/time_tot*1000
    df_tmp.rename(columns = {'time(ps)': 'freq(ns-1)'}, inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    Avogadro = 6.02214076e23
    molar = 1e27/Avogadro/volume
    ns2s = 1e9
    df_tmp = df_save['freq(ns-1)']*molar*ns2s
    df_tmp.rename('rate(M/s)', inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    count_tot = df_data[df_data['state']!='H2CO3']['time(ps)'].count()
    print(count_tot)
    df_tmp = gpby.count()/count_tot
    df_tmp.rename(columns = {'time(ps)': 'freqfrac'}, inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    print(file_save)
    print(df_save)
    if file_save:
        df_save.to_csv(file_save)

def carbonic_statistic_(
    time_tot: float, # ps
    volume: float, # ang3
    file_data: str = 'carbonic_lifedata.csv',
    file_save: str = 'carbonic_statistic.csv',
):

    '''
    From lifedata to calulate state fraction and frequency.
    '''

    print(file_data)
    df_data = pd.read_csv(file_data).loc[:, ['state', 'time(ps)']]
    print(df_data)
    gpby = df_data.groupby('state')

    df_tmp = gpby.sum()/time_tot
    df_tmp.rename(columns={'time(ps)': 'frac'}, inplace=True)
    df_save = df_tmp

    df_tmp = gpby.count()/time_tot*1000
    df_tmp.rename(columns = {'time(ps)': 'freq(ns-1)'}, inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    Avogadro = 6.02214076e23
    molar = 1e27/Avogadro/volume
    ns2s = 1e9
    df_tmp = df_save['freq(ns-1)']*molar*ns2s
    df_tmp.rename('rate(M/s)', inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    count_tot = df_data[df_data['state']!='H2CO3']['time(ps)'].count()
    print(count_tot)
    df_tmp = gpby.count()/count_tot
    df_tmp.rename(columns = {'time(ps)': 'freqfrac'}, inplace=True)
    df_save = pd.concat([df_save, df_tmp], axis=1)

    print(file_save)
    print(df_save)
    if file_save:
        df_save.to_csv(file_save)

def carbonic_lifedata(
    timestep: float,
    list_header: list = None,
    intermit_time: float = 0,
    file_data: str = 'carbonic_state.product.csv',
    file_save: str = 'carbonic_lifedata.csv',
    timelong_save: str = 'timelong.json',
):

    '''
    From carbonic state to life data. 
    '''

    intermit_frame = intermit_time/timestep

    df_life = pd.DataFrame()
    timelong_tot = 0
    dict_timelong = {}

    print(file_data)
    df_data = pd.read_csv(file_data)
    print(df_data)
    list_header = df_data.columns[1:]
    dict_timelong['timelong(ps)'] = len(df_data)*timestep

    list_life = []
    for header in list_header:
        ser_data = df_data[header]
        life = 0
        intermit = 0
        list_life.append([np.nan, header, np.nan, np.nan, np.nan])
        for idx, val in enumerate(ser_data):
            if pd.notnull(val):
                intermit = 0
                life += 1
                if life == 1:
                    start = idx
            else:
                intermit += 1
                if life <= intermit_frame:
                    life = 0
                elif intermit > intermit_frame:
                    list_life.append((start, header, life, 1, idx-intermit))
                    life = 0
        if life > 0:
            list_life.append((start, header, life, 0, idx))

    df_life = pd.DataFrame(list_life, columns=['start', 'state', 'time(ps)', 'event', 'end'])
    df_life['time(ps)'] = df_life['time(ps)']*timestep
    df_life['start'] = df_life['start']*timestep
    df_life['end'] = df_life['end']*timestep
    df_life.loc[df_life['state']=='H2CO3', 'start'] = np.nan
    df_life.sort_values(by=['start'], inplace=True)
    print(file_save)
    print(df_life)
    if file_save:
        df_life.to_csv(file_save, index=False)

    print(timelong_save)
    print(dict_timelong)
    with open(timelong_save, 'w') as fp:
        json.dump(dict_timelong, fp)

def carbonic_rolling_plt(
    ax,
    float_xscale: float = 1,
    file_data: str = 'carbonic_state.csv',
    int_window: int = 1,
    list_header: list = None,
    list_ypos: list = None,
    dict_color: dict = None,
):

    df_data = data_rolling(
        int_window = int_window,
        file_data = file_data,
    )
    print(df_data)
    
    df_new = df_data.where(df_data.isnull(), 1)
    df_new['frame'] = df_data['frame']
    
    if dict_color is None:  
        dict_color = {}

    for idx, header in enumerate(list_header):
        df_data_tmp = df_data[df_data[header].notnull()]
        df_new_tmp = df_new[df_new[header].notnull()]
        if len(df_data_tmp)==0:
            continue
        color = None
        if header in dict_color:
            color = dict_color[header]
        ax.scatter( df_new_tmp['frame']*float_xscale, df_new_tmp[header]*list_ypos[idx], s=2, edgecolors='none', alpha=df_data_tmp[header], rasterized=True, color=color)

    ax.set_yticks(list_ypos)

def data_rolling(
    int_window: int,
    file_data: str,
    file_save: str = None,
):

    print(file_data)
    df_data = pd.read_csv(file_data)
    #print(df_data)
    df_data = df_data.where(df_data.notnull(), 0)
    df_new = df_data.rolling(int_window, min_periods=1, center=True, step=int_window).mean()
    df_new = df_new.where( df_new!=0, None)
    print(file_save)
    #print(df_new)
    if file_save:
        df_new.to_csv(file_save, index=False)

    return df_new

def carbonic_state(
    file_data: str = 'carbonic.csv',
    file_save: str = 'carbonic_state.csv',
):
    start = time.time()

    print(file_data)
    df_data = pd.read_csv(file_data)
    print(df_data)

    df_new = df_data.apply(lambda x: carbonic_evalstate(x['ncarbonyl'], x['dh0(rad)'], x['dh1(rad)']), axis=1, result_type='expand')
    df_new.columns = ['CO3','HCO3','H2CO3','CC','CT','TT','H3CO3']
    df_new.insert(0, 'frame', df_data['frame'])
    
    print(file_save)
    print(df_new)
    df_new.to_csv(file_save, index=False)

    end = time.time()
    print(end-start)

def carbonic_evalstate(
    ncarbonyl,
    dh0,
    dh1,
):
    list_re = [None]*7

    if ncarbonyl == 3:
        # CO3
        list_re[0] = 1
    elif ncarbonyl == 2:
        list_re[1] = 1
    elif ncarbonyl == 1:
        list_re[2] = 1
        list_re[carbonic_conformer(dh0, dh1)] = 1
    elif ncarbonyl == 0:
        list_re[6] = 1
    return list_re

def carbonic_conformer(
    dh0,
    dh1
):
    pio2 = np.pi/2
    bool_dh0 = ((dh0 > -pio2) & (dh0 < pio2))
    bool_dh1 = ((dh1 > -pio2) & (dh1 < pio2))

    if bool_dh0:
        if bool_dh1:
            # CC
            return 3
        else:
            # TC
            return 4
    else:
        if bool_dh1:
            # CT
            return 4
        else:
            # TT
            return 5

class Carbonic(AnalysisBase):

    def __init__(
        self,
        carbonic_c,
        carbonic_o,
        atomg_h,
        water_o,
        cutoff = 1.3,
        ):

        trajectory = carbonic_c.universe.trajectory
        super(Carbonic, self).__init__(trajectory)

        self.u = carbonic_c.universe
        self.carbonic_c = carbonic_c
        self.carbonic_o = carbonic_o
        self.atomg_h = atomg_h
        self.water_o = water_o
        self.cutoff = cutoff

    def _prepare(self):

        self.results = np.zeros((self.n_frames, 11))
        self.carbonic_c2 = self.carbonic_c[[0,0]]

    def _single_frame(self):

        self.results[self._frame_index, 1:] = self._carbonic()
        self.results[self._frame_index, 0] = self._ts.frame

    def _conclude(self):

        columns = ['frame', 'ncarbonyl', 'dh0(rad)', 'dh1(rad)','roh0(ang)','roh1(ang)', 'dho0','dh0o','dh1o','dh0h','dh1h']
        self.df = pd.DataFrame(self.results, columns=columns)
        self.df = self.df.astype(dtype={
            'frame': 'int64',
            'ncarbonyl': 'int64',
            'dho0': 'Int64', 
            'dh0o': 'Int64', 
            'dh1o': 'Int64', 
            'dh0h': 'Int64', 
            'dh1h': 'Int64'})
        # Int64 for Nan
        
    def _carbonic(self):
    
        box = self._ts.dimensions

        carbonyl = self.u.atoms[[]]
        dh_o = self.u.atoms[[]]
        dh_h = self.u.atoms[[]]
        dh_o0 = self.u.atoms[[]]
        roh = []
        for atom_o in self.carbonic_o:
            np_id, np_distances = capped_distance(
                reference = mda.AtomGroup([atom_o]*2),
                configuration = self.atomg_h,
                max_cutoff = self.cutoff,
                box = box,
            )
            if np.size(np_distances) == 0:
                carbonyl += atom_o
                continue
            h_id = np_id[np.argmin(np_distances), 1]
            dh_h += self.atomg_h[h_id]
            dh_o += atom_o
            roh.append(min(np_distances))

        ncarbonyl = len(carbonyl)
        if ncarbonyl == 1:
            dh_o0 = carbonyl
        if ncarbonyl == 2:
            water_h = self.atomg_h - dh_h
            np_id, np_distances = capped_distance(
                reference = self.water_o,
                configuration = water_h,
                max_cutoff = self.cutoff,
                box = box,
            )
            df = pd.DataFrame(np_id)
            df_o = df.groupby(0)
            df_o_count = df_o.count()
            df_o_count = df_o_count[df_o_count[1]==3]
            if len(df_o_count) == 1:
                o_id = df_o_count.index.tolist()[0]
                h_proton = water_h[df_o.get_group(o_id)[1].tolist()]
            elif df_o_count.empty:
                # single H
                h_proton = water_h - water_h[np_id[:,1]]
            else:
                df_h = df.groupby(1)
                df_h_count = df_h.count()
                df_h_count = df_h_count[df_h_count[0]==2]
                o_list = []
                for h_id in df_h_count.index:
                    o_list.extend( df_h.get_group(h_id)[0].tolist())
                h_proton = self.u.atoms[[]]
                for o_id in df_o_count.index:
                    if o_id in o_list:
                        h_proton += self.water_o[o_id]
                    else:
                        h_proton += water_h[df_o.get_group(o_id)[1].tolist()]

            np_id, np_distances = capped_distance(
                reference = carbonyl,
                configuration = h_proton,
                max_cutoff = box[0],
                box = box,
            )
            roh.append(min(np_distances))
            h_id = np_id[np.argmin(np_distances), 1]
            dh_h += h_proton[h_id]
            o_id = np_id[np.argmin(np_distances), 0]
            dh_o += carbonyl[o_id]
            dh_o0 = carbonyl - carbonyl[o_id]
        
        if len(dh_o0) == 1:
            if roh[0] < roh[1]:
                dh_o = dh_o[1] + dh_o[0]
                dh_h = dh_h[1] + dh_h[0]
                roh = [roh[1], roh[0]]
            dh = calc_dihedrals(
                dh_o0 + dh_o0,
                self.carbonic_c2,
                dh_o,
                dh_h,
                box=box
            )
            idx = [dh_o0[0].index, dh_o[0].index, dh_o[1].index, dh_h[0].index, dh_h[1].index]
        else:
            dh = [None, None]
            roh = [None, None]
            idx = [None, None, None, None, None]

        return ncarbonyl, dh[0], dh[1], roh[0], roh[1], idx[0], idx[1], idx[2], idx[3], idx[4]

class Carbonic_(AnalysisBase):

    def __init__(
        self,
        carbonic_c,
        carbonic_o,
        atomg_h,
        water_o,
        cutoff = 1.3,
        ):

        trajectory = carbonic_c.universe.trajectory
        super(Carbonic, self).__init__(trajectory)

        self.u = carbonic_c.universe
        self.carbonic_c = carbonic_c
        self.carbonic_o = carbonic_o
        self.atomg_h = atomg_h
        self.water_o = water_o
        self.cutoff = cutoff

    def _prepare(self):

        self.results = np.zeros((self.n_frames, 6))
        self.carbonic_c2 = self.carbonic_c[[0,0]]

    def _single_frame(self):

        self.results[self._frame_index, 1:] = self._carbonic()
        self.results[self._frame_index, 0] = self._ts.frame

    def _conclude(self):

        columns = ['frame',
                   'ncarbonyl', 'dihedral0(rad)', 'dihedral1(rad)','roh0(ang)','roh1(ang)']
        self.df = pd.DataFrame(self.results, columns=columns)
    
    def _carbonic(self):
    
        box = self._ts.dimensions

        carbonyl = self.u.atoms[[]]
        dh_o = self.u.atoms[[]]
        dh_h = self.u.atoms[[]]
        dh_o0 = self.u.atoms[[]]
        list_dist = []
        for atom_o in self.carbonic_o:
            np_id, np_distances = capped_distance(
                reference = mda.AtomGroup([atom_o]*2),
                configuration = self.atomg_h,
                max_cutoff = self.cutoff,
                box = box,
            )
            if np.size(np_distances) == 0:
                carbonyl += atom_o
                continue
            h_id = np_id[np.argmin(np_distances), 1]
            dh_h += self.atomg_h[h_id]
            dh_o += atom_o
            list_dist.append(min(np_distances))

        ncarbonyl = len(carbonyl)
        if ncarbonyl == 1:
            dh_o0 = carbonyl
        if ncarbonyl == 2:
            water_h = self.atomg_h - dh_h
            np_id, np_distances = capped_distance(
                reference = self.water_o,
                configuration = water_h,
                max_cutoff = self.cutoff,
                box = box,
            )
            df = pd.DataFrame(np_id).groupby(0)
            df_count = df.count()
            df_count = df_count[df_count[1]==3]
            if df_count.empty:
                # single H
                h_proton = water_h - water_h[np_id[:,1]]
            else:
                h_proton = self.u.atoms[[]]
                for o_id in df_count.index:
                    h_proton += water_h[df.get_group(o_id)[1].tolist()]
            np_id, np_distances = capped_distance(
                reference = carbonyl,
                configuration = h_proton,
                max_cutoff = box[0],
                box = box,
            )
            list_dist.append(min(np_distances))
            h_id = np_id[np.argmin(np_distances), 1]
            dh_h += h_proton[h_id]
            o_id = np_id[np.argmin(np_distances), 0]
            dh_o += carbonyl[o_id]
            dh_o0 = carbonyl - carbonyl[o_id]
        
        if len(dh_o0) == 1:
            dihedrals = calc_dihedrals(
                dh_o0 + dh_o0,
                self.carbonic_c2,
                dh_o,
                dh_h,
                box=box
            )
        else:
            dihedrals = [None, None]
            list_dist = [None, None]

        return ncarbonyl, dihedrals[0], dihedrals[1], list_dist[0], list_dist[1]

class HbondLength(AnalysisBase):

    def __init__(
        self,
        hydroxyl_o,
        hydroxyl_h,
        water_o,
        f_cutoff = 4.0):

        trajectory = hydroxyl_o.universe.trajectory
        super(HbondLength, self).__init__(trajectory)

        self.hydroxyl_o = hydroxyl_o
        self.hydroxyl_h = hydroxyl_h
        self.water_o = water_o
        self.f_cutoff = f_cutoff

    def _prepare(self):

        self.results = np.zeros((self.n_frames, 5))

    def _single_frame(self):

        self.results[self._frame_index, 2:] = self._hbondlength()
        self.results[self._frame_index, 0] = self._ts.frame
        self.results[self._frame_index, 1] = self._trajectory.time

    def _conclude(self):

        columns = ['Frame', 'Time(ps)',
                   'OhWaterDistance',
                   'OhWaterAngle',
                   'WaterWaterDistance']
        self.df = pd.DataFrame(self.results, columns=columns)
    
    def _hbondlength(self):
    
        box = self._ts.dimensions

        np_indices, np_distances = capped_distance(
            reference = self.hydroxyl_o,
            configuration = self.water_o,
            max_cutoff = self.f_cutoff,
            min_cutoff = 1.0,
            box = box,
            return_distances = True,
        )
        if np.size(np_indices) == 0:
            raise RuntimeError('No 1st neighbor')
    
        i_indice = np.argmin(np_distances)
        f_1st_distance = np_distances[i_indice]
        
        water_1st = self.water_o[[np_indices[i_indice, 1]]]
        f_1st_angle = calc_angles(
            self.hydroxyl_h.positions,
            self.hydroxyl_o.positions,
            water_1st.positions,
            box=box
        )[0]

        np_indices, np_distances = capped_distance(
            reference = water_1st,
            configuration = self.water_o,
            max_cutoff = self.f_cutoff,
            min_cutoff = 1.0,
            box = box,
            return_distances = True,
        )
        if np.size(np_indices) == 0:
            raise RuntimeError('No 2st neighbor')
    
        i_indice = np.argmin(np_distances)
        f_2st_distance = np_distances[i_indice]
    
        return f_1st_distance, f_1st_angle, f_2st_distance

def error(
    list_file: str,
    str_save: str,
) -> None:

    print('data: ', list_file[0])
    np_0 = np.loadtxt(list_file[0])

    np_ave = np_0[:,1]
    np_ave2 = np_ave**2

    for str_file in list_file[1:]:
        print('data: ', str_file)
        np_tmp = np.loadtxt(str_file)
        np_ave += np_tmp[:,1]
        np_ave2 += np_tmp[:,1]**2

    int_n = len(list_file)
    np_ave /= int_n
    np_ave2 /= int_n

    # for the value inf in grid file
    np.seterr( invalid='ignore' )

    np_ssd = np.sqrt(int_n/(int_n-1)*(np_ave2 - np_ave**2))
    np_error = np_ssd/np.sqrt(int_n)

    np_final = np.zeros( shape=(np_0.shape[0], 3) )
    np_final[:,0] = np_0[:,0]
    np_final[:,1] = np_ave
    np_final[:,2] = np_error

    with open(list_file[0], 'r') as fp:
        str_header = fp.readline()[1:-1]

    print('save: ', str_save)
    np.savetxt(
        X = np_final,
        fname = str_save,
        header = str_header+' error'
    )

def msd(
    mda_u,
    str_select: str,
    float_dt: float,
    str_save: str,
    tup_snaprange: tuple = None,
):

    mda_msd = MDAnalysis.analysis.msd.EinsteinMSD(
        u = mda_u,
        select = str_select
    )

    if tup_snaprange is None:
        tup_snaprange = (0, len(mda_u.trajectory))

    mda_msd.run(
        start = tup_snaprange[0],
        stop = tup_snaprange[1],
        verbose = True,
    )
    
    int_nframes = mda_msd.n_frames
    np_final = np.zeros( shape=int_nframes, dtype=[('time', 'f4'), ('msd', 'f4')])
    np_final['time'] = np.arange(int_nframes)*float_dt # make the lag-time axis
    np_final['msd'] = mda_msd.results.timeseries
    np.savetxt(
        fname = str_save,
        X = np_final,
        header = ' '.join(np_final.dtype.names)
    )

def hbonds_collect(
    dict_file: dict,
    str_save: str,
) -> None:
    
    with open(list(dict_file.values())[0], 'r') as fp:
        list_header = fp.readline().split()[1:]

    for str_header in list_header:
        print(str_header)
        np_data = np.zeros( shape=(len(dict_file),3) )
        for int_id, int_key in enumerate(dict_file):
            print(f'data: {dict_file[int_key]}')
            np_tmp = np.genfromtxt(dict_file[int_key], names=True)
            np_data[int_id, 0] = int_key
            np_data[int_id, 1] = np_tmp[str_header][0]
            np_data[int_id, 2] = np_tmp[str_header][2]

        np.savetxt(
            fname = f'{str_save}.{str_header}.csv',
            X = np_data,
            header = f'temperature {str_header} error'
        )

def hbonds_ave(
    list_file: str,
    str_save: str,
) -> None:

    str_file = list_file[0]
    np_tmp = np.loadtxt(str_file)
    print(f'data: {str_file}')
    print(np_tmp)
    np_data = np_tmp[:,2:]

    for str_file in list_file[1:]:
        np_tmp = np.loadtxt(str_file)
        print(f'data: {str_file}')
        print(np_tmp)
        np_data = np.append(np_data, np_tmp[:,2:], axis=0)

    np_mean = np.mean(np_data, axis=0)
    np_ssd = np.std(np_data, axis=0, ddof=1)
    np_error = np_ssd/np.sqrt(np_data.shape[0])

    np_final = np.zeros( shape=(3, np_data.shape[1]) )
    np_final[0] = np_mean
    np_final[1] = np_ssd
    np_final[2] = np_error

    with open(str_file, 'r') as fp:
        list_header = fp.readline().split()[3:]
    np.savetxt(
        X = np_final,
        fname = str_save,
        header = ' '.join(list_header)
    )

def hbonds_status(
    str_file: str,
    list_range: list[tuple],
    str_save: str,
) -> np.array:

    print(f'data: {str_file}')
    np_data = np.genfromtxt(str_file, names=True)

    np_save = np.zeros(
        shape=(len(list_range)), 
        dtype=[
            ('rangel', 'i4'), 
            ('ranger', 'i4'), 
            ('nhbonds', 'f4'), 
            ('chain0', 'f4'), 
            ('chain1', 'f4'), 
            ('chain2', 'f4'), 
            ('chain3x', 'f4'), 
            ('chain1x', 'f4'),
            ('chain2x', 'f4'), 
        ]
    )

    for int_id, tup_range in enumerate(list_range):
        np_save['rangel'][int_id] = tup_range[0]
        np_save['ranger'][int_id] = tup_range[1]
        np_save['nhbonds'][int_id] = np.mean(np_data['nhbonds'][tup_range[0]:tup_range[1]])
        np_save['chain0'][int_id] = np.mean(np_data['chain0'][tup_range[0]:tup_range[1]])
        np_save['chain1'][int_id] = np.mean(np_data['chain1'][tup_range[0]:tup_range[1]])
        np_save['chain2'][int_id] = np.mean(np_data['chain2'][tup_range[0]:tup_range[1]])
        np_save['chain3x'][int_id] = np.mean(np_data['chain3x'][tup_range[0]:tup_range[1]])
        np_save['chain1x'][int_id] = np.mean(np_data['chain1'][tup_range[0]:tup_range[1]]) \
            + np.mean(np_data['chain2'][tup_range[0]:tup_range[1]]) \
            + np.mean(np_data['chain3x'][tup_range[0]:tup_range[1]])
        np_save['chain2x'][int_id] = np.mean(np_data['chain2'][tup_range[0]:tup_range[1]]) \
            + np.mean(np_data['chain3x'][tup_range[0]:tup_range[1]])

    if str_save:
        np.savetxt(
            fname = str_save,
            X = np_save,
            header = ' '.join(np_save.dtype.names)
        )

    return np_save

def hbonds(
    universe,
    hydrogens_sel: str,
    acceptors_sel: str,
    str_save: str,
    donors_sel: str = None,
    update_selections: bool = False,
    float_ave: float = 1.0,
    d_a_cutoff: float = 3.0,
    h_d_a_angle_cutoff: float = 30,
    int_sperate: int = 10000,
) -> None:

    print(str_save)

    mda_hba = HydrogenBondAnalysis(
        universe = universe,
        hydrogens_sel = hydrogens_sel,
        acceptors_sel = acceptors_sel,
        donors_sel = donors_sel,
        d_a_cutoff = d_a_cutoff,
        h_d_a_angle_cutoff = h_d_a_angle_cutoff,
        update_selections = update_selections,
    )

    tup_snaprange = (0, len(universe.trajectory))

    np_final = np.zeros( shape=(tup_snaprange[1], 6) )

    for start in range(tup_snaprange[0], tup_snaprange[1], int_sperate):
        stop = start + int_sperate
        if stop > tup_snaprange[1]:
            stop = tup_snaprange[1]
        print(f'{start} -> {stop}/{tup_snaprange[1]}')
        mda_hba.run(
            start = start,
            stop = stop,
            verbose = True,
        )

        np_final[start:stop, 0] = mda_hba.times
        np_final[start:stop, 1] = mda_hba.count_by_time()/float_ave
        np_final[start:stop, 2] = mda_hba.nchain_count_by_time(0)/float_ave
        np_final[start:stop, 3] = mda_hba.nchain_count_by_time(1)/float_ave
        np_final[start:stop, 4] = mda_hba.nchain_count_by_time(2)/float_ave
        np_final[start:stop, 5] = mda_hba.nchain_count_by_time(3, 10)/float_ave

    np.savetxt(
        fname = str_save,
        X = np_final,
        header = ' '.join(['snap', 'nhbonds', 'chain0', 'chain1', 'chain2', 'chain3x'])
    )

def rdf(
    mda_atomgroup_0,
    mda_atomgroup_1,
    tup_snaprange: tuple,
    str_save: str,
    tup_rrange: tuple,
    int_nbins: int = 200  
) -> None:

    print(str_save)
    
    mda_rdf = InterRDF( 
        mda_atomgroup_0,
        mda_atomgroup_1, 
        nbins = int_nbins,
        range = tup_rrange,
    )

    mda_rdf.run(
        start = tup_snaprange[0],
        stop = tup_snaprange[1],
        verbose = True,
    )
    array_final = np.zeros( shape=(int_nbins), dtype=[('r(Angstrom)', 'f4'), ('RDF', 'f4')])
    array_final['r(Angstrom)'] = mda_rdf.results.bins
    array_final['RDF'] = mda_rdf.results.rdf
    
    np.savetxt(
        fname = str_save,
        X = array_final,
        header = ' '.join(array_final.dtype.names)
    )

