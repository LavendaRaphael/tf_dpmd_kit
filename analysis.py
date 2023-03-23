import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from tf_dpmd_kit import plm
from MDAnalysis.analysis.rdf import InterRDF
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
import MDAnalysis.analysis.msd
from MDAnalysis.analysis.base import AnalysisBase
import pandas as pd
from MDAnalysis.lib.distances import capped_distance, calc_angles, calc_dihedrals
import time

def carbonic_rolling_plt(
    float_xscale: float = 1,
    str_xlabel: str = None,
    tup_ylim: tuple = None,
    str_file: str = 'carbonic_rolling.csv',
    list_header: list = None,
    list_ypos: list = None,
    list_yticklabels: list = None,
):

    df_data = pd.read_csv(str_file)
    print(df_data)
    
    fig, ax = plt.subplots()
    
    df_new = df_data.where(df_data.isnull(), 1)
    df_new['frame'] = df_data['frame']
    
    for idx, header in enumerate(list_header):
        df_data_tmp = df_data[df_data[header].notnull()]
        df_new_tmp = df_new[df_new[header].notnull()]
        if len(df_data_tmp)==0:
            continue
        ax.scatter( df_new_tmp['frame']*float_xscale, df_new_tmp[header]*list_ypos[idx], s=2, edgecolors='none', c='tab:blue', alpha=df_data_tmp[header])

    ax.set_xlabel(str_xlabel)
    ax.set_yticks(list_ypos)
    ax.set_yticklabels(list_yticklabels)
    ax.set_ylim(tup_ylim)

    return fig, ax

def carbonic_state(
    str_file: str = 'carbonic.csv',
    str_save: str = 'carbonic_state.csv',
):
    start = time.time()

    df_data = pd.read_csv(str_file)
    print(df_data)

    df_new = df_data.apply(lambda x: carbonic_evalstate(x['ncarbonyl'], x['noho'], x['dihedral0(rad)'], x['dihedral1(rad)']), axis=1, result_type='expand')
    df_new.columns = ['CO3','0.5','HCO3','1.5','TT','CT','CC','2.5','H3CO3']
    df_new.insert(0, 'frame', df_data['frame'])
    print(df_new)

    df_new.to_csv(str_save, index=False)

    end = time.time()
    print(end-start)

def carbonic_evalstate(
    ncarbonyl, 
    noho,
    alpha,
    beta,
):
    list_re = [None]*9

    if ncarbonyl == 3:
        # 300 CO3
        int_state = 0
    elif ncarbonyl == 2:
        if noho != 0:
            # 201
            int_state = 1
        else:
            # 210 HCO3
            int_state = 2
    elif ncarbonyl == 1:
        int_state = carbonic_conformer(alpha, beta)
        list_re[int_state] = 1
        if noho != 0:
            # 102 111
            int_state = 3
    elif ncarbonyl == 0:
        if noho != 0:
            # 003 012 021
            int_state = 7
        else:
            # 003 H3CO3
            int_state = 8

    list_re[int_state] = 1
    return list_re

def carbonic_conformer(
    alpha,
    beta
):
    pio2 = np.pi/2
    bool_alpha = ((alpha > -pio2) & (alpha < pio2))
    bool_beta = ((beta > -pio2) & (beta < pio2))

    if bool_alpha:
        if bool_beta:
            # TT
            int_state = 4
        else:
            # TC
            int_state = 5
    else:
        if bool_beta:
            # CT
            int_state = 5
        else:
            # CC
            int_state = 6

    return int_state

class Carbonic(AnalysisBase):

    def __init__(
        self,
        carbonic_c,
        carbonic_o,
        atomg_h,
        water_o,
        cutoff = 1.3):

        trajectory = carbonic_c.universe.trajectory
        super(Carbonic, self).__init__(trajectory)

        self.carbonic_c = carbonic_c
        self.carbonic_o = carbonic_o
        self.atomg_h = atomg_h
        self.water_o = water_o
        self.cutoff = cutoff

    def _prepare(self):

        self.results = np.zeros((self.n_frames, 5))
        self.carbonic_c2 = self.carbonic_c[[0,0]]

    def _single_frame(self):

        self.results[self._frame_index, 1:] = self._carbonic()
        self.results[self._frame_index, 0] = self._ts.frame

    def _conclude(self):

        columns = ['frame',
                   'ncarbonyl', 'noho', 'dihedral0(rad)', 'dihedral1(rad)']
        self.df = pd.DataFrame(self.results, columns=columns)
    
    def _carbonic(self):
    
        box = self._ts.dimensions

        list_carbonyl = []
        list_hydroxyl_o = []
        list_hydroxyl_h = []
        list_oho = []
        for o_id in range(3):
            np_id, np_distances = capped_distance(
                reference = self.carbonic_o[[o_id]],
                configuration = self.atomg_h,
                max_cutoff = self.cutoff,
                box = box,
            )
            if np.size(np_distances) == 0:
                list_carbonyl.append(o_id)
                continue
            h_id = np_id[np.argmin(np_distances), 1]
            list_hydroxyl_o.append(o_id)
            list_hydroxyl_h.append(h_id)

            np_id, _ = capped_distance(
                reference = self.atomg_h[[h_id]],
                configuration = self.water_o,
                max_cutoff = self.cutoff,
                box = box,
            )
            if np.size(np_id) != 0:
                list_oho.append(o_id)

        len_carbonyl = len(list_carbonyl)
        len_oho = len(list_oho)

        if len_carbonyl != 1:
            return len_carbonyl, len_oho, None, None

        dihedrals = calc_dihedrals(
            self.carbonic_o[ list_carbonyl*2 ],
            self.carbonic_c2,
            self.carbonic_o[list_hydroxyl_o],
            self.atomg_h[list_hydroxyl_h],
            box=box
        )
        return len_carbonyl, len_oho, dihedrals[0], dihedrals[1]

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
