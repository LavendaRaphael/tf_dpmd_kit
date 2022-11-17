from MDAnalysis.analysis.rdf import InterRDF
import numpy as np

def rdf(
    mda_atomgroup_0,
    mda_atomgroup_1,
    list_snaprange: list[tuple],
    str_save: str,
    tup_rrange: tuple,
    int_nbins: int = 200  
) -> None:
    
    mda_rdf = InterRDF( 
        mda_atomgroup_0,
        mda_atomgroup_1, 
        nbins = int_nbins,
        range = tup_rrange,
    )
    for tup_snaprange in list_snaprange:
        str_save = f'{str_save}{tup_snaprange[0]:07d}_{tup_snaprange[1]:07d}.csv'
        print(str_save)

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
