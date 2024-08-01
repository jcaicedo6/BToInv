import ROOT
import pandas as pd

def RDF_to_pandas(key:str, inputRootFiles:str, columns:list=None, cut=None, fraction=None):
    """
    :param key: Tree name of RootFiles
    :param inputRootFiles: {path}/files.root
    :param columns: Variables of interest
    :param cut: root-like cut to be applied.
    :param fraction: fraction of sampling
    :return: Pandas dataframe
    """

    if cut: df = ROOT.RDataFrame(key, inputRootFiles).Filter(cut)
    else: df = ROOT.RDataFrame(key, inputRootFiles)

    if columns: df = df.AsNumpy(columns=columns)
    else: df = df.AsNumpy()
    
    if fraction: 
        if (fraction >1):
            pdf =pd.DataFrame(df)
            pdf = pdf.sample(frac=fraction, replace=True, random_state=1)
        else:
            pdf =pd.DataFrame(df)
            pdf = pdf.sample(frac=fraction, random_state=1)
    else: pdf =pd.DataFrame(df)

    return pdf
#------#


#df = RDF_to_pandas(key='b', inputRootFiles='../ntuple.root', columns=['Mbc', 'deltaE'], cut='Mbc > 5.2 && B_p > 0.2')