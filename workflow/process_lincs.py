import numpy as np
from collections import defaultdict
from pathlib import Path
import pandas as pd
import sys
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
from typing import Tuple, List, Dict
import tqdm
import warnings
from cmapPy.pandasGEXpress.parse import parse
from damply import dirs

import utils
import time

warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)

def update_colData(
        compound:str,
        colData:Dict,
        cpd_metadata:pd.DataFrame,
        database:pd.DataFrame
    ) ->int:

    
    
    cpd_subset = cpd_metadata[cpd_metadata['cmap_name']==compound]
    # num_inchis = len(
    inchis = pd.unique(pd.unique(cpd_subset['inchi_key']))
    lincs_targets = list(cpd_subset['target'])
    db_subset = database[database['inchikey'].isin(inchis)]
    # print(db_subset)
    # db_subset=db_subset[db_subset['name']==compound.title()]
    if db_subset.shape[0]==0:
        return -1
    else:
        inchi_key  = list(pd.unique(cpd_subset['inchi_key']))[0]
        lincs_moa  =list(pd.unique(cpd_subset['moa']))[0]

        
        # cpd_info= database[database['inchikey']==inchi_key]
        cid = list(db_subset['cid'].values)[0]
        name = list(db_subset['name'].values)[0]
        inchi = list(db_subset['inchikey'].values)[0]
        colData['cid'].append(cid)
        colData['name'].append(name)
        colData['inchikey'].append(inchi)
        colData['lincs_moa'].append(lincs_moa)
        colData['lincs_targets'].append(lincs_targets)
        colData['lincs_aliases'].append(list(cpd_subset['compound_aliases']))
        colData['adb_aliases'].append(list(db_subset['mapped_name']))
        
    
    
    return cid
    
        
 
def main(
    gene_file:str = "geneinfo_beta.txt",
    cell_file:str = "cellinfo_beta.txt",
    cpd_file:str = "compoundinfo_beta.txt",
    sig_file:str = "siginfo_beta.txt",
    data_file:str = "level5_beta_trt_cp_n720216x12328.gctx",
    adb_file:str = "all_adb_compounds.csv",
    only_named_cells:bool = True,
    only_cancer_cells:bool = True,
    valid_gene_types:List[str]=['landmark'],
    valid_dose_units:List[str]=['uM'],
    only_nonneg_signatures:bool = True
    ) -> None:

    cpd_metadata, signature_data, gene_metadata = utils.load_metadata(
        gene_file = gene_file,
        cell_file = cell_file,
        cpd_file = cpd_file,
        sig_file = sig_file, 
        only_named_cells = only_named_cells,
        only_cancer_cells = only_cancer_cells,
        valid_gene_types = valid_gene_types,
        valid_dose_units = valid_dose_units
    )
    

    cell_names = np.array(signature_data['cell_iname']).reshape(-1,1)
    batch_names = np.array(signature_data['bead_batch']).reshape(-1,1)
    cell_encoder = OneHotEncoder().fit(cell_names)
    batch_encoder = OneHotEncoder().fit(batch_names)
    colData = defaultdict(list)
    
    adb = pd.read_csv(dirs.RAWDATA / adb_file,index_col=0)
    
    missing_list = []
    signature_assay = defaultdict(list)
    
    for compound in (cpd_pbar:=tqdm.tqdm(pd.unique(signature_data['cmap_name']) )):
        cpd_pbar.set_description(f"Working on {compound}")
        compound_specific_signatures = signature_data[signature_data["cmap_name"] == compound]
        compound_signature_ids = list(compound_specific_signatures['sig_id'])
        cid = update_colData(
            compound = compound,
            colData = colData,
            cpd_metadata=cpd_metadata,
            database = adb)
        if cid==-1:
            continue
        
        s = time.time()
        data = parse(
            str(dirs.RAWDATA/ data_file),
            cid=compound_signature_ids,
            rid = gene_metadata['gene_id']
        )
        e = time.time()
        
        cpd_pbar.set_description(f"Working on {compound} - Data Loaded in  {np.round(e-s,3)} seconds.")
        data = data.data_df
        
        signatures = []
        for gene in (gene_pbar:=tqdm.tqdm(data.index,leave=False)): # iterate over genes 
            gene_pbar.set_description(f"Working on gene {gene}")
            signature_vector = utils.compute_gene_specific_signature(
                expression_vector = data.loc[gene],
                signature_info = compound_specific_signatures,
                cell_encoder = cell_encoder,
                batch_encoder = batch_encoder)
            
            signature_vector = np.array(signature_vector).reshape(1,-1)
            signatures.append(signature_vector)
        all_signatures = np.vstack(signatures)
        
        signature_assay[cid] = np.mean(all_signatures,axis=0).flatten()
        # mean_signature = np.mean(all_signatures,axis=0,keepdims=True)
        
    colData = pd.DataFrame(colData)
    signature_assay = pd.DataFrame(signature_assay)
    colData.to_csv(dirs.RESULTS/ "colData.csv")
    signature_assay.to_csv(dirs.RESULTS/ "signatures.csv")
        
            
        
          
           
            
               
                
                    
if __name__=="__main__":
     main()
    
    
       
        
        
        

    





