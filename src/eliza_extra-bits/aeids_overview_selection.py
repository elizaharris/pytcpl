
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
pd.options.mode.chained_assignment = None  # default='warn'

parent_dir = "/data/eliza/pytcpl/"
sys.path.append(parent_dir)

from src.pipeline.pipeline_constants import ROOT_DIR, PROFILER_PATH, LOG_DIR_PATH, RAW_DIR_PATH
from src.pipeline.pipeline_helper import load_config, init_config, query_db, get_metadata

config, config_path = load_config()
global CONFIG
CONFIG = config

#%% Look at the assays and chemicals

assay_component_endpoint = query_db(f"SELECT * FROM assay_component_endpoint;") # lists 1499 assays (unique aeid)
chemicals = query_db(f"SELECT * FROM chemical;") # 9559 unique dsstox substance ids
samples = query_db(f"SELECT * FROM sample;") # 9559 unique dsstox substance ids

chemicals_assays_intersect = np.zeros((chemicals.shape[0],assay_component_endpoint.shape[0])) # rows = chemicals, cols = assays, 0 = not tested, 1 = tested
chemicals_assays_ntests = chemicals_assays_intersect.copy() # as above but number of assays tested per chemical
chemicals_assays_nconcs = chemicals_assays_intersect.copy() # as above but number of dif concs tested per chemical
assay_component_endpoint["n_chemicals"] = -1 # space for the number of chemicals tested per assay
assay_component_endpoint["n_chemicals_no_chid"] = -1 # space for the number of chemicals with no chid (?) tested per assay

if 1: # Run only if recalculating data
    for n, aeid in enumerate(assay_component_endpoint["aeid"]):
        # Get the raw data for an assay endpoint
        print(str(n)+", "+str(aeid))
        global AEID
        AEID = int(aeid)
        path = os.path.join(RAW_DIR_PATH, f"{AEID}{CONFIG['file_format']}")
        if not os.path.exists(path):
            select_cols = ['spid', 'aeid', 'logc', 'resp', 'cndx', 'wllt']
            table_mapping = {'mc0.m0id': 'mc1.m0id', 'mc1.m0id': 'mc3.m0id'}
            join_string = ' AND '.join([f"{key} = {value}" for key, value in table_mapping.items()])
            qstring = f"SELECT {', '.join(select_cols)} " \
                    f"FROM mc0, mc1, mc3 " \
                    f"WHERE {join_string} AND aeid = {AEID};"
            df = query_db(query=qstring)
            df.to_parquet(path, compression='gzip')
        else:
            df = pd.read_parquet(path)
        # Map to the chemical info
        df = df.merge(samples, on='spid', how='left')
        df = df.merge(chemicals, on='chid', how='left')
        n_spid = len(np.unique(df["spid"])) # Number of samples
        n_chid = len(np.unique(df["chid"])) # Number of unique chemicals for this assay (all nans are one category)
        n_chid_nan = sum(np.isnan(np.array(df["chid"]))) # Number of nan chemicals for this assay (eg. mixtures, etc.)
        assay_component_endpoint["n_chemicals"].iloc[n] = n_chid
        assay_component_endpoint["n_chemicals_no_chid"].iloc[n] = n_chid_nan
        for chem in np.unique(df["chid"]):
            if np.isnan(chem):
                continue
            index_this_chem = np.where(chemicals["chid"]==chem)[0][0]
            counts_this_chem = sum(df["chid"]==chem)
            diff_concs_this_chem = len(np.unique(round(df["logc"][df["chid"]==chem],2)))
            chemicals_assays_intersect[index_this_chem,n] = 1 # This chemical - assay combination was tested
            chemicals_assays_ntests[index_this_chem,n] = counts_this_chem 
            chemicals_assays_nconcs[index_this_chem,n] = diff_concs_this_chem 

        # Save after every 25 (in case...)
        if n%25 == 0:
            assay_component_endpoint.to_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/assays_overview.csv'))
            df2 = pd.DataFrame(chemicals_assays_intersect, columns = assay_component_endpoint["aeid"], index=chemicals["chid"])
            df2.to_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_intersect.csv'))
            df2 = pd.DataFrame(chemicals_assays_ntests, columns = assay_component_endpoint["aeid"], index=chemicals["chid"])
            df2.to_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_ntests.csv'))
            df2 = pd.DataFrame(chemicals_assays_nconcs, columns = assay_component_endpoint["aeid"], index=chemicals["chid"])
            df2.to_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_nconcs.csv'))

#%% Plots 

if 0:
    # Read in again
    assay_component_endpoint = pd.read_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/assays_overview.csv')).drop("Unnamed: 0",axis=1)
    chemicals_assays_intersect = pd.read_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_intersect.csv'))
    chemicals_assays_intersect.index = chemicals_assays_intersect["chid"]; 
    chemicals_assays_intersect = chemicals_assays_intersect.drop("chid",axis=1)
    chemicals_assays_ntests = pd.read_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_ntests.csv'))
    chemicals_assays_ntests.index = chemicals_assays_ntests["chid"]; 
    chemicals_assays_ntests = chemicals_assays_ntests.drop("chid",axis=1)
    chemicals_assays_nconcs = pd.read_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/df_chemicals_assays_nconcs.csv'))
    chemicals_assays_nconcs.index = chemicals_assays_nconcs["chid"]; 
    chemicals_assays_nconcs = chemicals_assays_nconcs.drop("chid",axis=1)

    # Add info about cell-based or not
    assays = query_db(f"SELECT * FROM assay;") 
    assay_components = query_db(f"SELECT * FROM assay_component;")
    assay_component_info_full = assay_component_endpoint.merge(assay_components, on='acid', how='left')
    assay_component_info_full = assay_component_info_full.merge(assays, on='aid', how='left')

    # Save full assay info
    assay_component_info_full.to_csv(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/assays_info_full.csv'))

    # Plot bar charts of relevant columns
    for c in assay_component_info_full.columns:
        data = assay_component_info_full[c]
        data = np.array([ np.nan if d is None else d for d in data  ])
        values = np.unique(data)
        # values = values[values != "NA"]
        if len(values)<2:
            print(c+": Only one category")
            continue
        if len(values)>40:
            print(c+": Too many categories (",str(len(values)),")")
            continue
        counts = np.array([ sum(data==v) for v in values ])
        order = np.argsort(-counts)
        fig, ax = plt.subplots()
        x_pos = np.arange(len(values))
        plt.bar(x_pos, counts[order],color = "slateblue")
        plt.xticks(x_pos, values[order], rotation=90,fontsize=8)
        plt.ylabel("Number of assays")
        plt.title(c)
        plt.tight_layout()
        plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/bar_by_'+c+'.png'))  
        plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/bar_by_'+c+'.pdf'))  

    # Plot: Histogram, chemicals per assay
    fig, ax = plt.subplots()
    ax.hist(assay_component_info_full['n_chemicals'], 50, histtype="bar",color = "slateblue",
                                cumulative=False, label="Histogram")
    ax.set_ylabel("n(chemicals per assay)")
    ax2 = ax.twinx()
    ax2.hist(assay_component_info_full['n_chemicals'], 50, histtype="step",color = "midnightblue",
                                cumulative=True, label="Cumulative histogram")
    ax2.set_ylabel("n(chemicals per assay), cumulative")
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/hist_chemicals_per_assay.png'))  
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/hist_chemicals_per_assay.pdf'))  

    # Plot: Histogram, assays per chemical (sum on chemicals_assays_intersect)
    data = np.nansum(chemicals_assays_intersect,axis=1)
    fig, ax = plt.subplots()
    ax.hist(data, 50, histtype="bar",color = "slateblue",cumulative=False, label="Histogram")
    ax.set_ylabel("n(assays per chemical)")
    ax2 = ax.twinx()
    ax2.hist(data, 50, histtype="step",color = "midnightblue",cumulative=True, label="Cumulative histogram")
    ax2.set_ylabel("n(assays per chemical), cumulative")
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/hist_assays_per_chemical.png'))  
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/hist_assays_per_chemical.pdf'))  

    # Colour map: Chemical-assay intersection, ordered
    order_assay = np.argsort(-np.nansum(chemicals_assays_intersect,axis=0))
    order_chemicals = np.argsort(-np.nansum(chemicals_assays_intersect,axis=1))
    data = np.array(chemicals_assays_intersect)[:,order_assay]
    data = data[order_chemicals,:]
    fig, ax = plt.subplots(1,3,figsize=(6,8))
    cmap = colors.ListedColormap(['white', 'blue'])
    ax[0].imshow(data, cmap=cmap)
    ax[0].set_title("Assay-chemical pairs",fontsize=8)
    # Add colour map: Chemical-assay intersection: Colour by number of assays run
    data = np.array(chemicals_assays_ntests)[:,order_assay]
    data = data[order_chemicals,:]
    ax[1].imshow(data, cmap="GnBu",vmin=0,vmax=np.nanmax(data))
    ax[1].set_title("Coloured by n(assays run) \n max = "+str(int(np.nanmax(data))),fontsize=8)
    # Colour map: Chemical-assay intersection: Colour by number of concs tested
    data = np.array(chemicals_assays_nconcs)[:,order_assay]
    data = data[order_chemicals,:]
    ax[2].imshow(data, cmap="GnBu",vmin=0,vmax=np.nanmax(data))
    ax[2].set_title("Coloured by n(concs tested) \n max = "+str(int(np.nanmax(data))),fontsize=8)
    for i in np.arange(0,3):
        ax[i].set_ylabel("Chemical",fontsize=8)
        ax[i].tick_params(axis='y', labelsize=6)
        ax[i].set_xlabel("Assay",fontsize=8)
        ax[i].tick_params(axis='x', labelsize=6)
    plt.tight_layout()
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/map_assay-chemical-intersect.png'))  
    plt.savefig(os.path.join(ROOT_DIR, '../../figs/aeids_overview_selection/map_assay-chemical-intersect.pdf'))  

#%% Identify assays to be used and save this info

assays_selected_info = assay_component_info_full.copy()

# 1. Filter for cell-based only
filter = assays_selected_info["assay_format_type"] == "cell-based"
print("Filtering for cell-based assays removes "+str(sum(~filter))+" leaving "+str(sum(filter))+"; initial n(assays) = "+str(assay_component_info_full.shape[0]))
assays_selected_info = assays_selected_info[filter]

# 2. Filter for assays with a dtxsid
filter = assays_selected_info["assay_format_type"] == "cell-based"
print("Filtering for cell-based assays removes "+str(sum(~filter))+" leaving "+str(sum(filter))+"; initial n(assays) = "+str(assay_component_info_full.shape[0]))
assays_selected_info = assays_selected_info[filter]

# 3. Filter for assays with the most chemicals
cutoff = np.percentile(assays_selected_info["n_chemicals"],75) 
filter = assays_selected_info["n_chemicals"] > cutoff
print("Filtering for number of chemicals removes "+str(sum(~filter))+" leaving "+str(sum(filter))+"; initial n(assays) = "+str(assay_component_info_full.shape[0]))
assays_selected_info = assays_selected_info[filter]

# Save
assays_selected_info.to_csv(os.path.join(ROOT_DIR, '../../data/assays_selected_info.csv'))
aeids = assays_selected_info["aeid"]
with open(os.path.join(ROOT_DIR, '../../config/aeid_list.in'), 'w') as fp: # Full aeid list for processing at this stage
    for item in aeids:
        fp.write("%s\n" % item) # write each item on a new line
    print('Done')
aeids_short = assays_selected_info["aeid"].iloc[np.round(np.random.uniform(0,len(aeids),3),0)] # Make a short list for troubleshooting also
with open(os.path.join(ROOT_DIR, '../../config/aeid_list_short.in'), 'w') as fp:
    for item in aeids_short:
        fp.write("%s\n" % item) # write each item on a new line
    print('Done')
