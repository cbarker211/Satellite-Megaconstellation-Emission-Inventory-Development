#!/usr/bin/python
""" Save out rocket emissions as COARDS compliant files for input to HEMCO/GEOS-Chem.

    Input Files:
        launchfile     (str): The path of the launch database, a netcdf file containing information for each rocket launch.
        reentryfile    (str): The path of the reentry database, a netcdf file containing information for each object reentry.
        rocketinfofile (str): The path of the rocket database, a netcdf file containing information for each rocket.
        peifile        (str): The path of the pei database, a csv file containing primary emission indices of major emission species for each propellant type.

    Returns:
        launch_propellant (csv) :   A file to view the vertical distribution of propellant consumption.
                                        Each launch or reentry adds three lines with len(altitude):
                                            1) Prop consumed in layer (0 for reentries).
                                            2) mid altitude of layer (from box heights).
                                            3) % prop consumed in layer with Ross14 study regridded onto desired model. 
        emis_distribution  (csv):   A file to view the vertical distribution of emissions per species.
                                        Each launch or reentry adds nine lines with len(altitude), containing species emissions in each layer in order of: 
                                            BC, CO, NOx, H2O, Al2O3, Cl, HCl, Cl2, CO2
                                        Last line is mid altitude of layer (from box heights).
        HEMCO Emission Files  (netcdf): Daily emission files with hourly time bins (dim = time, lon, lat, lev), containing the following variables:
                                            rocket_launch_nox
                                            rocket_fuel_nox
                                            rocket_reentry_nox
                                            rocket_h2o
                                            rocket_bc
                                            rocket_launch_al
                                            rocket_reentry_al
                                            rocket_hcl
                                            rocket_cl
                                            rocket_cl2
                                            rocket_co  
                                            rocket_co2
"""

# Import modules:
import argparse

import os
from os import path
from netCDF4 import Dataset
import numpy as np
import pandas as pd
pd.set_option('display.max_colwidth', None)
import matplotlib.pyplot as plt
import geopandas as gpd
import fiona

import sys
sys.path.append('./python_modules/')
from distribute_emis_func import make_grid_LL, read_gc_box_height, get_ross_profiles, interp_prop_mass
from alt_emis_func import calculate_bc_ei, calculate_nox_ei, calculate_co_ei, calculate_cl_ei

class RocketData:
    '''Read rocket launch and re-entry emissions'''
    def __init__(self,launchfile,reentryfile, rocketinfofile, peifile):
        
        self.read_launch_emissions(launchfile)
        self.read_reentry_emissions(reentryfile)
        self.read_rocket_info(rocketinfofile)
        self.define_pei(peifile)
    
    def read_launch_emissions(self, launchfile):
        """Read launch emissions from the database.

        Args:
            launchfile (str): The path of the launch database.
        """        

        #Open Dataset.
        fh = Dataset(launchfile, mode='r')
        
        #Extract variables.
        self.launch_lons    = fh.variables['Longitude'][:]
        self.launch_lats    = fh.variables['Latitude'][:]
        self.launch_time    = fh.variables['Time(UTC)'][:]  # UTC hours
        self.launch_datestr = fh.variables['Date'][:]
        self.launch_id      = fh.variables['COSPAR_ID'][:].astype(str)
        self.launch_rocket  = fh.variables['Rocket_Name'][:]
        self.launch_smc     = fh.variables['Megaconstellation_Flag'][:]
        
        self.launch_year, self.launch_month, self.launch_day = [], [], []
        for datestring in self.launch_datestr:
            self.launch_year.append(np.int64(datestring[:4]))
            self.launch_month.append(np.int64(datestring[4:6]))
            self.launch_day.append(np.int64(datestring[6:]))
            
        #Close Dataset.
        fh.close()
        
    def read_reentry_emissions(self, reentryfile):
        """Read re-entry emissions from the database.

        Args:
            reentryfile (str): The path of the reentry database.
        """        

        # Open Dataset
        fg = Dataset(reentryfile, mode='r')
        
        # Extract variables:
        self.reentry_id           = fg.variables['COSPAR_ID'][:]
        self.reentry_name         = fg.variables['Object_Name'][:]
        self.reentry_category     = fg.variables['Category'][:]
        self.reentry_lats         = fg.variables['Latitude'][:]
        self.reentry_lons         = fg.variables['Longitude'][:]
        self.reentry_datestr      = fg.variables['Date'][:]
        self.reentry_time         = fg.variables['Time (UTC)'][:]
        self.reentry_abl_mass     = fg.variables['Ablatable_Mass'][:]
        self.reentry_abl_deg      = fg.variables['Ablation_Degree'][:]
        self.reentry_per_alu      = fg.variables['Percent_Aluminium'][:]
        self.reentry_other_mass   = fg.variables['Other_Mass'][:]
        self.reentry_smc          = fg.variables['Megaconstellation_Flag'][:]
        
        self.reentry_year, self.reentry_month, self.reentry_day = [], [], []
        for datestring in self.reentry_datestr:
            self.reentry_year.append(np.int64(datestring[:4]))
            self.reentry_month.append(np.int64(datestring[4:6]))
            self.reentry_day.append(np.int64(datestring[6:]))
        
        #Close Dataset.
        fg.close()
              
    def read_rocket_info(self, rocketinfofile):
        """Read information about the propellant mass and type for each rocket stage.
        The stage mass is also included in the database for reentry purposes, but not required here.

        Args:
            rocketinfofile (str): The path of the rocket database.
        """        

        # Open Dataset.
        fh = Dataset(rocketinfofile, mode='r')
        
        # Extract variables:
        self.rocket_name         = fh.variables['Rocket_Name'][:]
        self.booster_prop_mass   = fh.variables['Booster_PropMass'][:]
        self.booster_prop_type   = fh.variables['Booster_Fuel_Type'][:]
        self.stage1_prop_mass    = fh.variables['Stage1_PropMass'][:]
        self.stage1_prop_type    = fh.variables['Stage1_Fuel_Type'][:]
        self.stage2_prop_mass    = fh.variables['Stage2_PropMass'][:]
        self.stage2_prop_type    = fh.variables['Stage2_Fuel_Type'][:]
        self.stage3_prop_mass    = fh.variables['Stage3_PropMass'][:]
        self.stage3_prop_type    = fh.variables['Stage3_Fuel_Type'][:]
        self.stage4_prop_mass    = fh.variables['Stage4_PropMass'][:]
        self.stage4_prop_type    = fh.variables['Stage4_Fuel_Type'][:]
        
        # Close Dataset.
        fh.close()
    
    def define_pei(self, peifile):
        """ Define the primary emission indices.
            Each array should contain 5 values, corresponding to the emission index for :
                Hypergolic, Kerosene, Solid, Hydrogen, Methane

        Args:
            peifile (str): The path of the pei database.
        """        
        
        # Define propellant types.
        self.pei_fuel_type     = np.array(['Hypergolic', 'Kerosene', 'Solid', 'Hydrogen', 'Methane'])
        
        # Load pei values.
        pei_file = pd.read_csv(peifile)
        self.h2o_pei   = pei_file["H2O"].to_numpy() # Each variable must be a numpy array in order to make the np.where routine work later.
        self.h2_pei    = pei_file["H2"].to_numpy()          
        self.co_pei    = pei_file["CO"].to_numpy()
        self.co2_pei   = pei_file["CO2"].to_numpy()
        self.bc_pei    = pei_file["BC"].to_numpy()
        self.nox_pei   = pei_file["NOx"].to_numpy()
        self.al2o3_pei = pei_file["Al2O3"].to_numpy()
        self.cly_pei   = pei_file["Cly"].to_numpy()
              
class OutputEmis:
    def __init__(self, rocketdata, res, ts, levs, booster_levs, meco_levs, sei_levs, seco_levs, 
                 reentry_levs, vert_filepath, months, dataset, year, ground_landing_list,
                 stage_alt_dict, vertical_source, save_diag, save_emis, outdir, model_alt):

        ####################################
        # Initial Setup
        ####################################
        
        # Setup input variables.
        self.launch_id = rocketdata.launch_id # Extract launch IDs from the imported rocket data.
        self.ts = ts
        self.nlev = levs
        self.booster_levs = booster_levs
        self.meco_levs = meco_levs
        self.sei_levs = sei_levs
        self.seco_levs = seco_levs
        self.bot_reenter = min(reentry_levs)
        self.top_reenter = max(reentry_levs)
        self.n_reenter_levs = self.top_reenter - self.bot_reenter + 1
        self.nhours = 24
        self.include_landing = True
        self.dataset = dataset
        self.year = year
        self.ground_landing_list = ground_landing_list
        self.stage_alt_dict = stage_alt_dict
        self.total_landing_prop = 0
        self.included_prop = 0
        self.missing_prop_total = 0
        self.model_alt = model_alt
        
        # Define Molecular Weights:
        self.mw_h   = 1.008
        self.mw_h2  = self.mw_h * 2
        self.mw_h2o = self.mw_h2 + 16.00

        # Get Ross vertical profiles of propellant burned on a fine grid.
        # The smallest GEOS-Chem grid box has a height of 0.129 km, so setting find grid res to 0.1 km.
        # Setting maximum to 100 km, this is the top of the highest bin in the Ross distribution. 
        self.ross_alt_edge, self.ross_prop_mass, self.ross_cumulative_mass = get_ross_profiles()
        self.fine_grid_bot_alt = np.arange(0,100,0.1)*1e3
        self.fine_grid_mid_alt = np.arange(0.05,100,0.1)*1e3
        self.fine_grid_top_alt = np.arange(0.1,100.1,0.1)*1e3
        self.prop_in_fine_grid, self.fine_grid_mass = interp_prop_mass(self.fine_grid_bot_alt, self.fine_grid_mid_alt, self.fine_grid_top_alt, 
                                                                         self.ross_alt_edge, self.ross_cumulative_mass)

        # Create variables for totals for a future sanity check:
        self.nox_launch_total, self.bc_total, self.h2o_total = 0.0, 0.0, 0.0
        self.co_total,self.al2o3_launch_total, self.cl_total, self.co2_total = 0.0, 0.0, 0.0, 0.0
        self.hcl_total, self.cl2_total, self.nox_reentry_total, self.al2o3_reentry_total = 0.0, 0.0, 0.0, 0.0
        self.mass_survive = 0
        self.missing_emis = np.zeros(9)
        self.booster_prop_consumed, self.stage1_prop_consumed, self.stage2_prop_consumed = np.zeros(len(self.launch_id)), np.zeros(len(self.launch_id)), np.zeros(len(self.launch_id))
        self.total_prop_consumed = np.zeros((len(rocket_data.h2o_pei),3))
        self.launch_count = 0
        
        # Make grid at input resolution (no bounds supplied so full global).
        # Can use this with different in and out edges.
        if res == "2x25":
            grid = make_grid_LL("2x2.5")  
        else:
            grid = make_grid_LL(res)  
        #Extract useful variables from the grid.    
        self.nlon = len(grid['lon'])
        self.nlat = len(grid['lat'])
        self.lon =  grid['lon']
        self.lat =  grid['lat']
        
        # Create the arrays for the output file to check prop/species distribution.
        self.csv_count, self.csv_count_2 = 0, 0
        if self.dataset == 1:
            launch_length  = np.sum(np.where((rocket_data.launch_year == year)&(rocket_data.launch_smc==False), 1, 0))
            reentry_length = np.sum(np.where((rocket_data.reentry_year == year)&(rocket_data.reentry_smc==False), 1, 0))
        elif self.dataset == 2:
            launch_length  = np.sum(np.where((rocket_data.launch_year == year)&(rocket_data.launch_smc==True), 1, 0))
            reentry_length = np.sum(np.where((rocket_data.reentry_year == year)&(rocket_data.reentry_smc==True), 1, 0))
        elif self.dataset == 3:
            launch_length  = np.sum(np.where((rocket_data.launch_year == year), 1, 0))
            reentry_length = np.sum(np.where((rocket_data.reentry_year == year), 1, 0))
        self.output_csv_launch_prop  = np.zeros((launch_length*3,self.nlev))
        self.output_csv_emis  = np.zeros(((launch_length+reentry_length)*10,self.nlev))
        
        # Build the ocean landing list for the year in question.
        ocean_landings = raul_data[(raul_data["Name"].str.contains("ASDS")) & (raul_data["Description"].str.contains("Landing -"))].reset_index(drop=True)
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planned")].reset_index(drop=True)
        ocean_landings = ocean_landings[~ocean_landings["Name"].str.contains("planed")].reset_index(drop=True)
        date_col = []
        for landing in range(len(ocean_landings["Name"])):
            year_ind = ocean_landings["Description"][landing].index("Landing -")
            datestr = ocean_landings["Description"][landing][year_ind:].replace(" ","")[8:17].replace("-","")
            date_col.append(datestr)
        ocean_landings["Date"] = date_col
        ocean_landings['Date'] = pd.to_datetime(ocean_landings['Date'],format='%d%b%y')
        self.ocean_landings = ocean_landings[(ocean_landings["Date"] < f'{year+1}-01-01') & (ocean_landings["Date"] >= f'{year}-01-01')].reset_index(drop=True)
        
        #################################
        # Loop over each day and month.
        #################################
        #for m in range(0,12):
        for m in months:
            
            #Setup days in month
            if m in [0,2,4,6,7,9,11]:
                ndays = 31
            elif m in [3,5,8,10]:
                ndays = 30
            elif m == 1:
                ndays = 29
            else:
                sys.exit("Invalid month selected.")

            # Get month string:
            if m+1<=9:
                self.strmon='0'+str(m+1)
            else:
                self.strmon=str(m+1)

            print('MONTH = ', self.strmon)             
            
            #Loop over a single day
            #for d in range(6,9):
            
            #Loop over all days:
            for d in range(ndays):
                
                # Get day string:
                if d+1<=9:
                    self.strday='0'+str(d+1)
                else:
                    self.strday=str(d+1)

                if vertical_source == "StateMet":
                    self.gc_file = path.join(vert_filepath, f'{self.year}/GEOSChem.StateMet.{self.year}{self.strmon}{self.strday}_0000z.nc4')
                    file_exists = os.path.exists(self.gc_file)
                    if not file_exists: continue
                    self.box_height, self.area, self.bot_alt, self.mid_alt, self.top_alt = read_gc_box_height(self.gc_file,self.nlev)
                elif vertical_source == "custom":
                    layer_data = np.loadtxt("./input_files/" + vert_filepath ,delimiter=",")
                    self.top_alt = np.zeros((self.nlev,self.nlat,self.nlon))
                    self.bot_alt = np.zeros((self.nlev,self.nlat,self.nlon))
                    self.mid_alt = np.zeros((self.nlev,self.nlat,self.nlon))
                    self.area = np.zeros((self.nlat,self.nlon))
                    for lat in range(self.nlat):
                        for lon in range(self.nlon):
                            self.bot_alt[:,lat,lon] = layer_data[::-1,0]*1000
                            self.mid_alt[:,lat,lon] = layer_data[::-1,1]*1000
                            self.top_alt[:,lat,lon] = layer_data[::-1,2]*1000
                            self.area[lat,lon]      = 1 # This is set to 1 as its only needed to output the final netcdf files, not the stats.
                
                # Set up the data arrays which will be saved to file.
                self.rocket_launch_nox  = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_fuel_nox    = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_reentry_nox = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_h2o         = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_bc          = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_co          = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_co2         = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_launch_al   = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_reentry_al  = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_hcl         = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_cl          = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))
                self.rocket_cl2         = np.zeros((self.nhours,self.nlev,self.nlat,self.nlon))

                # Initialize:
                self.pmin,self.pmax=np.nan,np.nan
                self.qmin,self.qmax=np.nan,np.nan

                # Define output file name for this day
                if self.dataset == 1:
                    l_ind=np.where((rocket_data.launch_month==(m+1))&(rocket_data.launch_day==np.int64(d+1))&(rocket_data.launch_year==self.year)&(~np.isnan(rocket_data.launch_time))&(rocket_data.launch_smc==False))[0]
                    r_ind=np.where((rocket_data.reentry_month==(m+1))&(rocket_data.reentry_day==np.int64(d+1))&(rocket_data.reentry_year==self.year)&(~np.isnan(rocket_data.reentry_time))&(rocket_data.reentry_smc==False))[0]
                    self.outfile_rocket_emis = path.join(outdir,self.strmon,f'byproduct_emis_nonsmc_{res}_{str(self.nlev)}L_{self.year}{self.strmon}{self.strday}.nc4')
                elif self.dataset == 2:
                    l_ind=np.where((rocket_data.launch_month==(m+1))&(rocket_data.launch_day==np.int64(d+1))&(rocket_data.launch_year==self.year)&(~np.isnan(rocket_data.launch_time))&(rocket_data.launch_smc==True))[0]
                    r_ind=np.where((rocket_data.reentry_month==(m+1))&(rocket_data.reentry_day==np.int64(d+1))&(rocket_data.reentry_year==self.year)&(~np.isnan(rocket_data.reentry_time))&(rocket_data.reentry_smc==True))[0]
                    self.outfile_rocket_emis = path.join(outdir,self.strmon,f'byproduct_emis_smc_{res}_{str(self.nlev)}L_{self.year}{self.strmon}{self.strday}.nc4')
                elif self.dataset == 3:
                    l_ind=np.where((rocket_data.launch_month==(m+1))&(rocket_data.launch_day==np.int64(d+1))&(rocket_data.launch_year==self.year)&(~np.isnan(rocket_data.launch_time)))[0]
                    r_ind=np.where((rocket_data.reentry_month==(m+1))&(rocket_data.reentry_day==np.int64(d+1))&(rocket_data.reentry_year==self.year)&(~np.isnan(rocket_data.reentry_time)))[0]
                    self.outfile_rocket_emis = path.join(outdir,self.strmon,f'byproduct_emis_all_{res}_{str(self.nlev)}L_{self.year}{self.strmon}{self.strday}.nc4')

                ########################################################################
                # Call grid_emis function to calculate distribution from launches
                ########################################################################
                
                if len(l_ind)>0:
                    #print('----> Day processing: ',self.strday,f"{len(l_ind)} Launches.")
                    self.grid_emis(l_ind,
                                   rocket_data.launch_lons[l_ind],
                                   rocket_data.launch_lats[l_ind],
                                   rocket_data.launch_time[l_ind],
                                   'launch',
                                   rocket_data.launch_id[l_ind],
                                   rocket_data.launch_rocket[l_ind],
                                   d,
                                   m,
                                   )

                ########################################################################
                # Call grid_emis function to calculate distribution from re-entries
                ########################################################################              
                
                if len(r_ind)>0:
                    #print('----> Day processing: ',self.strday,f"{len(r_ind)} Re-entries.")
                    self.grid_emis(r_ind,
                                   rocket_data.reentry_lons[r_ind],
                                   rocket_data.reentry_lats[r_ind],
                                   rocket_data.reentry_time[r_ind],
                                   'reentry',
                                   ' ',
                                   ' ',
                                   d,
                                   m,
                                   )

                ####################################
                # Save file as NetCDF
                ####################################
                
                if save_emis == True:
                    if len(l_ind)>0 or len(r_ind)>0:
                        # If there are launches or reentries, then save file normally.
                        self.save_netcdf()
                    else:
                        # If there are no launches, create a dummy file to accommodate HEMCO file reading.
                        lat_ind=np.argmin(abs(self.lat-0))
                        self.qmin,self.qmax=int(lat_ind),int(lat_ind)
                        lon_ind=np.argmin(abs(self.lon-0))
                        self.pmin,self.pmax=int(lon_ind),int(lon_ind)
                        self.save_netcdf()
                       
        ####################################
        # Conservation of mass check
        ####################################
        
        # At the end of the year, check that the emissions going out to the netcdf files are the same as those from the fine grid.
        # This is done by comparing the emis_distribution output (equivalent to netcdf without the conversion to kgm-2s-1) and the totals.
        final_emis = np.asarray([self.bc_total, self.co_total, self.nox_reentry_total+self.nox_launch_total, self.h2o_total, 
                                 self.al2o3_launch_total+self.al2o3_reentry_total, self.cl_total, self.hcl_total, self.cl2_total, self.co2_total])
        
        spec_names = ["BC", "CO", "NOx", "H2O", "Al2O3", "Cl", "HCl", "Cl2", "CO2"]
        diff_out = np.zeros(9)
        for event in range(int(np.shape(self.output_csv_emis)[0]/10)):
            for spec in range(9):
                diff_out[spec] += np.sum(self.output_csv_emis[event*10+spec,:]) 
        
        for i_emis in range(9):
            with np.errstate(divide='ignore', invalid='ignore'):
                if np.abs(diff_out[i_emis]-final_emis[i_emis])/final_emis[i_emis]*100 > 0.1:
                    print(f"Error with {spec_names[i_emis]} emissions - {np.abs(diff_out[i_emis]-final_emis[i_emis])/final_emis[i_emis]*100:.2f}%.") 
        
        df   = pd.DataFrame(self.output_csv_launch_prop)
        df_2 = pd.DataFrame(self.output_csv_emis)
        if save_diag == True:
            df.to_csv(f"./emission_diag/launch_propellant_{self.year}_{dataset}_{res}_{self.nlev}.csv", float_format='%.3f')
            df_2.to_csv(f"./emission_diag/emis_distribution_{self.year}_{dataset}_{res}_{self.nlev}.csv", float_format='%.3f')
             
    def process_launch_event_altitudes(self, valid_index, p, q, launch_rocket,launch_id):
        
        # Define the propellant saved for landing for reusable rockets.
        if launch_id in self.ground_landing_list:
            landing_prop_percent = 100-5.6-5.6-1.2
        else:
            landing_prop_percent = 100-5.6-1.2
                    
        stage_alt_beco = stage_alt_dict[f"{launch_rocket} BECO"] 
        stage_alt_meco = stage_alt_dict[f"{launch_rocket} MECO"]
        stage_alt_sei  = stage_alt_dict[f"{launch_rocket} SEI1"] 
        stage_alt_seco = stage_alt_dict[f"{launch_rocket} SECO"]
        
        if np.isnan(stage_alt_meco) and np.isnan(stage_alt_sei):
            # Define the rocket configuration.
            if rocket_data.booster_prop_type[valid_index] != '':
                if rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] == '' and rocket_data.stage3_prop_type[valid_index] == '' and rocket_data.stage4_prop_type[valid_index] == '':
                    rocket_config_type = 0
                elif rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] == '' and rocket_data.stage4_prop_type[valid_index] == '':
                    rocket_config_type = 0
                elif rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] != '' and rocket_data.stage4_prop_type[valid_index] == '':
                    rocket_config_type = 1
                elif rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] != '' and rocket_data.stage4_prop_type[valid_index] != '':
                    rocket_config_type = 2
                else:
                    sys.exit("Incorrect rocket configuration.")
            else:
                if rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] == '' and rocket_data.stage4_prop_type[valid_index] == '':
                    rocket_config_type = 3
                elif rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] != '' and rocket_data.stage4_prop_type[valid_index] == '':
                    rocket_config_type = 4
                elif rocket_data.stage1_prop_type[valid_index] != '' and rocket_data.stage2_prop_type[valid_index] != '' and rocket_data.stage3_prop_type[valid_index] != '' and rocket_data.stage4_prop_type[valid_index] != '':
                    rocket_config_type = 5
                else:
                    sys.exit("Incorrect rocket configuration.")
                    
            stage_alt_beco = self.booster_levs[rocket_config_type]
            stage_alt_meco = self.meco_levs[rocket_config_type]
            stage_alt_sei  = self.sei_levs[rocket_config_type]    
            stage_alt_seco = self.seco_levs[rocket_config_type]
        
        ###########
        # Boosters
        ###########
        
        # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
        
        if rocket_data.booster_prop_type[valid_index] != '':
            
            # If BECO occurs above the fine grid, then we only use some of the booster emissions.
            if stage_alt_beco * 1e3 > self.fine_grid_top_alt[-1]:
                self.fine_grid_mass_booster = np.asarray(self.fine_grid_mass)
                self.booster_alt_index = None
            # If BECO occurs below the fine grid, then we can use all of the booster emissions.
            else:
                # Find out which layer BECO occurs in.
                self.booster_alt_index = np.argmax(self.fine_grid_top_alt > stage_alt_beco*1e3) + 1
                if not (int(self.fine_grid_bot_alt[self.booster_alt_index-1]) <= stage_alt_beco*1e3 <= self.fine_grid_top_alt[self.booster_alt_index-1]):
                    sys.exit(f"Booster indexing error: {self.fine_grid_bot_alt[self.booster_alt_index-1]} {stage_alt_beco*1e3} {self.fine_grid_top_alt[self.booster_alt_index-1]}")
                fine_grid_mass_booster = np.asarray(self.fine_grid_mass[:self.booster_alt_index])
                
                # Normalize the vertical emission profile for the booster emissions.
                if launch_rocket == "Falcon Heavy" and self.include_landing == True:
                    self.fine_grid_mass_booster = fine_grid_mass_booster / np.sum(fine_grid_mass_booster) * landing_prop_percent
                    self.total_landing_prop += (100 - landing_prop_percent) * rocket_data.booster_prop_mass[valid_index] / 100
                else:
                    self.fine_grid_mass_booster = fine_grid_mass_booster / np.sum(fine_grid_mass_booster) * 100
                
        ###########
        # Stage 1
        ###########              
        # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 100km.           
        if launch_rocket == "GSLV Mk III":
            fei_alt = 46.534
        elif launch_rocket == "Pegasus XL":
            fei_alt = 11.79
        elif launch_rocket == "LauncherOne":
            fei_alt = 12.37
        else:
            fei_alt = 0     
        self.fei_alt_index = np.argmax(np.abs(self.fine_grid_top_alt > fei_alt*1e3))   
        
        # If MECO occurs above the fine grid, then we only use some of the stage 1 emissions. 
        if stage_alt_meco * 1e3 > self.fine_grid_top_alt[-1]:
            self.MECO_alt_index = None
        # If MECO occurs below the fine grid, then we can use all of the stage 1 emissions.
        else:    
            # Find out which layer MECO occurs in.
            self.MECO_alt_index = np.argmax(np.abs(self.fine_grid_top_alt > stage_alt_meco*1e3)) + 1
            if not (int(self.fine_grid_bot_alt[self.MECO_alt_index-1]) <= stage_alt_meco*1e3 <= self.fine_grid_top_alt[self.MECO_alt_index-1]):
                    sys.exit(f"Stage 1 indexing error: {self.fine_grid_bot_alt[self.MECO_alt_index-1]} {stage_alt_meco*1e3} {self.fine_grid_top_alt[self.MECO_alt_index-1]}") 
                      
        self.fine_grid_mass_stage1 = np.asarray(self.fine_grid_mass[self.fei_alt_index:self.MECO_alt_index])
        
        if stage_alt_meco * 1e3 > self.fine_grid_top_alt[-1]:
            if launch_rocket == "GSLV Mk III":
                # Now interpolate the cumulative mass to optimise what percentage of propellant is within the model.
                fei_mass = np.interp(fei_alt,self.ross_alt_edge, self.ross_cumulative_mass)
                model_percent = (self.prop_in_fine_grid - fei_mass) / (100 - fei_mass) * 100
                self.fine_grid_mass_stage1 = self.fine_grid_mass_stage1 / np.sum(self.fine_grid_mass_stage1) * model_percent
        else:
            # If the rocket is Falcon v1.2 Block 5, then MECO is always in model, and 93.2% of the emissions are used for ascent.
            # There is an extra 'boostback' burn for ground landings. No information of % of propellant, so assumed to be equal to the reentry burn.
            if launch_rocket == "Falcon 9 v1.2" and self.include_landing == True:
                self.fine_grid_mass_stage1 = self.fine_grid_mass_stage1 / np.sum(self.fine_grid_mass_stage1) * landing_prop_percent 
                self.total_landing_prop += (100 - landing_prop_percent) * rocket_data.stage1_prop_mass[valid_index] / 100                
            # For all other rockets, 100% is used for ascent.
            else:
                self.fine_grid_mass_stage1 = self.fine_grid_mass_stage1 / np.sum(self.fine_grid_mass_stage1) * 100  
        
        ###########
        # Stage 2
        ###########
                               
        # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 100km.
        # If SEI occurs above the fine grid, then we can ignore the second stage.
        if stage_alt_sei * 1e3 > self.fine_grid_top_alt[-1]:
            self.fine_grid_mass_stage2 = np.asarray([0])
            self.SEI_alt_index = None
            self.seco_alt_index = None
        
        # If SEI is within fine grid, we need to differentiate based on SECO.
        else:
            self.SEI_alt_index = np.argmax(np.abs(self.fine_grid_top_alt > stage_alt_sei*1e3))
            
            # If SECO occurs above the fine grid, then we only use some of the stage 2 emissions.
            if stage_alt_seco * 1e3 > self.fine_grid_top_alt[-1]:
                
                self.fine_grid_mass_stage2 = np.asarray(self.fine_grid_mass[self.SEI_alt_index:])
                # Now interpolate the cumulative mass to optimise what percentage of propellant is within GEOS-Chem.
                sei_mass = np.interp(stage_alt_sei,self.ross_alt_edge, self.ross_cumulative_mass)
                model_percent = (self.prop_in_fine_grid - sei_mass) / (100 - sei_mass) * 100
                self.fine_grid_mass_stage2 = self.fine_grid_mass_stage2 / np.sum(self.fine_grid_mass_stage2) * model_percent
                self.seco_alt_index = None
            
            # If SECO occurs below the fine grid, then we can use all of the stage 2 emissions.
            else:
                self.seco_alt_index = np.argmax(np.abs(self.fine_grid_top_alt > stage_alt_seco*1e3)) + 1
                self.fine_grid_mass_stage2 = np.asarray(self.fine_grid_mass[self.SEI_alt_index:self.seco_alt_index])
                self.fine_grid_mass_stage2 = self.fine_grid_mass_stage2 / np.sum(self.fine_grid_mass_stage2) * 100
            
        ###########
        # Stage 3
        ###########
                
        # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.       
        # TEI data not collected here, would need to check if any other rocket third stages start below 100 km.        
        # Calculate third stage emissions for Minotaur 1 only.    
        if launch_rocket == "Minotaur 1":
            # Minotaur 1 third stage ignition occurs below the model boundary, so we can need to add some of the stage 3 emissions.
            self.TEI_alt_index = np.argmax(np.abs(self.fine_grid_top_alt > 73.47*1e3))
            self.fine_grid_mass_stage3 = np.asarray(self.fine_grid_mass[self.TEI_alt_index:])
            # Now interpolate the cumulative mass to optimise what percentage of propellant is within GEOS-Chem.
            tei_mass = np.interp(73.47,self.ross_alt_edge, self.ross_cumulative_mass)
            model_percent = (self.prop_in_fine_grid - tei_mass) / (100 - tei_mass) * 100
            self.fine_grid_mass_stage3 = self.fine_grid_mass_stage3 / np.sum(self.fine_grid_mass_stage3) * model_percent
         
        ##################
        # Falcon Reusable
        ##################
           
        # Finally create an extra stage just for the landing emissions of Falcon 9 v1.2. 
        # NOTE: Need to add boostback emissions if top is above 80.
        if launch_rocket in ["Falcon 9 v1.2","Falcon Heavy"]:
            # 5.6% are used in the entry burn, over 70-54.7 km.
            self.entry_top = np.argmax(np.abs(self.fine_grid_top_alt > 70000))
            self.entry_bot = np.argmax(np.abs(self.fine_grid_top_alt > 54700))
            self.fine_grid_mass_entry = np.asarray(self.fine_grid_mass[self.entry_bot:self.entry_top+1])
            self.fine_grid_mass_entry = self.fine_grid_mass_entry / np.sum(self.fine_grid_mass_entry) * 5.6
            
            # 1.2% are used in the landing burn, over 3.3-0 km.
            self.landing_top = np.argmax(np.abs(self.fine_grid_top_alt > 3300))
            self.fine_grid_mass_landing = np.asarray(self.fine_grid_mass[:self.landing_top+1])
            self.fine_grid_mass_landing = self.fine_grid_mass_landing / np.sum(self.fine_grid_mass_landing) * 1.2
        else:
            self.fine_grid_mass_landing = 0
            self.fine_grid_mass_entry = 0
        
        return stage_alt_beco,stage_alt_meco,stage_alt_sei 
    
    def calc_emis(self,start_ind,stop_ind,pei_index,prop_mass,vertical_profile,time_index, q, p, kg_to_kgm2s, total_vertical_propellant,stage):
        
        # Calculate the emissions over the range of the stag within the fine grid (0-100 km).
        ei_bc = calculate_bc_ei(self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.bc_pei[pei_index])
        ei_co, ei_co2 = calculate_co_ei(self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.co_pei[pei_index], rocket_data.co2_pei[pei_index])
        ei_sec_nox = calculate_nox_ei(self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3)
        ei_h2o =  rocket_data.h2o_pei[pei_index] + rocket_data.h2_pei[pei_index] * self.mw_h2o / self.mw_h2 
        ei_cl, ei_hcl, ei_cl2 = calculate_cl_ei(self.fine_grid_mid_alt[start_ind:stop_ind]*1e-3, rocket_data.cly_pei[pei_index])
        bc_emis         = ei_bc * prop_mass * 1e-2 * vertical_profile
        co_emis         = ei_co * prop_mass * 1e-2 * vertical_profile
        co2_emis        = ei_co2 * prop_mass * 1e-2 * vertical_profile
        launch_nox_emis = ei_sec_nox * prop_mass * 1e-2 * vertical_profile
        fuel_nox_emis   = rocket_data.nox_pei[pei_index] * prop_mass * 1e-2 * vertical_profile
        h2o_emis        = ei_h2o * prop_mass * 1e-2 * vertical_profile
        al2o3_emis      = rocket_data.al2o3_pei[pei_index] * prop_mass * 1e-2 * vertical_profile
        cl_emis         = ei_cl * prop_mass * 1e-2 * vertical_profile
        hcl_emis        = ei_hcl * prop_mass * 1e-2 * vertical_profile
        cl2_emis        = ei_cl2 * prop_mass * 1e-2 * vertical_profile
        
        # Place the emissions into a larger array covering the whole fine grid (0-100km). 
        emis_full = np.zeros((len(self.fine_grid_mid_alt),11))
        emis_full[start_ind:stop_ind,0]  = bc_emis
        emis_full[start_ind:stop_ind,1]  = co_emis
        emis_full[start_ind:stop_ind,2]  = co2_emis
        emis_full[start_ind:stop_ind,3]  = launch_nox_emis
        emis_full[start_ind:stop_ind,4]  = fuel_nox_emis
        emis_full[start_ind:stop_ind,5]  = h2o_emis
        emis_full[start_ind:stop_ind,6]  = al2o3_emis
        emis_full[start_ind:stop_ind,7]  = cl_emis
        emis_full[start_ind:stop_ind,8]  = hcl_emis
        emis_full[start_ind:stop_ind,9]  = cl2_emis
        emis_full[start_ind:stop_ind,10] = vertical_profile
        
        vertical_propellant = np.zeros(len(self.mid_alt[:,q,p]))
        selected_alts = []
        bc_prop_initial = np.sum(total_vertical_propellant[:,1])
        # Regrid the vertical_profile to the desired model profile.
        for i, alt in enumerate(self.mid_alt[:,q,p]):
            if i == 0:
                bot_ind = 0
            else:
                bot_ind = np.argmin(np.abs(self.fine_grid_mid_alt - self.bot_alt[i,q,p])) + 1
            top_ind = np.argmin(np.abs(self.fine_grid_mid_alt - self.top_alt[i,q,p])) + 1
            selected_alts.extend(np.arange(bot_ind,top_ind))
            
            # Sum the emissions in this range for each species and place in an array.
            self.rocket_bc[time_index,i,q,p]         += (np.sum(emis_full[bot_ind:top_ind,0]) * 1e-3 / kg_to_kgm2s)
            self.rocket_co[time_index,i,q,p]         += (np.sum(emis_full[bot_ind:top_ind,1]) * 1e-3 / kg_to_kgm2s)
            self.rocket_co2[time_index,i,q,p]        += (np.sum(emis_full[bot_ind:top_ind,2]) * 1e-3 / kg_to_kgm2s)
            self.rocket_launch_nox[time_index,i,q,p] += (np.sum(emis_full[bot_ind:top_ind,3]) * 1e-3 / kg_to_kgm2s)
            self.rocket_fuel_nox[time_index,i,q,p]   += (np.sum(emis_full[bot_ind:top_ind,4]) * 1e-3 / kg_to_kgm2s)
            self.rocket_h2o[time_index,i,q,p]        += (np.sum(emis_full[bot_ind:top_ind,5]) * 1e-3 / kg_to_kgm2s)
            self.rocket_launch_al[time_index,i,q,p]  += (np.sum(emis_full[bot_ind:top_ind,6]) * 1e-3 / kg_to_kgm2s)
            self.rocket_cl[time_index,i,q,p]         += (np.sum(emis_full[bot_ind:top_ind,7]) * 1e-3 / kg_to_kgm2s)
            self.rocket_hcl[time_index,i,q,p]        += (np.sum(emis_full[bot_ind:top_ind,8]) * 1e-3 / kg_to_kgm2s)
            self.rocket_cl2[time_index,i,q,p]        += (np.sum(emis_full[bot_ind:top_ind,9]) * 1e-3 / kg_to_kgm2s)
            total_vertical_propellant[i,0]           += np.sum(emis_full[bot_ind:top_ind,10]) * prop_mass
            total_vertical_propellant[i,1]           += np.sum(emis_full[bot_ind:top_ind,0]) 
            total_vertical_propellant[i,2]           += np.sum(emis_full[bot_ind:top_ind,1]) 
            total_vertical_propellant[i,3]           += np.sum(emis_full[bot_ind:top_ind,3:5]) 
            total_vertical_propellant[i,4]           += np.sum(emis_full[bot_ind:top_ind,5]) 
            total_vertical_propellant[i,5]           += np.sum(emis_full[bot_ind:top_ind,6]) 
            total_vertical_propellant[i,6]           += np.sum(emis_full[bot_ind:top_ind,7]) 
            total_vertical_propellant[i,7]           += np.sum(emis_full[bot_ind:top_ind,8]) 
            total_vertical_propellant[i,8]           += np.sum(emis_full[bot_ind:top_ind,9]) 
            total_vertical_propellant[i,9]           += np.sum(emis_full[bot_ind:top_ind,2]) 
           
        if len(list(set(selected_alts))) != len(selected_alts):
            sys.exit("Error in fine grid indexing.")
            
        self.bc_total           += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,0])       
        self.co_total           += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,1])         
        self.co2_total          += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,2])   
        self.nox_launch_total   += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,3:5])     
        self.h2o_total          += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,5])   
        self.al2o3_launch_total += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,6])   
        self.cl_total           += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,7])   
        self.hcl_total          += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,8])   
        self.cl2_total          += np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,9])   
        self.included_prop      += (np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10]) * prop_mass * 1e-2)
        
        if stage == 5:
            self.missing_prop[stage] += np.round((np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10])),2)
        else:
            self.missing_prop[stage] += np.round((100-np.sum(emis_full[selected_alts[0]:selected_alts[-1]+1,10])),2)
        
        return total_vertical_propellant
    
    def calc_missing_emis(self, pei_index, prop_mass, percent_included):
        ei_bc = calculate_bc_ei([self.model_alt], rocket_data.bc_pei[pei_index])
        ei_co, ei_co2 = calculate_co_ei([self.model_alt], rocket_data.co_pei[pei_index], rocket_data.co2_pei[pei_index])
        ei_sec_nox = calculate_nox_ei([self.model_alt])
        ei_h2o =  rocket_data.h2o_pei[pei_index] + rocket_data.h2_pei[pei_index] * self.mw_h2o / self.mw_h2 
        ei_cl, ei_hcl, ei_cl2 = calculate_cl_ei([self.model_alt], rocket_data.cly_pei[pei_index]) 
        
        t_bc_emis = ei_bc * prop_mass * 1e-2 * percent_included
        t_co_emis = ei_co * prop_mass * 1e-2 * percent_included
        t_co2_emis = ei_co2 * prop_mass * 1e-2 * percent_included
        t_launch_nox_emis = ei_sec_nox * prop_mass * 1e-2 * percent_included
        t_fuel_nox_emis = rocket_data.nox_pei[pei_index] * prop_mass * 1e-2 * percent_included
        t_h2o_emis = ei_h2o * prop_mass * 1e-2 * percent_included
        t_al2o3_emis = rocket_data.al2o3_pei[pei_index] * prop_mass * 1e-2 * percent_included
        t_cl_emis = ei_cl * prop_mass * 1e-2 * percent_included
        t_hcl_emis = ei_hcl * prop_mass * 1e-2 * percent_included
        t_cl2_emis = ei_cl2 * prop_mass * 1e-2 * percent_included
        
        self.missing_prop_total += prop_mass * 1e-2 * percent_included
        self.missing_emis[0] += np.sum(t_bc_emis) 
        self.missing_emis[1] += np.sum(t_co_emis)
        self.missing_emis[2] += np.sum(t_launch_nox_emis)
        self.missing_emis[2] += np.sum(t_fuel_nox_emis)
        self.missing_emis[3] += np.sum(t_h2o_emis)
        self.missing_emis[4] += np.sum(t_al2o3_emis)
        self.missing_emis[5] += np.sum(t_cl_emis)
        self.missing_emis[6] += np.sum(t_hcl_emis)
        self.missing_emis[7] += np.sum(t_cl2_emis)
        self.missing_emis[8] += np.sum(t_co2_emis)            
                  
    def grid_emis(self,index,lon,lat,hour,emis_type,launch_id,launch_rocket,d,m):
        """Grid the data onto the GEOS-Chem horizontal and vertical grid"""
        
        #Loop over each launch/reentry.
        for w in range(len(lon)):
            
            ################
            # Setup grid.
            ################
            # Get grid horizontal indices (p is lon index; q is lat index).
            # This works out the nearest latitude and longitude on the grid.
            p,q = np.argmin(abs(lon[w]-self.lon)), np.argmin(abs(lat[w]-self.lat))
            
            # Usually, any reentries from launches are placed at the same lat/lon as the launch itself.
            # This code is only used where a non launch related reentry occurs on the same day as a launch.
            # It updates the grid range (lat and lon edges) to account for the new location and to prepare for saving to file:
            if np.isnan(self.pmin):
                self.pmin,self.pmax=p,p
                self.qmin,self.qmax=q,q
            if ~np.isnan(self.pmin) and self.pmin>p:
                self.pmin=p
            if ~np.isnan(self.pmax) and self.pmax<p:
                self.pmax=p
            if ~np.isnan(self.qmin) and self.qmin>q:
                self.qmin=q
            if ~np.isnan(self.qmax) and self.qmax<q:
                self.qmax=q
                                
            # Get factor to convert from kg to kg/m2/s:
            kg_to_kgm2s = self.area[q,p]*self.ts
            
            # Format the time of the event.
            time_index = int(hour[w]*60*60)//self.ts
            
            #########################################################
            # Work out vertical distribution for launch emissions.
            #########################################################

            if emis_type=='launch':
                
                ############################################
                # Interpolate propellant mass distribution.
                ############################################

                # Linearly interpolate Ross profile proportion of propellant mass to GEOS-Chem vertical grid.
                # This is just to double check that the propellant consumption profile is calculated correctly.
                self.propellant_in_model, self.gc_relative_mass = interp_prop_mass(self.bot_alt[:,q,p], self.mid_alt[:,q,p], self.top_alt[:,q,p],
                                                                                    self.ross_alt_edge, self.ross_cumulative_mass)

                ############################################
                # Check the rocket type 
                ############################################
                
                rocket_found = False
                valid_index = 0
                
                # Loop over rockets to find matching rocket and associated propellant mass:
                for ll in range(len(rocket_data.rocket_name)):
                    if launch_rocket[w] == rocket_data.rocket_name[ll]:
                        rocket_found = True
                        valid_index=ll
                        #print(f"Rocket name is {launch_rocket[w]}")
                        break
                    
                # Catch any missing rocket propellant masses:
                if rocket_found == False:
                    sys.exit(f"No propellant mass found for {launch_rocket[w]}. Exiting programme.")
                    
                ############################################################################
                # Start processing the launch event altitudes and deal with failed launches.
                ############################################################################
                
                if launch_id[w] in ['2020-F04','2021-F04']:
                    print(f'Launch {launch_id[w]} failed. Not processing.')
                    self.csv_count += 3
                    self.csv_count_2 += 10
                    continue
                                                    
                # Next find index of fuel type in primary emissions index data:
                pei_booster_index = np.where( rocket_data.pei_fuel_type == rocket_data.booster_prop_type[valid_index] )[0]
                pei_stage1_index  = np.where( rocket_data.pei_fuel_type == rocket_data.stage1_prop_type[valid_index] )[0]
                pei_stage2_index  = np.where( rocket_data.pei_fuel_type == rocket_data.stage2_prop_type[valid_index] )[0]
                pei_stage3_index  = np.where( rocket_data.pei_fuel_type == rocket_data.stage3_prop_type[valid_index] )[0]
                pei_stage4_index  = np.where( rocket_data.pei_fuel_type == rocket_data.stage4_prop_type[valid_index] )[0]
                
                stage_alt_beco,stage_alt_meco,stage_alt_sei = self.process_launch_event_altitudes(valid_index, p, q, launch_rocket[w],launch_id[w])
                # NOTE: Need to recheck which failures are within model limits if wanting to run for different vertical heights above 80km.
                # Most failures are for upper stages, and so can be treated as normal here.
                # Full information is provided in source_info/failed_launch_info.txt.
                if launch_id[w] in ['2020-F02','2020-F05',
                                      '2021-F02','2021-F07','2021-F08',
                                      '2022-F01','2022-F02']:
                    self.fine_grid_mass_stage2 = np.asarray([0])
                    self.SEI_alt_index = None
                elif launch_id[w] in ['2020-F07','2021-F01','2021-F07']:
                    if launch_id[w] == '2020-F07': # Astra Rocket 3 launch, rocket shut off at 0.9km.
                       failed_alt = 900
                    elif launch_id[w] == '2021-F01': # Shuang Quxian-1 launch, rocket disintegrated at Max-Q. Approximating altitude using Proton-M and Minotaur-IV max-q alts.
                        failed_alt = 10700
                    elif launch_id[w] == '2021-F07': # Astra Rocket 3 launch, rocket shut off at 31km.
                        failed_alt = 31000
                    cutoff_ind = np.argmin(np.abs(failed_alt - self.fine_grid_mid_alt))
                    self.fine_grid_mass_stage1[:cutoff_ind+1] = np.zeros((len(self.fine_grid_mass_stage1[:cutoff_ind+1])))
                
                ##################################################
                # Sanity Checks for Propellant Mass Distributions
                ##################################################
                
                # Each time, use 0.001% as a the maximum error from rounding / floating point errors.
                error_lim = 0.001
                
                # For all rockets with boosters, the propellant consumed for the boosters should never be bigger than the total booster propellant.
                # Suppressed for Falcon Heavy, as this has reusable boosters. 
                with np.errstate(divide='ignore', invalid='ignore'):
                    if rocket_data.booster_prop_type[valid_index] != '':
                        self.booster_prop_consumed[self.launch_count] = np.sum(self.fine_grid_mass_booster * rocket_data.booster_prop_mass[valid_index] * 1e-2)
                        if stage_alt_meco < self.model_alt:
                            if ((np.abs(self.booster_prop_consumed[self.launch_count] - rocket_data.booster_prop_mass[valid_index]) / rocket_data.booster_prop_mass[valid_index] * 100) > error_lim) and (rocket_data.rocket_name[valid_index] != "Falcon Heavy"):
                                print(np.abs(self.booster_prop_consumed[self.launch_count]),rocket_data.booster_prop_mass[valid_index],np.abs(self.booster_prop_consumed[self.launch_count] - rocket_data.booster_prop_mass[valid_index]) / rocket_data.booster_prop_mass[valid_index] * 100)
                                sys.exit(f"Error with booster emissions -2. {launch_id[w]}") 
                            
                    self.stage1_prop_consumed[self.launch_count]  = np.sum(self.fine_grid_mass_stage1 * rocket_data.stage1_prop_mass[valid_index] * 1e-2)
                    self.stage2_prop_consumed[self.launch_count]  = np.sum(self.fine_grid_mass_stage2 * rocket_data.stage2_prop_mass[valid_index] * 1e-2)

                    # For all rockets, the propellant consumed for each stage should never be bigger than the total propellant in each stage.
                    if ((self.stage1_prop_consumed[self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100 ) > error_lim:
                        sys.exit(f"Error with Stage 1 emissions - 1.")
                    if (rocket_data.rocket_name[valid_index] != "Long March (CZ) 5B"):
                        if (((self.stage2_prop_consumed[self.launch_count] - rocket_data.stage2_prop_mass[valid_index]) / rocket_data.stage2_prop_mass[valid_index] * 100 ) > error_lim):
                            sys.exit("Error with Stage 2 emissions.")

                    # When MECO occurs in the model, the consumed propellant should be within 1% of the total propellant mass of stage 1.  
                    # The error is suppressed for 2020-F07, 2021-F01, and 2021-F07 where rockets failed early. 
                    # Also suppressed for Falcon 9, as this has a reusable first stage.
                    # A check that the Falcon landing distribution is no greater than 7% of the stage 1 emissions is undertaken later in the grid_emis function.
                    if stage_alt_meco < self.model_alt:
                        if (launch_id[w] not in ['2020-F07','2021-F01','2021-F07']) and rocket_data.rocket_name[valid_index] != "Falcon 9 v1.2" and ((np.abs(self.stage1_prop_consumed[self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100) > 1):
                            print(launch_id[w],self.stage1_prop_consumed[self.launch_count],rocket_data.stage1_prop_mass[valid_index],(np.abs(self.stage1_prop_consumed[self.launch_count] - rocket_data.stage1_prop_mass[valid_index]) / rocket_data.stage1_prop_mass[valid_index] * 100))
                            sys.exit("Error with Stage 1 emissions - 2.")  
                        
                self.total_prop_consumed[pei_booster_index,0] += self.booster_prop_consumed[self.launch_count]
                self.total_prop_consumed[pei_stage1_index,1]  += self.stage1_prop_consumed[self.launch_count]
                self.total_prop_consumed[pei_stage2_index,2]  += self.stage2_prop_consumed[self.launch_count]                  

                total_prop_mass = rocket_data.booster_prop_mass[valid_index] + rocket_data.stage1_prop_mass[valid_index] + rocket_data.stage2_prop_mass[valid_index]
                total_prop_mass += rocket_data.stage3_prop_mass[valid_index] + rocket_data.stage4_prop_mass[valid_index]
                
                ##############################################
                # Calculate the emissions for each species.
                ##############################################     
                
                # Creating a 2d array for the prop output.
                # Total prop, bc, co, nox, h2o, al2o3, cl, hcl, cl2, co2
                total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10))    
                self.missing_prop = np.zeros(6)  
                                
                # Check whether there is a booster:
                if ( rocket_data.booster_prop_type[valid_index] != '' ):
                    if np.sum(1e-2 * self.fine_grid_mass_booster) > 1.01:
                        sys.exit("Error with Boosters. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(None,
                                   self.booster_alt_index,
                                   pei_booster_index,
                                   rocket_data.booster_prop_mass[valid_index],
                                   self.fine_grid_mass_booster,
                                   time_index, 
                                   q, 
                                   p, 
                                   kg_to_kgm2s,
                                   total_vertical_propellant,
                                   0)
                        
                # Every rocket has a first stage.
                if np.sum(1e-2 * self.fine_grid_mass_stage1) > 1.01:
                    sys.exit("Error with Stage 1. Propellant distribution exceeds unity.")  
                total_vertical_propellant = self.calc_emis(self.fei_alt_index,
                                self.MECO_alt_index,
                                pei_stage1_index,
                                rocket_data.stage1_prop_mass[valid_index],
                                self.fine_grid_mass_stage1,
                                time_index, 
                                q, 
                                p, 
                                kg_to_kgm2s,
                                total_vertical_propellant,
                                1)

                # Check whether there is a second stage:
                if rocket_data.stage2_prop_type[valid_index] != '' and self.SEI_alt_index != None:
                    if np.sum(1e-2 * self.fine_grid_mass_stage2) > 1.01:
                        sys.exit("Error with Stage 2. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(self.SEI_alt_index,
                                self.seco_alt_index,
                                pei_stage2_index,
                                rocket_data.stage2_prop_mass[valid_index],
                                self.fine_grid_mass_stage2,
                                time_index, 
                                q, 
                                p, 
                                kg_to_kgm2s,
                                total_vertical_propellant,
                                2)
                    
                # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
                # Check for more rockets that have third stage emissions within model.    
                # Add third stage emissions for Minotaur 1.
                if rocket_data.rocket_name[valid_index] == "Minotaur 1":
                    if np.sum(1e-2 * self.fine_grid_mass_stage3) > 1.01:
                        sys.exit("Error with Stage 3 for Minotaur 1. Propellant distribution exceeds unity.")
                    total_vertical_propellant = self.calc_emis(self.TEI_alt_index,
                                None,
                                pei_stage3_index,
                                rocket_data.stage3_prop_mass[valid_index],
                                self.fine_grid_mass_stage3,
                                time_index, 
                                q, 
                                p, 
                                kg_to_kgm2s,
                                total_vertical_propellant,
                                3)
                    
                # If the rocket is a Falcon 9, then add the landing emissions. Kerosene, so no Al2O3 or Cly.     
                if rocket_data.rocket_name[valid_index] in ["Falcon 9 v1.2","Falcon Heavy"] and self.include_landing == True:
                    if launch_id[w] in self.ground_landing_list:
                        falcon_p = p
                        falcon_q = q
                    else:
                        matching = self.ocean_landings[self.ocean_landings["Date"].isin([f"{self.year}-{self.strmon}-{self.strday}"])].reset_index(drop=True)
                        # There is a typo in the databse, a 2023 launch is listed as 2022.      
                        if (matching.shape[0] == 1) or (f"{self.year}-{self.strmon}-{self.strday}" == "2022-04-27"):
                            falcon_lat = np.array(matching["geometry"].y)[0]
                            falcon_lon = np.array(matching["geometry"].x)[0]
                            falcon_p, falcon_q = np.argmin(abs(falcon_lon-self.lon)), np.argmin(abs(falcon_lat-self.lat))
                        # There are two launches on the same day.
                        elif launch_id[w] == "2022-124":
                            falcon_lat = np.array(matching["geometry"].y)[0]
                            falcon_lon = np.array(matching["geometry"].x)[0]
                            falcon_p, falcon_q = np.argmin(abs(falcon_lon-self.lon)), np.argmin(abs(falcon_lat-self.lat))
                        elif launch_id[w] == "2022-125":
                            falcon_lat = np.array(matching["geometry"].y)[1]
                            falcon_lon = np.array(matching["geometry"].x)[1]
                            falcon_p, falcon_q = np.argmin(abs(falcon_lon-self.lon)), np.argmin(abs(falcon_lat-self.lat))
                        elif matching.shape[0] == 0:
                            # The database hasn't been well updated for 2022, so lets just fill in based on most common geolocation for all other 2020-2022 launches.
                            if lon[w] == -81.0 and lat[w] == 28.5:
                                falcon_p = np.argmin(abs(-75-self.lon))
                                falcon_q = np.argmin(abs(32-self.lat))
                            elif lon[w] == -120.6 and lat[w] == 34.7:
                                falcon_p = np.argmin(abs(-122.5-self.lon))
                                falcon_q = np.argmin(abs(30-self.lat))
                            else:
                                sys.exit("Launch not from assigned site.")
                            print(f"Warning: Allocating Falcon Stage 1 landing to lonxlat {self.lon[falcon_p]}x{self.lat[falcon_q]} for {launch_id[w]}.") 
                        elif matching.shape[0] > 1: 
                            sys.exit(f"Multiple ocean entries for Falcon Stage 1 landing for {launch_id[w]}.")
                        else:
                            sys.exit(f"Problem geolocating Falcon Stage 1 landing for {launch_id[w]}") 
                        #print(self.lon[falcon_p],self.lat[falcon_q])
                            
                        if falcon_p > self.pmax:
                            self.pmax = falcon_p
                            #print(f"Falcon9 stage 1 reentry longitude out of bounds for {launch_id[w]}. Updating bounds.")
                        if falcon_q > self.qmax:
                            self.qmax = falcon_q
                            #print(f"Falcon9 stage 1 reentry latitude out of bounds for {launch_id[w]}. Updating bounds.")                  
                    
                    if np.sum(self.fine_grid_mass_entry) + np.sum(self.fine_grid_mass_landing) > 7:
                        sys.exit("Error with Stage 1 Ocean Landing. Propellant distribution exceeds what is expected.")
                        
                    if rocket_data.rocket_name[valid_index] == "Falcon 9 v1.2":
                        pei = pei_stage1_index
                        prop_mass = rocket_data.stage1_prop_mass[valid_index]
                    elif rocket_data.rocket_name[valid_index] == "Falcon Heavy":
                        pei = pei_booster_index
                        prop_mass = rocket_data.booster_prop_mass[valid_index]
                    
                    # NOTE: This section needs to be tweaked if wanting to run for different vertical heights above 80km.
                    # Should implement the boostback burn if the vertical height is increased to 100 km.
                    # First the entry emissions.
                    total_vertical_propellant = self.calc_emis(self.entry_bot,
                                self.entry_top+1,
                                pei,
                                prop_mass,
                                self.fine_grid_mass_entry,
                                time_index, 
                                falcon_q, 
                                falcon_p, 
                                self.area[falcon_q,falcon_p]*self.ts,
                                total_vertical_propellant,
                                5)
                    # Now the landing emissions.
                    total_vertical_propellant = self.calc_emis(None,
                                self.landing_top+1,
                                pei,
                                prop_mass,
                                self.fine_grid_mass_landing,
                                time_index, 
                                falcon_q, 
                                falcon_p, 
                                self.area[falcon_q,falcon_p]*self.ts,
                                total_vertical_propellant,
                                5)
                
                # NOTE: This section may need to be tweaked if wanting to run for different vertical heights above 80km. 
                # Calculate the missing emissions from above model.
                
                # Check if there are missing emissions from the boosters stage.
                if rocket_data.booster_prop_type[valid_index] != '' and self.missing_prop[0] > 0:
                    if 100 - np.sum(self.fine_grid_mass_booster) > 100:
                        sys.exit("Error with Booster emissions above model.") 
                    self.calc_missing_emis(pei_booster_index,
                        rocket_data.booster_prop_mass[valid_index],
                        self.missing_prop[0],                    
                        )
            
                # Check if there are missing emissions from the first stage.
                if rocket_data.stage1_prop_type[valid_index] != '' and self.missing_prop[1] > 0 and launch_id[w] not in ['2020-F04','2020-F07',
                                                                                                                         '2021-F01','2021-F07']:
                    if 100 - np.sum(self.fine_grid_mass_stage1) > 100:
                        sys.exit("Error with Stage 1 emissions above model.") 
                    self.calc_missing_emis(pei_stage1_index,
                        rocket_data.stage1_prop_mass[valid_index],
                        self.missing_prop[1] - self.missing_prop[5],                    
                        )
                    
                # Check if there are missing emissions from the second stage.
                if rocket_data.stage2_prop_type[valid_index] != '' and self.missing_prop[2] > 0 and launch_id[w] not in ['2020-F02','2020-F04','2020-F05','2020-F07',
                                                                                                                                               '2021-F01','2021-F02','2021-F07','2021-F08',
                                                                                                                                               '2022-F01','2022-F02']:
                    if 100 - np.sum(self.fine_grid_mass_stage2) > 100:
                        sys.exit("Error with Stage 2 emissions above model.") 
                    
                    if launch_id[w] == "2022-F03":
                        included_emis = 240/315*100 - (100-self.missing_prop[2])
                    else:
                        included_emis = self.missing_prop[2]
                        
                    self.calc_missing_emis(pei_stage2_index,
                        rocket_data.stage2_prop_mass[valid_index],
                        included_emis,                    
                        )
                    
                # Check if there are missing emissions from the third stage.
                if rocket_data.stage3_prop_type[valid_index] != '' and launch_id[w] not in ['2020-F02','2020-F03','2020-F05','2020-F06',
                                                                                            '2021-F01','2021-F06','2021-F10',
                                                                                            '2022-F02','2022-F05','2022-F07']:
                    if launch_id[w] == "2021-F09":
                        included_emis = 475/521*100
                    elif rocket_data.rocket_name[valid_index] == "Minotaur 1":
                        included_emis = self.missing_prop[3]
                    else:
                        included_emis = 100
                        
                    self.calc_missing_emis(pei_stage3_index,
                        rocket_data.stage3_prop_mass[valid_index],
                        included_emis,                    
                        )
                
                # Check if there are missing emissions from the fourth stage.
                if rocket_data.stage4_prop_type[valid_index] != '' and launch_id[w][-3:] not in ['2020-F08','2020-F09',
                                                                                                 '2021-F01',
                                                                                                 '2022-F02','2022-F04','2022-F05','2022-F07']:
                    self.calc_missing_emis(pei_stage4_index,
                        rocket_data.stage4_prop_mass[valid_index],
                        100,                    
                        )
                
                ##############################################
                # Output the emissions to a file for viewing.
                ##############################################                  
                #with np.errstate(divide='ignore', invalid='ignore'):
                #    total_vertical_propellant[:,0] = total_vertical_propellant[:,0] / total_prop_mass
                self.output_csv_launch_prop[self.csv_count,:] = total_vertical_propellant[:,0]
                self.csv_count += 1
                self.output_csv_launch_prop[self.csv_count,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count += 1  
                self.output_csv_launch_prop[self.csv_count,:] = self.gc_relative_mass
                self.csv_count += 1 
                
                if np.sum(self.gc_relative_mass) == 0 and np.sum(total_vertical_propellant[:,0]) == 0:
                    print(rocket_data.rocket_name[valid_index])
                
                for i in range(1,10):
                    self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
                    self.csv_count_2 += 1 
                self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count_2 += 1

                self.launch_count += 1
                
            #########################################################
            # Work out vertical distribution for reentry emissions.
            #########################################################
            
            if emis_type=='reentry':
                
                total_vertical_propellant = np.zeros((len(self.mid_alt[:,q,p]),10))
                if np.ma.is_masked(rocket_data.reentry_abl_mass[index[w]]) or np.ma.is_masked(rocket_data.reentry_other_mass[index[w]]):
                    pass
                else:
                
                    reentry_ei_al2o3   = rocket_data.reentry_abl_deg[index[w]] * rocket_data.reentry_per_alu[index[w]]
                    if rocket_data.reentry_abl_deg[index[w]] == 0:
                        reentry_ei_nox = 0.175
                    else:
                        reentry_ei_nox = 0.4
                        
                    # Calculate the total mass surviving re-entry in kg
                    if rocket_data.reentry_abl_deg[index[w]] != 0:
                        self.mass_survive += ((rocket_data.reentry_abl_mass[index[w]] + rocket_data.reentry_other_mass[index[w]]) * (1-rocket_data.reentry_abl_deg[index[w]]))
                     
                    # For consistency with launch emissions, the totals are kept in g units. 
                               
                    t_nox_reentry = (rocket_data.reentry_abl_mass[index[w]] + rocket_data.reentry_other_mass[index[w]]) * reentry_ei_nox * 1000
                    self.nox_reentry_total += t_nox_reentry
                    self.rocket_reentry_nox[time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_nox_reentry / self.n_reenter_levs * 1e-3 / kg_to_kgm2s
                    
                    t_al2o3_reentry = rocket_data.reentry_abl_mass[index[w]] * reentry_ei_al2o3 * 1000
                    self.al2o3_reentry_total += t_al2o3_reentry
                    self.rocket_reentry_al[time_index,self.bot_reenter:self.top_reenter+1,q,p] += t_al2o3_reentry / self.n_reenter_levs * 1e-3 / kg_to_kgm2s
                
                    total_vertical_propellant[self.bot_reenter:self.top_reenter+1,3]   += np.full((self.n_reenter_levs),t_nox_reentry/self.n_reenter_levs)
                    total_vertical_propellant[self.bot_reenter:self.top_reenter+1,5]   += np.full((self.n_reenter_levs),t_al2o3_reentry/self.n_reenter_levs)
                    
                for i in range(1,10):
                    self.output_csv_emis[self.csv_count_2,:] = total_vertical_propellant[:,i]
                    self.csv_count_2 += 1
                self.output_csv_emis[self.csv_count_2,:] = self.mid_alt[:,q,p]*1e-3
                self.csv_count_2 += 1
                     
    def save_netcdf(self):
        """Save rocket emissions as COARDS compliant NetCDF files"""
        
        ###################
        # Set up range.
        ###################
        # Change (shift) indices on edges of grid:
        if self.pmin==0: self.pmin=1
        if self.pmax==0: self.pmax=1
        if self.pmax==self.nlon: self.pmax=self.nlon-1
        
        if self.qmin==0: self.qmin=1
        if self.qmin-1==0: self.qmin=1
        if self.qmax==0: self.qmax=1
        if self.qmax==self.nlat: self.qmax=self.nlat-1
        
        emis_initial =    np.sum(self.rocket_launch_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_fuel_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_reentry_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_h2o[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_bc[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_co[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_co2[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_launch_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_reentry_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_hcl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_cl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                        + np.sum(self.rocket_cl2[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                                                
        updated_coords = False
        if (self.qmax+2)>=self.nlat:
            #print('Latitude index out of bounds: updating')  
            updated_coords = True
            if ( self.qmax==self.qmin ):
                self.qmax,self.qmin=self.nlat-2,self.nlat-2
            else:
                self.qmax=self.nlat-2  
                          
        if (self.pmax+2)>=self.nlon:
            #print('Longitude index out of bounds: updating')
            updated_coords = True
            if ( self.pmax==self.pmin ):
                self.pmax,self.pmin=self.nlon-2,self.nlon-2
            else:
                self.pmax=self.nlon-2   
        
        if updated_coords == True:  
            emis_final =    np.sum(self.rocket_launch_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_fuel_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_reentry_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_h2o[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_bc[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_co[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_co2[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_launch_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_reentry_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_hcl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_cl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
                            + np.sum(self.rocket_cl2[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]) \
            
            with np.errstate(divide='ignore', invalid='ignore'):                                                
                if np.abs(emis_initial-emis_final) / emis_initial * 100 > 0.01:
                    sys.exit('Error: index correction removes relevant data')                

        # Check that output range isn't out-of-bounds:
        if self.qmin-1<0 or self.pmin-1<0:
            print('index is less than 0')
        if self.qmax+1>self.nlat:
            print('lat index exceeds no. of lats')
        if self.pmax+1>self.nlon:
            print('lon index exceeds no. of lons')

        ########################################
        # Open Dataset and set up dimensions.
        ########################################
        
        # define dataset:
        ncout = Dataset(self.outfile_rocket_emis, mode='w', format='NETCDF4')

        # Set dimensions (lat, lon, time, lev):
        # Lat and lon data are trimmed to just where the launch occurs.
        ncout.createDimension('lat', self.qmax-self.qmin+3)
        ncout.createDimension('lon', self.pmax-self.pmin+3)
        ncout.createDimension('lev', self.nlev)
        ncout.createDimension('time', self.nhours)
        
        ######################################
        # Create variables for the dimensions.
        ######################################

        # Create individual variables:
        # Time:
        time = ncout.createVariable('time', np.int64, ('time',))
        time.long_name = 'time'
        time.units = 'hours since '+str(self.year)+'-'+self.strmon+'-'+self.strday+' 00:00:00'
        time.calendar = 'standard'
        time.axis = 'T'
        time[:] = np.arange(0,self.nhours)

        # Longitude:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.long_name = 'Longitude'
        lon.units = 'degrees_east'
        lon.axis = 'X'
        lon[:] = self.lon[self.pmin-1:self.pmax+2]

        # Latitude:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.long_name = 'Latitude'
        lat.units = 'degrees_north'
        lat.axis = 'Y'
        lat[:] = self.lat[self.qmin-1:self.qmax+2]
        
        # Level:
        lev = ncout.createVariable('lev', np.int64, ('lev',))
        lev.long_name = 'GEOS-Chem level on '+str(self.nlev)+'-layer model grid'
        lev.units = 'level'
        lev.positive = 'up'
        lev.axis = 'Z'
        lev[:] = np.arange(0, self.nlev)+1

        ######################################
        # Create variables for each species.
        ######################################

        # (1a) Launch NOx:
        rocket_launch_nox = ncout.createVariable('launch_nox_thermal', np.float32, ('time','lev','lat','lon',))
        rocket_launch_nox.long_name = 'Indirect launch emissions of NOx'
        rocket_launch_nox.short_name = 'Indirect launch NOx'
        rocket_launch_nox.units = 'kg/m2/s'
        rocket_launch_nox[:] = self.rocket_launch_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]
        
        # (1b) Fuel NOx:
        rocket_fuel_nox = ncout.createVariable('launch_nox_fuel', np.float32, ('time','lev','lat','lon',))
        rocket_fuel_nox.long_name = 'Direct launch emissions of NOx from propellant combustion'
        rocket_fuel_nox.short_name = 'Direct launch NOx'
        rocket_fuel_nox.units = 'kg/m2/s'
        rocket_fuel_nox[:] = self.rocket_fuel_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]
        
        # (1c) Reentry NOx:
        rocket_reentry_nox = ncout.createVariable('reentry_nox', np.float32, ('time','lev','lat','lon',))
        rocket_reentry_nox.long_name = 'Reentry NOx emissions'
        rocket_reentry_nox.short_name = 'Reentry NOx emissions'
        rocket_reentry_nox.units = 'kg/m2/s'
        rocket_reentry_nox[:] = self.rocket_reentry_nox[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (2) H2O:
        rocket_h2o = ncout.createVariable('launch_h2o', np.float32, ('time','lev','lat','lon',))
        rocket_h2o.long_name = 'Launch H2O emissions'
        rocket_h2o.short_name = 'Launch H2O emissions'
        rocket_h2o.units = 'kg/m2/s'
        rocket_h2o[:] = self.rocket_h2o[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (3) BC:
        rocket_bc = ncout.createVariable('launch_bc', np.float32, ('time','lev','lat','lon',))
        rocket_bc.long_name = 'Launch BC emissions'
        rocket_bc.short_name = 'Launch BC emissions'
        rocket_bc.units = 'kg/m2/s'
        rocket_bc[:] = self.rocket_bc[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (4a) Launch Alumina:
        rocket_launch_al = ncout.createVariable('launch_al', np.float32, ('time','lev','lat','lon',))
        rocket_launch_al.long_name = 'Launch Al2O3 emissions'
        rocket_launch_al.short_name = 'Launch Al2O3 emissions'
        rocket_launch_al.units = 'kg/m2/s'
        rocket_launch_al[:] = self.rocket_launch_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]
        
        # (4b) Reentry Alumina:
        rocket_reentry_al = ncout.createVariable('reentry_al', np.float32, ('time','lev','lat','lon',))
        rocket_reentry_al.long_name = 'Reentry alumina emissions'
        rocket_reentry_al.short_name = 'Reentry alumina emissions'
        rocket_reentry_al.units = 'kg/m2/s'
        rocket_reentry_al[:] = self.rocket_reentry_al[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (5a) HCl:
        rocket_hcl = ncout.createVariable('launch_hcl', np.float32, ('time','lev','lat','lon',))
        rocket_hcl.long_name = 'Launch HCl emissions'
        rocket_hcl.short_name = 'Launch HCl emissions'
        rocket_hcl.units = 'kg/m2/s'
        rocket_hcl[:] = self.rocket_hcl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (5b) Cl:
        rocket_cl = ncout.createVariable('launch_cl', np.float32, ('time','lev','lat','lon',))
        rocket_cl.long_name = 'Launch Cl emissions'
        rocket_cl.short_name = 'Launch Cl emissions'
        rocket_cl.units = 'kg/m2/s'
        rocket_cl[:] = self.rocket_cl[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]
        
        # (5c) Cl2:
        rocket_cl2 = ncout.createVariable('launch_cl2', np.float32, ('time','lev','lat','lon',))
        rocket_cl2.long_name = 'Launch Cl2 emissions'
        rocket_cl2.short_name = 'Launch Cl2 emissions'
        rocket_cl2.units = 'kg/m2/s'
        rocket_cl2[:] = self.rocket_cl2[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2]

        # (6) CO:
        rocket_co = ncout.createVariable('launch_co', np.float32, ('time','lev','lat','lon',))
        rocket_co.long_name = 'Launch CO emissions'
        rocket_co.short_name = 'Launch CO emissions'
        rocket_co.units = 'kg/m2/s'
        rocket_co[:] = self.rocket_co[:,:,self.qmin-1:self.qmax+2,self.pmin-1:self.pmax+2] 

        # Close the file:
        ncout.close()

def check_total_emissions(year,dataset,save,res,levels):
    """ Print total emissions of each species.
        The total emissions including all afterburning are compared to a scenario where only primary emission indices are used.
    """
    total_missing_emis = np.sum(emis_data.missing_emis)*1e-9 
    total_inc_emis = np.sum([emis_data.bc_total*1e-9, emis_data.nox_launch_total*1e-9,   emis_data.h2o_total*1e-9, 
                                    emis_data.co_total*1e-9, emis_data.co2_total*1e-9, emis_data.al2o3_launch_total*1e-9, 
                                    emis_data.hcl_total*1e-9, emis_data.cl_total*1e-9, emis_data.cl2_total*1e-9, 
                                    emis_data.nox_reentry_total*1e-9, emis_data.al2o3_reentry_total*1e-9])
    percent_missing_emis = total_missing_emis / (total_inc_emis + total_missing_emis) * 100
    percent_inc_emis     = total_inc_emis / (total_inc_emis + total_missing_emis) * 100
    total_missing_prop   = emis_data.missing_prop_total*1e-6
    total_inc_prop       = emis_data.included_prop*1e-6
    percent_missing_prop = total_missing_prop / (total_inc_prop + total_missing_prop) * 100
    percent_inc_prop     = total_inc_prop / (total_inc_prop + total_missing_prop) * 100  
    reusable_prop        = np.round(emis_data.total_landing_prop * 1e-6 / (total_inc_prop + total_missing_prop) * 100,2)   
     
    data = {'Species'        : ['BC', 'NOx (Launch)', 'H2O', 'CO', 'CO2', 'Al2O3 (Launch)', 'HCl', 'Cl', 'Cl2', 'NOx (Reentry)', 'Al2O3 (Reentry)',
                                'Total Emis', 'Total Percent Emis','Total Prop', 'Total Percent Prop','Surviving Mass', 'Percent Reusable Propellant'],
            'Emissions 0-80 km [Gg]': [emis_data.bc_total*1e-9, emis_data.nox_launch_total*1e-9,   emis_data.h2o_total*1e-9, 
                                        emis_data.co_total*1e-9, emis_data.co2_total*1e-9, emis_data.al2o3_launch_total*1e-9, emis_data.hcl_total*1e-9, 
                                        emis_data.cl_total*1e-9, emis_data.cl2_total*1e-9, emis_data.nox_reentry_total*1e-9, 
                                        emis_data.al2o3_reentry_total*1e-9, total_inc_emis, percent_inc_emis, 
                                        total_inc_prop, percent_inc_prop,emis_data.mass_survive * 1e-6, reusable_prop],
            'Emissions >80 km [Gg]': [emis_data.missing_emis[0] * 1e-9, emis_data.missing_emis[2] * 1e-9, emis_data.missing_emis[3] * 1e-9, 
                                        emis_data.missing_emis[1] * 1e-9, emis_data.missing_emis[8] * 1e-9, emis_data.missing_emis[4] * 1e-9, emis_data.missing_emis[5] * 1e-9,
                                        emis_data.missing_emis[6] * 1e-9, emis_data.missing_emis[7] * 1e-9, 0, 0, 
                                        total_missing_emis,percent_missing_emis, total_missing_prop, percent_missing_prop,0,0]}
    df = pd.DataFrame(data)
    print(df.round(4))  
    if save == True:  
        df.to_csv(f"./emission_diag/emis_stats_{year}_{dataset}_{res}_{levels}.csv",sep=',')
    
# Main section of the program
if __name__ == "__main__":

    # Configure the script arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-sm', "--start_month", default = "1", choices=str(np.arange(1,13)), help='Start Month (will override final month if greater than final month).')
    parser.add_argument('-fm', "--final_month", default = "12", choices=str(np.arange(1,13)), help='Final Month.')
    parser.add_argument('-sd', "--start_dataset", default = "1", choices=str(np.arange(1,4)), help='Dataset. 1=Non-SMC, 2=SMC, 3=All')
    parser.add_argument('-fd', "--final_dataset", default = "3", choices=str(np.arange(1,4)), help='Dataset. 1=Non-SMC, 2=SMC, 3=All')
    parser.add_argument('-sy', "--start_year", default = "2020", choices=str(np.arange(2020,2023)), help='Start Year.')
    parser.add_argument('-fy', "--final_year", default = "2022", choices=str(np.arange(2020,2023)), help='Final Year.')
    parser.add_argument('-diag', "--save_diagnostics", action='store_true', help='Save diagnostic info.')
    parser.add_argument('-emis', "--save_emissions", action='store_true', help='Save emissions.')
    args = parser.parse_args()

    # Sort out the year range.
    start_year, final_year = int(args.start_year), int(args.final_year)
    if start_year > final_year:
        final_year = start_year + 1
    year_range = np.arange(start_year,final_year+1)  
    
    # Sort out the month range.
    start_month, final_month = int(args.start_month), int(args.final_month)
    if start_month > final_month:
        final_month = start_month + 1
    months = np.arange(start_month-1,final_month)
    
    # Sort out the dataset range.
    start_dataset, final_dataset = int(args.start_dataset), int(args.final_dataset)
    if start_dataset > final_dataset:
        final_dataset = start_dataset + 1
    dataset_range = np.arange(start_dataset,final_dataset+1)
    
    print(f"Years: {start_year}-{final_year}. Months: {start_month}-{final_month}. Datasets: {start_dataset}-{final_dataset}")
      
    #################################
    # Define launch event altitudes.
    #################################
    
    # BECO, MECO and SEI values from literature sources are used wherever possible. 
    # When not available, the average of other rockets with the same configuration is used.

    stage_alt_dict = {}
    stage_alt_rockets = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=str,skip_header=1,usecols=[0],delimiter=",")
    stage_alt_data = np.genfromtxt("./input_files/launch_event_altitudes.csv",dtype=np.float64,skip_header=1,usecols=[1,2,3,4],delimiter=",")
    for i in range(len(stage_alt_data)):
        if stage_alt_data[i][0] == '':
            stage_alt_dict[f"{stage_alt_rockets[i]} BECO"] = None
        else:
            stage_alt_dict[f"{stage_alt_rockets[i]} BECO"] = np.float64(stage_alt_data[i][0])
        if stage_alt_data[i][1] == '':
            stage_alt_dict[f"{stage_alt_rockets[i]} MECO"] = None
        else:
            stage_alt_dict[f"{stage_alt_rockets[i]} MECO"] = np.float64(stage_alt_data[i][1])
        if stage_alt_data[i][2] == '':
            stage_alt_dict[f"{stage_alt_rockets[i]} SEI1"] = None
        else:
            stage_alt_dict[f"{stage_alt_rockets[i]} SEI1"] =  np.float64(stage_alt_data[i][2])
        if stage_alt_data[i][2] == '':
            stage_alt_dict[f"{stage_alt_rockets[i]} SECO"] = None
        else:
            stage_alt_dict[f"{stage_alt_rockets[i]} SECO"] =  np.float64(stage_alt_data[i][3])
    
    # B+1/2S, B+3S, B+4S, 2S, 3S, 4S,    
    booster_cutoff_alts = [66,55,29,0,0,0]
    MECO_alts           = [220,120,64,90,56,52]
    SEI_alts            = [229,120,64,103,61,59]
    SECO_alts           = [356,232,216,312,176,149]
        
    ground_landing_list = ["2020-016","2020-059","2020-086","2020-101", "2021-059", "2022-002","2022-008","2022-009",
                           "2022-040","2022-057","2022-063","2022-144","2022-166","2022-168","2022-173","2022-179"]
    
    ######################################   
    # Define resolutions and set timings.
    ###################################### 
    
    with open("./input_files/define_grid_resolution_timestep.txt","r") as input_file:
        inputs = []
        for line in input_file:
            if line[0] != "#":
                inputs.append(line.strip())
    
    settings = ["timestep","lat_res","lon_res","vertical_levels","vertical_extent","vertical_source","vertical_path"]
    outputs = []
    for setting in settings:
        setting_list = [s for s in inputs if setting in s]
        if len(setting_list) > 1:
            sys.exit("Duplicate settings found")
        else:
            outputs.append(setting_list[0][setting_list[0].find(":")+1:])
    try:
        timestep = int(outputs[0])
        levels   = int(outputs[3])
        model_alt = int(outputs[4])
    except:
        sys.exit("timestep, vertical_levels, and vertical_extent must be integers.")    
    if (model_alt < 0) or (model_alt > 100):
        print("Model altitude must be within 0-100")  
    grid_res  = f'{outputs[1]}x{outputs[2]}'  
    vertical_source = outputs[5]
    if vertical_source == "StateMet":
        if levels not in [47,72]:
            sys.exit("StateMet option requires 47 or 72 levels.")
        vert_filepath = outputs[6] + f"{grid_res}/{levels}/"
    elif vertical_source == "custom":
        vert_filepath = f"geoschem_vertical_grid_{levels}.csv"
    else:
        sys.exit("Invalid vertical_source.")  
        
    print(f"Timestep: {timestep}s. Resolution: {grid_res}x{levels}. Vertical Range: 0-{model_alt} km using {vertical_source} grid.")
    
    # Re-entry emissions are placed in the GEOS-Chem model at ~68 km (reentry burnup location).
    # The emissions are evenly distributed across two vertical layers in the model.
    # 46-47 for the 47-layer model, 65-72 for the 72-layer model (equivalent to 60-80 km).
    # This means greater smoothing for the 47-layer model.
    # WARNING: Changing vertical grid will require adjusting this.
    if levels==47:
        reentry_levs = [45,46] 
    elif levels==72:
        reentry_levs = [64,71]
    else:
        sys.exit("Invalid number of levels.")
        
    ################
    # Import files.
    ################
    
    fiona.drvsupport.supported_drivers['KML'] = 'rw'
    raul_data = gpd.read_file('./databases/reentry/General_SpaceX_Map_Raul.kml', driver='KML', layer =2) # Falcon landing data.   
    launch_path       = './databases/launch_activity_data_2020-2022.nc'
    rocket_info_path  = './databases/rocket_attributes_2020-2022.nc'
    reentry_path      = './databases/reentry_activity_data_2020-2022.nc'
    pei_path          = './input_files/primary_emission_indices.csv'  
    
    #Loop over all years and run functions depending on input arguments.
    for year in year_range:
        for dataset in dataset_range:
            print(f"Year: {year} Dataset: {dataset}")            
            outdir            = f'/home/ucfacb0/Scratch/Emission_Inv_Final/{grid_res}/{levels}/{year}/'    
            rocket_data       = RocketData(launch_path, reentry_path, rocket_info_path, pei_path)
            print("Successfully loaded input databases.")
            # Go through process of gridding and saving rocket emissions:
            emis_data = OutputEmis(rocket_data, 
                                grid_res, 
                                timestep, 
                                levels, 
                                booster_cutoff_alts, 
                                MECO_alts, 
                                SEI_alts, 
                                SECO_alts,
                                reentry_levs, 
                                vert_filepath, 
                                months,
                                dataset,
                                year,
                                ground_landing_list,
                                stage_alt_dict,
                                vertical_source,
                                args.save_diagnostics,
                                args.save_emissions,
                                outdir,
                                model_alt
                                )

            check_total_emissions(year,dataset,args.save_diagnostics,grid_res,levels)