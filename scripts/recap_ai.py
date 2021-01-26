#!/usr/bin/env python3

import msprime, pyslim
import os, glob, sys
import numpy as np

def main():

        # directory = '/exports/eddie/scratch/s1854903/starlike/s_weak'
        # directory_out = '/exports/eddie/scratch/s1854903/starlike/s_weak'

        directory = "/ceph/users/gbisshop/output/adaptive_introgression/trees"
        directory_out = "/scratch/gbisshop/ai_recap"

        print(os.listdir(directory))

        file_list = [file for file in os.listdir(directory) if file.endswith('.trees')]
        #file_list = [file for file in os.listdir(directory) if file.endswith('.recap')]


        
        Ne = 1e4
        recombination_rate = 1e-6
        burnInTime = 1000
        Tdonor = (4)*Ne
        sample_size = 4


        population_configurations = [
                msprime.PopulationConfiguration(initial_size = 1), # empty pop for pyslim
                msprime.PopulationConfiguration(initial_size = Ne), #recipient jpopulation
                msprime.PopulationConfiguration(initial_size = 1), #donor population
         ]
        migration_matrix = [
                [0,0,0],
                [0,0,0],
                [0,0,0],
        ]
        # demographic_events = [
        #     msprime.MassMigration(time=Tdonor,source = 2, destination = 1, proportion = 1.0),
        # ]
        
        # tree_sequence = msprime.simulate(recombination_rate=recombination_rate,
        #         population_configurations=population_configurations,
        #         migration_matrix=migration_matrix,
        #         demographic_events=demographic_events
        #         )
        # quit()
        #create numpy array
        #sites = [250*i for i in range(40)]+[10000+1000*i for i in range(15)]
        #results = np.zeros((len(sites), sample_size, len(file_list)))

        for idx, treefile in enumerate(file_list):
                #print("recapitating")
                ts = pyslim.load(os.path.join(directory, treefile))
                
                Tdonor = Tdonor + ts.slim_generation

                demographic_events = [
                        msprime.MassMigration(time=Tdonor,source = 2, destination = 1, proportion = 1.0),
                ]
                ts = ts.recapitate(    
                        recombination_rate=recombination_rate,
                        population_configurations=population_configurations,
                        migration_matrix = migration_matrix,
                        demographic_events=demographic_events
                        )
                ts.dump(os.path.join(directory_out,'tree'+str(idx)+'.recap'))
                
if __name__ == '__main__':
        main()
