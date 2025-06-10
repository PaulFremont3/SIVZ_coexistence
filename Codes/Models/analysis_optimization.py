import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math as mt
import os
from matplotlib.backends.backend_pdf import PdfPages
from generic_functions import phivs_encounter, load_vector, target_concentrations
import sys
from scipy.stats import pearsonr


def load_vector_str(name, sep):
    with open(name, 'r') as f:
        line =f.readline()
        line =line.rstrip()
        l = [num for num in line.split(sep)]
    f.close()
    return(l)

def load_matrix_str(name, sep):
    with open(name, 'r') as f:
        l = [[str(num) for num in line.rstrip().split(sep)] for line in f]
    f.close()
    return(l)

if __name__ == '__main__':
    POM=str(sys.argv[1]) # 0
    indices=[0,1,2,3]
    names_alga=['Diatom', 'Eukaryote', 'Synechococcus', 'Prochlorochoccus']
    r_virus=[20,80,35,35]
    Vols=load_vector('../trait_data/Vs_5.txt', sep=' ')

    ind1=np.argmin(np.array(Vols[0:100]))
    ind2=np.argmin(np.array(Vols[100:200]))
    ind3=np.argmin(np.array(Vols[200:300]))
    ind4=np.argmin(np.array(Vols[300:400]))

    V1=Vols[0:100][ind1]
    V2=Vols[100:200][ind2]
    V3=Vols[200:300][ind3]
    V4=Vols[300:400][ind4]
    Vs=[V1, V2, V3, V4]

    type_rs=['intracellular']
    to_loop=[(ind, type_r) for ind in indices for type_r in type_rs]
    # lopp through all 4 phytoplankton types
    for k in to_loop:
        print(k)
        indice=k[0]
        type_r=k[1]
        type_a=names_alga[indice]
        if POM=='0':
            filename='SIVZ_'+type_r+'_res_optimization_'+type_a+'.txt'
        else:
            filename='SIVZ_'+type_r+'_res_optimization_'+type_a+'_POM.txt'
        if not os.path.exists(filename):
            continue
        type_a=names_alga[indice]
        # load result of optimization procedure
        data = load_matrix_str(filename, ' ')
        print(len(data))
        # load header of optimization procedure
        if POM=='0':
            columns = load_vector_str('SIVZ_'+type_r+'_res_optimization_header.txt', ' ')
        else:
            columns = load_vector_str('SIVZ_'+type_r+'_res_optimization_header_POM.txt', ' ')

        data=np.array(data)
        nrows=data.shape[0]
        venvs=data[:,0]

        virus_radius=r_virus[indice]
        Te=20 # 20 C
        Host_volume=Vs[indice]

        phi_th, phi_th_a=phivs_encounter(virus_radius, Host_volume, Te,indice)
        if indice==1:
            print(virus_radius)
            print(Host_volume)
            print(Te)
            print(indice)
            print(phi_th)

        data=data.astype('float')
        print(columns)
        df = pd.DataFrame(data, columns=columns)
        df['valid_envs']=venvs


        unique_venvs = df['valid_envs'].unique()

        to_opt='ER_av'
        # top 200 minimum average absolute error
        df=df.nsmallest(200, to_opt)
        df.loc[:,'phi/phir']=phi_th/df['phi']
        #Find the row with the minimum value of average absolute error
        min_row = df[df[to_opt] == df[to_opt].min()]

        # names of variable to extract
        if type_r=='intracellular':
            vars_name=['phi','phi/phir','eps' ,'dv', 'dv2', 'dz2', 'cost', 'ED_u', 'ED_m', 'ED_o', 'ED_av', 'ER_u', 'ER_m', 'ER_o', 'ER_av']
            vars_name_bis=['phi', 'eps' ,'dv', 'dv2', 'dz2', 'cost']
            if POM=='1':
                vars_name=['phi','phi/phir','eps' ,'dv', 'dv2', 'dz2', 'cost', 'Pc','ED_u', 'ED_m', 'ED_o', 'ED_av', 'ER_u', 'ER_m', 'ER_o', 'ER_av']
                vars_name_bis=['phi', 'eps' , 'dv2', 'dz2', 'cost', 'Pc']
        else:
            vars_name=['phi','phi/phir','dv', 'dv2', 'dz2', 'cost', 'ED_u', 'ED_m', 'ED_o', 'ED_av', 'ER_u', 'ER_m', 'ER_o', 'ER_av']
            vars_name_bis=['phi','dv', 'dv2', 'dz2', 'cost']
            if POM=='1':
                vars_name=['phi','phi/phir','dv', 'dv2', 'dz2', 'cost', 'Pc','ED_u', 'ED_m', 'ED_o', 'ED_av', 'ER_u', 'ER_m', 'ER_o', 'ER_av']
                vars_name_bis=['phi', 'dv2', 'dz2', 'cost', 'Pc']
        data_ranges=[]
        data_ranges_cols=[]
        for v in vars_name:
            mi=min(df[v])
            mx=max(df[v])
            data_ranges.append(mi)
            data_ranges.append(mx)
            data_ranges_cols.append('min '+v)
            data_ranges_cols.append('max '+v)


        data_ranges_tow=pd.DataFrame([data_ranges], columns=data_ranges_cols)

        if type_r=='intracellular':
            if POM=='0':
                best_match=min_row.iloc[:, np.r_[0:7, 8:16,min_row.shape[1]-1]]
            else:
                best_match=min_row.iloc[:, np.r_[0:8, 9:17,min_row.shape[1]-1]]
        else:
            if POM=='0':
                best_match=min_row.iloc[:, np.r_[0:6, 7:15, min_row.shape[1]-1]]
            else:
                best_match=min_row.iloc[:, np.r_[0:7, 8:16, min_row.shape[1]-1]]
        
        # save best parameters and ranges of the 200 best performing parameter combination
        if POM=='0':
            best_match.to_csv('SIVZ_optimum_fit_'+type_a+'_'+type_r+'.txt', sep=' ', index=False)
            data_ranges_tow.to_csv('SIVZ_range_fit_'+type_a+'_'+type_r+'.txt', sep=' ', index=False)
        else:
            best_match.to_csv('SIVZ_optimum_fit_'+type_a+'_'+type_r+'_POM.txt', sep=' ', index=False)
            data_ranges_tow.to_csv('SIVZ_range_fit_'+type_a+'_'+type_r+'_POM.txt', sep=' ', index=False)

        # errors and errors ranges to target concentrations for the best and the 200 best paramter combination
        target_conc_u, target_conc_m, target_conc_o=target_concentrations(indice)
        target_concs=[target_conc_u, target_conc_m, target_conc_o]
        trs=['P', 'V', 'Z', 'Inf', 'PK']
        envs=['_u', '_m', '_o']
        col_names=[] 
        err_to_tracers=[]
        for i in range(3):
            for j in range(5):
                ev=envs[i]
                tr=trs[j]
                if j<3:
                    err_to_tracers.append(min_row[tr+ev].iloc[0])
                    col_names.append(tr+ev)
                    err=(min_row[tr+ev].iloc[0]-target_concs[i][j])*100/target_concs[i][j]
                    err_to_tracers.append(err)
                    col_names.append('ER_'+tr+ev)
                else:
                    err_to_tracers.append(min_row[tr+ev].iloc[0])
                    col_names.append(tr+ev)
        err_to_tracers=pd.DataFrame([err_to_tracers], columns=col_names)
        if POM=='0':
            err_to_tracers.to_csv('SIVZ_optimum_fit_errors_tracers_'+type_a+'_'+type_r+'.txt', sep=' ', index=False)
        else:
            err_to_tracers.to_csv('SIVZ_optimum_fit_errors_tracers_'+type_a+'_'+type_r+'_POM.txt', sep=' ', index=False)

        # pdf of histograms of parameters distribution
        # Define the number of bins (20 breaks)
        num_bins = 20
        if POM=='0':
            pp = PdfPages('Histogram_optimization_'+type_a+'_'+type_r+'.pdf')
        else:
            pp = PdfPages('Histogram_optimization_'+type_a+'_'+type_r+'_POM.pdf')
        for v in vars_name:
            data=df[v]
            if len(data)>0:
                if v=='eps':
                    data=1/df[v]
                data_min, data_max = min(data), max(data)

                # Create bin edges based on the min and max of the data
                bins = np.linspace(data_min, data_max, num_bins + 1)

                # Extended x-axis range
                if v=='phi':
                    plot_xmin = 1e-12
                    plot_xmax = 1e-8
                elif v=='eps':
                    plot_xmin = 1
                    plot_xmax = 10000
                elif v=='dv':
                    plot_xmin = 0.01
                    plot_xmax = 1
                elif v=='dv2':
                    plot_xmin = 10
                    plot_xmax = 2000
                elif v=='dz2':
                    plot_xmin = 10
                    plot_xmax = 30
                elif v=='cost':
                    plot_xmin =0.5
                    plot_xmax = 1
                elif v=='phi/phir':
                    plot_xmin=1
                    plot_xmax=phi_th/1e-11
                else:
                    plot_xmin = 0
                    plot_xmax = data_max
                    if v=='ER_av':
                        plot_xmax=100

                fig, ax = plt.subplots(figsize=(4,4))
                # Plot the histogram
                sns.kdeplot(data, fill=True, color='grey', alpha=0.3)
            
                # Add a vertical line at a specific point (e.g., 0.5)
                if best_match.shape[0]==0:
                    vertical_line_x = best_match[v]
                    if v=='eps':
                        vertical_line_x = 1/best_match[v]
                    plt.axvline(vertical_line_x, color='red', linestyle='-', label=f'Line at {vertical_line_x}')
                else:
                    vertical_line_x =max(best_match[v])
                    vertical_line_x_bis =min(best_match[v])
                    if v=='eps':
                        vertical_line_x =max(1/best_match[v])
                        vertical_line_x_bis =min(1/best_match[v])
                    if vertical_line_x!=vertical_line_x_bis:
                        plt.axvline(vertical_line_x, color='red', linestyle='-', label=f'Line at {vertical_line_x}')
                        plt.axvline(vertical_line_x_bis, color='red', linestyle='-', label=f'Line at {vertical_line_x_bis}')
                        plt.axvspan(vertical_line_x_bis, vertical_line_x, color='red', alpha=0.3)
                    else:
                        plt.axvline(vertical_line_x, color='red', linestyle='-', label=f'Line at {vertical_line_x}')
                # Customize the plot
                plt.xlabel('Value')
                if v in ['phi', 'eps', 'dv2', 'dv', 'phi/phir']:
                    plt.xscale('log')
                plt.ylabel('Frequency')
                plt.title(v)
                plt.legend()

                # Set the x-axis limits
                plt.xlim(plot_xmin, plot_xmax)

               #  Show the plot
                pp.savefig()
        pp.close()

        # pdf correlograms between prameters for the 200 best paramters combination
        if POM=='0':
            pp = PdfPages('parameters_correlogram_optimization_'+type_a+'_'+type_r+'.pdf')
        else:
            pp = PdfPages('parameters_correlogram_'+type_a+'_'+type_r+'_POM.pdf')

        d3=df['ER_av']
        for i, p1 in enumerate(vars_name_bis):
            for j, p2 in enumerate(vars_name_bis):
                if i<=j:
                    d1=df[p1]
                    d2=df[p2]

                    r_value, p_value = pearsonr(d1, d2)
                    
                    fig, ax = plt.subplots(figsize=(4,4))
                
                    df0 = pd.DataFrame({'d1': d1, 'd2': d2, 'd3': d3})
                    df0_mean = df0.groupby(['d1', 'd2'], as_index=False).mean()

                    # Extract unique values
                    d1_unique = df0_mean['d1']
                    d2_unique = df0_mean['d2']
                    d3_mean = df0_mean['d3']


                    # Create figure and scatter plot
                    sc = plt.scatter(d1_unique, d2_unique, c=d3_mean, cmap="viridis")
                    
                    
                    #sc=plt.scatter(d1, d2, c=d3, cmap="viridis", alpha=0.3)

                    cbar = plt.colorbar(sc)
                    cbar.set_label("Average total error", fontsize=12)

                    # Add R and p-value as text
                    text = f"R = {r_value:.2f}\nP-value = {p_value:.3g}"
                    plt.title(text)

                    # Labels and title
                    plt.xlabel(p1)
                    plt.ylabel(p2)

                    # Show plot
                    pp.savefig()
        pp.close()
