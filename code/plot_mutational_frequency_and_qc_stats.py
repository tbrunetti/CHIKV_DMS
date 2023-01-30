import pandas
import argparse
import typing
import numpy
import os
from matplotlib.ticker import ScalarFormatter

# since ScalarFormatterClss does not allow flexibility in formatting
# of scales, override it by changing the class and the _set_format method
class ScalarFormatterClass(ScalarFormatter):
   def _set_format(self):
      self.format = "%1.1f" # returns format as a value with one digit before the decimal and one significant figure

def filter(df: pandas.DataFrame, minQ: float, minAlt: int, codonRange: str) -> pandas.DataFrame:
    '''
    Potential columns to make filtering decisions
    FWD_MEAN_MIN_QUAL
    REV_MEAN_MIN_QUAL
    FWD_STDDEV_MIN_QUAL
    REV_STDDEV_MIN_QUAL
    FWD_CNT
    REV_CNT
    FWD_DENOM
    REV_DENOM
    CNT
    DENOM
    codonRange = '9-397'
    minAlt = 100
    minQ = 24.0

    '''
    
    df['FWD_REV_MEAN_QUAL'] = df[['FWD_MEAN_MIN_QUAL', 'REV_MEAN_MIN_QUAL']].mean(axis = 1)

    # minQ filtering
    df = df.loc[df['FWD_REV_MEAN_QUAL'] >= minQ]

    # min count filtering
    df = df.loc[df['CNT'] >= minAlt]
    
    # remove perfect ref to alt codons and keep mutations only
    df = df.loc[df['REF_CODON'] != df['CODON']]
    
    # select range (inclusive)
    start, end = [int(region.strip()) for region in codonRange.split('-')]
    df = df.loc[(df['POSITION'] >= start) & (df['POSITION'] <= end)]
    df.reset_index(inplace = True, drop=True)
    
    return df


def get_per_codon_ntNum_mutational_freq(df_hash : dict, outdir : str) -> None:
    import matplotlib.pyplot as plt

    
    for key, value in df_hash.items():
        mut_freqs = dict()
        # summarize total number of nucleotide mutations per row/codon
        value['total_nt_mutations'] = 0
        for idx, row in value.iterrows():
            value.loc[idx,'total_nt_mutations'] = len([nt for nt in range(len(row['REF_CODON'])) if row['REF_CODON'][nt] != row['CODON'][nt]])
        
        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()
        
            # iterate through all nt changes (1, 2, 3)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for ntChanges in range(1,4):
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['total_nt_mutations'] == ntChanges)]['CNT'])
                freq_values['nt_change_{}_counts'.format(ntChanges)] = total_mutational_counts
                freq_values['nt_change_{}_freq'.format(ntChanges)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        
        # sum all nucleotide mutations per codon position
        mut_freqs_df['all_nt_muts_counts']  = mut_freqs_df['nt_change_1_counts'] + mut_freqs_df['nt_change_2_counts'] + mut_freqs_df['nt_change_3_counts']
        mut_freqs_df['all_nt_muts_freq'] = mut_freqs_df['all_nt_muts_counts']/mut_freqs_df['total_read_depth']
        
        combined_fig, axs = plt.subplots(2, 2)
        combined_fig.suptitle('Frequency of nucleotide changes per codon in {} sample'.format(key), fontweight = 'bold')
        
        # single nucleotide changes frequency per codon
        axs[0, 0].plot('CODON_POSITION', 'nt_change_1_freq', data = mut_freqs_df, color=colors[key], alpha=1,  linestyle='-' , linewidth=1)
        axs[0, 0].set_title('single nucleotide', fontsize = 10)
        axs[0, 0].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[0,0].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[0,0].yaxis.set_major_formatter(yScalarFormatter)
        axs[0,0].tick_params(axis='y', labelsize= 9)
        #ax1.title('Frequency of single nucleotide changes per codon', fontweight = 'bold')
        #ax1.xlabel('Codon position', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        #ax1.ylabel('Freqency of mutation', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        
        # double nucleotide changes frequency per codon
        axs[0, 1].plot('CODON_POSITION', 'nt_change_2_freq', data = mut_freqs_df, color=colors[key], alpha=1, linestyle = '-')
        axs[0, 1].set_title('double nucleotide', fontsize = 10)
        axs[0, 1].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[0,1 ].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[0,1].yaxis.set_major_formatter(yScalarFormatter)
        axs[0,1].tick_params(axis='y', labelsize= 9)
        #ax1.title('Frequency of single nucleotide cha
        #ax2.title('Frequency of double nucleotide changes per codon', fontweight = 'bold')
        #ax2.xlabel('Codon position', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        #ax2.ylabel('Freqency of mutation', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        
        # triple nucleotide changes frequency per codon
        axs[1, 0].plot('CODON_POSITION', 'nt_change_3_freq', data = mut_freqs_df, color=colors[key], alpha=1, linestyle='-')
        axs[1, 0].set_title('triple nucleotide', fontsize = 10)
        axs[1, 0].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[1,0].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[1,0].yaxis.set_major_formatter(yScalarFormatter)
        axs[1,0].tick_params(axis='y', labelsize= 9)
        #ax1.title('Frequency of single nucleotide cha
        #ax3.title('Frequency of triple nucleotide changes per codon', fontweight = 'bold')
        #ax3.xlabel('Codon position', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        #ax3.ylabel('Freqency of mutation', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')

        # triple nucleotide changes frequency per codon
        axs[1, 1].plot('CODON_POSITION', 'all_nt_muts_freq', data = mut_freqs_df, color=colors[key], alpha=1, linestyle='-')
        axs[1, 1].set_title('sum of all nucleotides', fontsize = 10)
        axs[1, 1].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[1, 1].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[1,1].yaxis.set_major_formatter(yScalarFormatter)
        axs[1,1].tick_params(axis='y', labelsize= 9)
        #ax4.title('Frequency of all nucleotide changes per codon', fontweight = 'bold')
        #ax4.xlabel('Codon position', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        #ax4.ylabel('Freqency of mutation', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')

        combined_fig.tight_layout()
        combined_fig.savefig(os.path.join(outdir, 'freq_of_nt_changes_per_codon_{}.png'.format(key)), dpi=600, format = 'png')
        

def get_per_codon_aaTypeChange_mutational_freq(df_hash: dict, colors: dict, outdir : str) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    from matplotlib.lines import Line2D
  
    combined_samples = dict()

    for key, value in df_hash.items():
        mut_freqs = dict()
        # summarize total number of nucleotide mutations per row/codon
        value['aaType'] = ''
        
        for idx, row in value.iterrows():
            if row['REF_AA'] == row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'synonymous'
            elif row['REF_AA'] != row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'nonsynonymous'
            elif row['AA'] == '.':
                value.loc[idx,'aaType'] = 'stop'
            else:
                print(row)

        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()

            # iterate through all aa changes (synonynous, nonsynonymous, stop)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for aaTypes in ['synonymous', 'nonsynonymous', 'stop']:
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['aaType'] == aaTypes)]['CNT'])
                freq_values['aa_type_change_{}_counts'.format(aaTypes)] = total_mutational_counts
                freq_values['aa_type_change_{}_freq'.format(aaTypes)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        combined_samples[key] = mut_freqs_df
        
        combined_fig, axs = plt.subplots(2, 2)
        combined_fig.suptitle('Frequency of mutational changes per codon in {} DNA'.format(key), fontweight = 'bold')
        
        # plot synonymous change codon frequency
        axs[0, 0].plot('CODON_POSITION', 'aa_type_change_synonymous_freq', data = mut_freqs_df, color='blue', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[0, 0].set_title('synonymous', fontsize = 12)
        axs[0, 0].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[0,0].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[0,0].yaxis.set_major_formatter(yScalarFormatter)
        axs[0,0].tick_params(axis='y', labelsize= 9)
        #axs[0, 0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))

        
        # plot nonsynonymous codon change freqeuncy
        axs[0, 1].plot('CODON_POSITION', 'aa_type_change_nonsynonymous_freq', data = mut_freqs_df, color='green', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[0, 1].set_title('nonsynonymous', fontsize = 12)
        axs[0, 1].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[0, 1].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[0,1].yaxis.set_major_formatter(yScalarFormatter)
        axs[0,1].tick_params(axis='y', labelsize= 9)
        
        # plot stop codon frequency
        axs[1, 0].plot('CODON_POSITION', 'aa_type_change_stop_freq', data = mut_freqs_df, color='gold', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[1, 0].set_title('stop', fontsize = 12)
        axs[1, 0].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[1, 0].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[1,0].yaxis.set_major_formatter(yScalarFormatter)
        axs[1,0].tick_params(axis='y', labelsize= 9)
        
        # overlay all line plots
        axs[1, 1].plot('CODON_POSITION', 'aa_type_change_stop_freq', data = mut_freqs_df, color='gold', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[1, 1].plot('CODON_POSITION', 'aa_type_change_nonsynonymous_freq', data = mut_freqs_df, color='green', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[1, 1].plot('CODON_POSITION', 'aa_type_change_synonymous_freq', data = mut_freqs_df, color='blue', alpha=0.5,  linestyle='-' , linewidth=1)
        axs[1, 1].set_title('overlay', fontsize = 12)
        axs[1, 1].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[1, 1].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[1,1].yaxis.set_major_formatter(yScalarFormatter)
        axs[1,1].tick_params(axis='y', labelsize= 9)
        
        combined_fig.tight_layout()
        combined_fig.savefig(os.path.join(outdir, 'freq_of_aa_changes_per_codon_lineplot_{}.png'.format(key)), dpi=600, format = 'png')


    combined_sample_fig, axs = plt.subplots(1, 3, figsize=(9, 3))
    combined_sample_fig.suptitle('Frequency of mutational changes per codon'.format(key), fontweight = 'bold')
    legend_labels=[]
    for colorPick, (combinedKey, combinedValue) in enumerate(combined_samples.items()):
        #plt.legend(title='Parameter where:')
        legend_labels.append(combinedKey)
        #colors = ['#CC79A7', 'teal']
        # plot synonymous change codon frequency
        axs[0].plot('CODON_POSITION', 'aa_type_change_synonymous_freq', data = combinedValue, color=colors[combinedKey], alpha=1,  linestyle='-' , linewidth=1)
        axs[0].set_title('synonymous', fontsize = 12)
        axs[0].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[0].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[0].yaxis.set_major_formatter(yScalarFormatter)
        axs[0].tick_params(axis='y', labelsize= 9)
        #axs[0, 0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
            
        # plot nonsynonymous codon change freqeuncy
        axs[1].plot('CODON_POSITION', 'aa_type_change_nonsynonymous_freq', data = combinedValue, color=colors[combinedKey], alpha=1,  linestyle='-' , linewidth=1)
        axs[1].set_title('nonsynonymous', fontsize = 12)
        axs[1].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[1].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[1].yaxis.set_major_formatter(yScalarFormatter)
        axs[1].tick_params(axis='y', labelsize= 9)

            
        # plot stop codon frequency
        axs[2].plot('CODON_POSITION', 'aa_type_change_stop_freq', data = combinedValue, color=colors[combinedKey], alpha=1,  linestyle='-' , linewidth=1)
        axs[2].set_title('stop', fontsize = 12)
        axs[2].ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs[2].set_ylim([0,0.065])
        yScalarFormatter = ScalarFormatterClass(useMathText=True)
        yScalarFormatter.set_powerlimits((0,0))
        axs[2].yaxis.set_major_formatter(yScalarFormatter)
        axs[2].tick_params(axis='y', labelsize= 9)
    # add manual legend
    axs[2].legend([Line2D([0], [0], color=colors[legend_labels[0]], lw=4),
                Line2D([0], [0], color=colors[legend_labels[1]], lw=4), 
                Line2D([0], [0], color=colors[legend_labels[2]], lw=4)], legend_labels)
    combined_sample_fig.tight_layout()
    combined_sample_fig.savefig(os.path.join(outdir, 'freq_of_aa_changes_per_codon_lineplot_sample_overlay.png'), dpi=600, format = 'png')


def get_per_sample_aaTypeChange_mutational_freq_stackedBarPlot(df_hash: dict, outdir : str) -> None:
    
    import matplotlib.pyplot as plt
    samples = []
    nonsyn = []
    syn = []
    stop = []
    
    # key is sample name, value is cleaned dictionary post filter method
    for key,value in df_hash.items():
        print(key)
    
        mut_freqs = dict()
        # summarize total number of nucleotide mutations per row/codon
        value['aaType'] = ''
        
        for idx, row in value.iterrows():
            if row['REF_AA'] == row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'synonymous'
            elif row['REF_AA'] != row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'nonsynonymous'
            elif row['AA'] == '.':
                value.loc[idx,'aaType'] = 'stop'
            else:
                print(row)

        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()

            # iterate through all aa changes (synonymous, nonsynonymous, stop)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for aaTypes in ['synonymous', 'nonsynonymous', 'stop']:
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['aaType'] == aaTypes)]['CNT'])
                freq_values['aa_type_change_{}_counts'.format(aaTypes)] = total_mutational_counts
                freq_values['aa_type_change_{}_freq'.format(aaTypes)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
    
        # get the average of the frequencies for each mutation across all codons
        samples.append(key)
        nonsyn.append(mut_freqs_df.aa_type_change_nonsynonymous_freq.mean(axis=0))
        syn.append(mut_freqs_df.aa_type_change_synonymous_freq.mean(axis=0))
        stop.append(mut_freqs_df.aa_type_change_stop_freq.mean(axis=0))
        

    nonsyn = numpy.array(nonsyn)
    syn = numpy.array(syn)
    stop = numpy.array(stop)
    
    
    plt.bar(samples, nonsyn, color='#E69F00')
    plt.bar(samples, syn, bottom=nonsyn, color='#009E73')
    plt.bar(samples, stop, bottom=nonsyn+syn, color='#CC79A7')
    plt.legend(['nonsynonymous', 'synonymous', 'stop'])
    plt.ylabel('average frequences')
    plt.xlabel('samples')
    plt.title('mean mutational frequency across all codons')
    plt.ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0), useMathText=True)
    yScalarFormatter = ScalarFormatterClass(useMathText=True)
    yScalarFormatter.set_powerlimits((0,0))
    plt.gca().yaxis.set_major_formatter(yScalarFormatter)
    plt.tick_params(axis='y', labelsize= 9)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'freq_of_aa_changes_across_codons_all_samples_stackedBarPlot.png'), dpi=600, format = 'png')


def get_per_sample_ntNum_mutational_freq_stackedBarPlot(df_hash: dict, nonsynOnly: bool, outdir : str) -> None:
    import matplotlib.pyplot as plt

    one = [] # frequencies of zeroes are included at every codon position
    two = [] # frequencies of zeroes are included at every codon position
    three = [] # frequencies of zeroes are included at every codon position
    samples= []
    
    for key, value in df_hash.items():
        mut_freqs = dict()
                # removes synonymous amino acid changes
        
        if nonsynOnly:
            print('Removing synonymous amino acid changes')
            value = value.loc[value['REF_AA'] != value['AA']]
        
        # summarize total number of nucleotide mutations per row/codon
        value['total_nt_mutations'] = 0
        for idx, row in value.iterrows():
            value.loc[idx,'total_nt_mutations'] = len([nt for nt in range(len(row['REF_CODON'])) if row['REF_CODON'][nt] != row['CODON'][nt]])
        
        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()
        
            # iterate through all nt changes (1, 2, 3)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for ntChanges in range(1,4):
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['total_nt_mutations'] == ntChanges)]['CNT'])
                freq_values['nt_change_{}_counts'.format(ntChanges)] = total_mutational_counts
                freq_values['nt_change_{}_freq'.format(ntChanges)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        
        # sum all nucleotide mutations per codon position
        mut_freqs_df['all_nt_muts_counts']  = mut_freqs_df['nt_change_1_counts'] + mut_freqs_df['nt_change_2_counts'] + mut_freqs_df['nt_change_3_counts']
        mut_freqs_df['all_nt_muts_freq'] = mut_freqs_df['all_nt_muts_counts']/mut_freqs_df['total_read_depth']
        
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
    
        # get the average of the frequencies for each mutation across all codons
        samples.append(key)
        one.append(mut_freqs_df.nt_change_1_freq.mean(axis=0))
        two.append(mut_freqs_df.nt_change_2_freq.mean(axis=0))
        three.append(mut_freqs_df.nt_change_3_freq.mean(axis=0))
                

    one = numpy.array(one)
    two = numpy.array(two)
    three = numpy.array(three)
    
    combined_fig, axs = plt.subplots(1, 1, figsize=(7,7))
    
    axs.bar(samples, one, color='lightsteelblue')
    axs.bar(samples, two, bottom=one, color='cornflowerblue')
    axs.bar(samples, three, bottom=one+two, color='royalblue')
    axs.legend(['1 nt', '2 nt', '3 nt'])
    axs.set_ylabel('average frequences')
    axs.set_xlabel('samples')
    axs.ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0), useMathText=True)
    yScalarFormatter = ScalarFormatterClass(useMathText=True)
    yScalarFormatter.set_powerlimits((0,0))
    axs.yaxis.set_major_formatter(yScalarFormatter)
    axs.tick_params(axis='y', labelsize= 9)
    if nonsynOnly:
        axs.set_title('Mean frequency of nucleotide nonsynonymous mutations \n across all codons')
    else:
        axs.set_title('Mean frequency of nucleotide mutation types \n across all codons')

    combined_fig.tight_layout()
    combined_fig.savefig(os.path.join(outdir, 'freq_of_nt_changes_across_codons_all_samples_stackedBarPlot.png'), dpi=600, format = 'png')

    
def get_per_codon_aa_mutational_information_logoplot(df_hash: dict, nonsynOnly: bool, outdir: str) -> None:
    import logomaker
    import matplotlib.pyplot as plt
    import math
    import matplotlib.ticker as mtick

    all_aa_counts = {}
        
    for key, value in df_hash.items():
        print(key)
        mut_freqs = dict()

        # removes synonymous amino acid changes
        if nonsynOnly:
            print('Removing synonymous amino acid changes')
            value = value.loc[value['REF_AA'] != value['AA']]

        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()

            # iterate through all aa changes
            freq_values = {}
            for aa in set(value['AA']):
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['AA'] == aa)]['CNT'])
                freq_values['{}'.format(aa)] = total_mutational_counts
            
            mut_freqs[codon] = freq_values
            
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        #mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        mut_freqs_df.rename(columns={'.':'X'}, inplace = True) # so elongated stop codon does not look like an I or L
                
        # Counts matrix -> Information matrix
        height_per_row = 2
        width_per_col = 7
        start=9
        end=50
        increment=end-start
        endPos=max(mut_freqs_df.index)
        num_pos=endPos-start+1
        num_rows = math.ceil(num_pos/increment)
        fig = plt.figure(figsize=[width_per_col * 1, height_per_row * num_rows])
        if nonsynOnly:
            fig.suptitle('Nonsynonymous mutations for {} library'.format(key), fontsize=15)
        else:
            fig.suptitle('All mutations for {} library'.format(key), fontsize=15)
        #fig = plt.figure()

        for i in range(0, num_rows):
            # TODO also check when the last row is plotted since we may have to scale width differently
            tmp=mut_freqs_df.loc[start:end]
            info_mat = logomaker.transform_matrix(tmp, 
                                      from_type='counts', 
                                      to_type='information')
            ax = plt.subplot2grid((num_rows, 1), (i, 0))
            logomaker.Logo(info_mat, ax=ax, color_scheme="skylign_protein", show_spines=True)
            #logomaker.Logo(info_mat, ax=ax, color_scheme="charge", show_spines=True)
            #logomaker.Logo(info_mat, ax=ax, color_scheme="dmslogo_funcgroup", show_spines=True)

            ax.set_ylim([0,4]) # make sure all axises are the same; be careful to not truncate to early though; 4 bits can represent 16 values
            ax.set_ylabel('bits')
            ax.set_xlabel('codon position')
            ax.xaxis.set_major_locator(mtick.MaxNLocator(integer=True))
            start=end+1
            end=end+increment

        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'logoplot_of_mutations_information_{}.png'.format(key)), dpi=600, format = 'png')


def get_coverage_per_codon(df_hash : dict, outdir : str) -> None:
    import matplotlib.pyplot as plt
    for key, value in df_hash.items():
        fig = plt.figure()
        depth = dict()
        # iterate through each codon position to generate summary of average depth per codon
        for codon in set(value['POSITION']):
            depth[codon] =  value.loc[value['POSITION'] == codon]['DENOM'].mean()
        
        depth_df = pandas.DataFrame.from_dict(depth, orient  = 'index', columns = ['avg_depth'])
        depth_df['codon_pos'] = depth_df.index
        
        plt.plot('codon_pos', 'avg_depth', data = depth_df, color='black', alpha=0.5,  linestyle='-' , linewidth=1)
        plt.title('Depth per codon in {} sample'.format(key), fontweight = 'bold')
        plt.xlabel('Codon position', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        plt.ylabel('Depth \n(avg reads per codon)', fontweight = 'bold', color = 'darkblue', fontsize = '10', horizontalalignment = 'center')
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, 'depth_of_coverage_{}.png'.format(key)), dpi=600, format = 'png')


def get_coverage_per_base():
    pass


def get_combined_mutational_frequencies_stacked_barplot(df_hash: dict, nonsynOnly: bool, includeStop: bool, outdir : str) -> None:
    import matplotlib.pyplot as plt
    import gc

    one = [] # frequencies of zeroes are included at every codon position
    two = [] # frequencies of zeroes are included at every codon position
    three = [] # frequencies of zeroes are included at every codon position
    samples= []
    
    for key, value in df_hash.items():
        mut_freqs = dict()
                # removes synonymous amino acid changes
        
        if nonsynOnly:
            print('Removing synonymous amino acid changes')
            value = value.loc[value['REF_AA'] != value['AA']]
        
        if includeStop == False:
            print('Removing stop codons')
            value = value.loc[value['AA'] != '.']

        # summarize total number of nucleotide mutations per row/codon
        value['total_nt_mutations'] = 0
        for idx, row in value.iterrows():
            value.loc[idx,'total_nt_mutations'] = len([nt for nt in range(len(row['REF_CODON'])) if row['REF_CODON'][nt] != row['CODON'][nt]])
        
        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()
        
            # iterate through all nt changes (1, 2, 3)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for ntChanges in range(1,4):
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['total_nt_mutations'] == ntChanges)]['CNT'])
                freq_values['nt_change_{}_counts'.format(ntChanges)] = total_mutational_counts
                freq_values['nt_change_{}_freq'.format(ntChanges)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        
        # sum all nucleotide mutations per codon position
        mut_freqs_df['all_nt_muts_counts']  = mut_freqs_df['nt_change_1_counts'] + mut_freqs_df['nt_change_2_counts'] + mut_freqs_df['nt_change_3_counts']
        mut_freqs_df['all_nt_muts_freq'] = mut_freqs_df['all_nt_muts_counts']/mut_freqs_df['total_read_depth']
        
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
    
        # get the average of the frequencies for each mutation across all codons
        samples.append(key)
        one.append(mut_freqs_df.nt_change_1_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
        two.append(mut_freqs_df.nt_change_2_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
        three.append(mut_freqs_df.nt_change_3_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
                
    one = numpy.array(one)
    two = numpy.array(two)
    three = numpy.array(three)
    
    samples = []
    nonsyn = []
    syn = []
    stop = []
    
    # key is sample name, value is cleaned dictionary post filter method
    for key,value in df_hash.items():
        print(key)
    
        mut_freqs = dict()
        # summarize total number of nucleotide mutations per row/codon
        value['aaType'] = ''
        
        for idx, row in value.iterrows():
            if row['REF_AA'] == row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'synonymous'
            elif row['REF_AA'] != row['AA'] and (row['AA'] != '.'):
                value.loc[idx,'aaType'] = 'nonsynonymous'
            elif row['AA'] == '.':
                value.loc[idx,'aaType'] = 'stop'
            else:
                print(row)

        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()

            # iterate through all aa changes (synonymous, nonsynonymous, stop)
            freq_values = {}
            freq_values['total_read_depth'] = avg_codon_depth
            for aaTypes in ['synonymous', 'nonsynonymous', 'stop']:
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['aaType'] == aaTypes)]['CNT'])
                freq_values['aa_type_change_{}_counts'.format(aaTypes)] = total_mutational_counts
                freq_values['aa_type_change_{}_freq'.format(aaTypes)] = total_mutational_counts/avg_codon_depth
        
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
    
        # get the average of the frequencies for each mutation across all codons
        samples.append(key)
        nonsyn.append(mut_freqs_df.aa_type_change_nonsynonymous_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
        syn.append(mut_freqs_df.aa_type_change_synonymous_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
        stop.append(mut_freqs_df.aa_type_change_stop_freq.sum(axis=0)/len(set(df_hash[key].POSITION)))
        

    nonsyn = numpy.array(nonsyn)
    syn = numpy.array(syn)
    stop = numpy.array(stop)
    
    fig, axs = plt.subplots(1)
    fig.suptitle('Mutational frequency per codon', fontweight = 'bold')
    
    x_axis = numpy.arange(len(samples))
    axs.bar(x_axis + 0.20, nonsyn, color='#E69F00', width = 0.2)
    axs.bar(x_axis + 0.20, syn, bottom=nonsyn, color='#009E73', width = 0.2)
    axs.bar(x_axis+ 0.20, stop, bottom=nonsyn+syn, color='#CC79A7', width = 0.2)   
    axs.bar(x_axis+ 0.20*2, one, color='lightsteelblue', width = 0.2)
    axs.bar(x_axis+ 0.20*2, two, bottom=one, color='cornflowerblue', width = 0.2)
    axs.bar(x_axis+ 0.20*2, three, bottom=one+two, color='royalblue', width = 0.2)
    
    axs.set_xticks(x_axis + 0.20*2-(0.2/2))  # Set text labels and properties.
    axs.set_xticklabels(samples, minor=False)
    axs.legend(['nonsynonymous', 'synonymous', 'stop', '1 nt', '2 nt', '3 nt'])
    axs.set_ylabel('per-codon mean frequency of mutation')
    axs.set_xlabel('samples')
    axs.ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
    yScalarFormatter = ScalarFormatterClass(useMathText=True)
    yScalarFormatter.set_powerlimits((0,0))
    axs.yaxis.set_major_formatter(yScalarFormatter)
    axs.tick_params(axis='y', labelsize= 9)

    fig.savefig(os.path.join(outdir, 'combined_mutational_frequencies_across_samples_and_mutations_grouped_stackedBarPlot.png'), dpi=600, format = 'png')

    '''
    x_axis = numpy.arange(len(samples))
    plt.bar(x_axis + 0.20, nonsyn, color='#E69F00', width = 0.2)
    plt.bar(x_axis + 0.20, syn, bottom=nonsyn, color='#009E73', width = 0.2)
    plt.bar(x_axis+ 0.20, stop, bottom=nonsyn+syn, color='#999999"', width = 0.2)   
    plt.bar(x_axis+ 0.20*2, one, color='lightsteelblue', width = 0.2)
    plt.bar(x_axis+ 0.20*2, two, bottom=one, color='cornflowerblue', width = 0.2)
    plt.bar(x_axis+ 0.20*2, three, bottom=one+two, color='royalblue', width = 0.2)
    plt.xticks(x_axis + 0.20*2-(0.2/2), samples, rotation=0)  # Set text labels and properties.
    plt.legend(['nonsynonymous', 'synonymous', 'stop', '1 nt (nonsyn)', '2 nt (nonsyn)', '3 nt (nonsyn)'])
    plt.ylabel('per-codon mean frequency of mutation')
    plt.xlabel('samples')
    plt.title('mutation frequency per codon')
    plt.ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
    yScalarFormatter = ScalarFormatterClass(useMathText=True)
    yScalarFormatter.set_powerlimits((0,0))
    axs.yaxis.set_major_formatter(yScalarFormatter)
    axs.tick_params(axis='y', labelsize= 9)

    '''    


def get_aa_diversity(df_hash : dict, outdir : str, colors : dict) -> None:
    import matplotlib.pyplot as plt
    
    
    combined_fig, axs = plt.subplots(1, 1,figsize=(20,7))
    combined_fig.suptitle('Amino acid diversity per codon position', fontweight = 'bold', fontsize=30)
    legend_order = []  
    for key, value in df_hash.items():
        mut_freqs = dict()
        legend_order.append(key)
        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            mut_freqs[codon] = len(set(value.loc[value['POSITION'] == codon]['AA']))
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame.from_dict(mut_freqs, orient = 'index')
        mut_freqs_df['CODON_POSITION'] = mut_freqs_df.index
        mut_freqs_df.rename({0:'total amino acids'}, axis = 1, inplace=True)

    
        print('{}: max diversity is {}, mean diversity is {}'.format(key, max(mut_freqs.values()), float(sum(mut_freqs.values()))/float(len(mut_freqs.values()))))

        # single nucleotide changes frequency per codon
        axs.plot('CODON_POSITION', 'total amino acids', data = mut_freqs_df, color=colors[key], alpha=1,  linestyle='-' , linewidth=1)
        axs.set_title('{}: max diversity is {}, mean diversity is {}'.format(key, max(mut_freqs.values()), float(sum(mut_freqs.values()))/float(len(mut_freqs.values()))), fontsize = 20)
        #axs.ticklabel_format(axis = 'y', style='scientific', scilimits=(0.0, 0.0))
        axs.set_ylim([0,22])
        axs.tick_params(axis='y', labelsize= 20)
        axs.tick_params(axis='x', labelsize= 20)

        axs.set_ylabel('total number of different \n amino acids represented', fontsize = 20)
        axs.set_xlabel('codon position', fontsize = 20)
        
    axs.legend(legend_order)
    combined_fig.savefig(os.path.join(outdir, 'amino_acid_diversity_at_each_codon_across_samples.png'), dpi=600, format = 'png')
   

def get_per_codon_aa_mutational_freq_logoplot(df_hash: dict, nonsynOnly: bool, includeStop: bool, annot: bool, outdir: str) -> None:
    import logomaker
    import matplotlib.pyplot as plt
    import math
    import matplotlib.ticker as mtick


    all_aa_counts = {}
    
    for key, value in df_hash.items():
        print(key)
        mut_freqs = dict()

        # removes synonymous amino acid changes
        if nonsynOnly:
            print('Removing synonymous amino acid changes')
            value = value.loc[value['REF_AA'] != value['AA']]

        # iterate through each codon position to generate summary df
        for codon in set(value['POSITION']):
            avg_codon_depth = value.loc[value['POSITION'] == codon]['DENOM'].mean()

            # iterate through all aa changes
            freq_values = {}
            for aa in set(value['AA']):
                total_mutational_counts = sum(value.loc[(value['POSITION'] == codon) & (value['AA'] == aa)]['CNT'])
                freq_values['{}'.format(aa)] = total_mutational_counts/avg_codon_depth
            
            mut_freqs[codon] = freq_values
        
        # convert key-value pairs to dataframe
        mut_freqs_df = pandas.DataFrame(mut_freqs).T
        mut_freqs_df.rename(columns={'.':'X'}, inplace = True) # so elongated stop codon does not look like an I or L
           
        # formatting logoplot iteration
        height_per_row = 2
        width_per_col = 7
        start=9
        end=50
        increment=end-start
        endPos=max(mut_freqs_df.index)
        num_pos=endPos-start+1
        num_rows = math.ceil(num_pos/increment)
        fig = plt.figure(figsize=[width_per_col * 1, height_per_row * num_rows])
        if includeStop == False:
            mut_freqs_df.drop(axis = 1, labels = 'X', inplace=True)
        if nonsynOnly:
            fig.suptitle('Nonsynonymous mutations for {} library'.format(key), fontsize=15)
        else:
            fig.suptitle('All mutations for {} library'.format(key), fontsize=15)
        #fig = plt.figure()

        for i in range(0, num_rows):
            # TODO also check when the last row is plotted since we may have to scale width differently
            tmp=mut_freqs_df.loc[start:end]
            #info_mat = logomaker.Logo(tmp)
            ax = plt.subplot2grid((num_rows, 1), (i, 0))
            logomaker.Logo(tmp, ax=ax, color_scheme="skylign_protein", show_spines=True)

            #logomaker.Logo(info_mat, ax=ax, color_scheme="charge", show_spines=True)
            #logomaker.Logo(info_mat, ax=ax, color_scheme="dmslogo_funcgroup", show_spines=True)

            ## TESTING ONLY -- NOT AUTOMATED ##
            if annot:
                E3_start = 9
                E3_end = 72
                E2_start = 73
                E2_end = 397 
                y = -0.005
                plot_coords = range(start, end+1)
                try:
                    e2_coords = range(E2_start, E2_end+1)
                    e2_plot_range_start = min(set(e2_coords).intersection(set(plot_coords)))
                    e2_plot_range_end = max(set(e2_coords).intersection(set(plot_coords)))
                except ValueError:
                    print("There is no E2!")
                    e2_plot_range_start = None
                try:
                    e3_coords = range(E3_start, E3_end+1)
                    e3_plot_range_start = min(set(e3_coords).intersection(set(plot_coords)))
                    e3_plot_range_end = max(set(e3_coords).intersection(set(plot_coords)))
                except ValueError:
                    print("There is no E3!")
                    e3_plot_range_start = None

                if e2_plot_range_start != None:
                    ax.plot([e2_plot_range_start, e2_plot_range_end+1], [y, y], color = "lightblue", linewidth = 10, solid_capstyle="butt")
                    ax.text(((e2_plot_range_start + e2_plot_range_end)/2), y*1.4,'E2',fontsize=12)
                    #if len(set([E2_start]).intersection(set(plot_coords))) == 1:
                    #    ax.plot(E2_start, y*2, '^k', markersize=10)
                    #    ax.text(E2_start, y*3.5, start, fontsize=10)
                    #if len(set([E2_end]).intersection(set(plot_coords))) == 1:
                    #    ax.plot(E2_end, y*2, '^k', markersize=10)
                    #    ax.text(E2_end, y*3.5, start, fontsize=10)
                if e3_plot_range_start != None:
                    ax.plot([e3_plot_range_start, e3_plot_range_end+1], [y, y], color = "lightgreen", linewidth = 10, solid_capstyle="butt")
                    ax.text(((e3_plot_range_start + e3_plot_range_end)/2), y*1.4,'E3',fontsize=12)
                    #if len(set([E3_start]).intersection(set(plot_coords))) == 1:
                    #    ax.plot(E3_start, y*2, '^k', markersize=10)
                    #    ax.text(E3_start, y*3.5, start, fontsize=10)
                    #if len(set([E3_end]).intersection(set(plot_coords))) == 1:
                    #    ax.plot(E3_end, y*2, '^k', markersize=10)
                    #    ax.text(E3_end, y*3.5, start, fontsize=10)
            ## END OF ANNOTATION TESTING ##        
                    
            ax.set_ylim([-0.01,0.04]) # make sure all axises are the same; be careful to not truncate to early though; 4 bits can represent 16 values
            ax.set_yticks([0, 0.04])
            ax.set_ylabel('frequency')
            ax.set_xlabel('codon position')
            ax.xaxis.set_major_locator(mtick.MaxNLocator(integer=True))
            start=end+1
            end=end+increment

        fig.tight_layout()
        if includeStop:
            fig.savefig(os.path.join(outdir, 'logoplot_of_mutations_freq_{}_including_stop_codons.png'.format(key)), dpi=600, format = 'png')
        else:
            fig.savefig(os.path.join(outdir, 'logoplot_of_mutations_freq_{}_excluding_stop_codons.png'.format(key)), dpi=600, format = 'png')
            

if __name__ == '__main__':
    
    #TO DO
    
    parser = argparse.ArgumentParser(description = 'Generates plots of mutational frequencies from deep mutational scanning experiments.')
    
    parser.add_argument('--data', required = True, action = 'extend', nargs = '+', type = str, help = 'list of all codon tables per sample as generated by virVar')
    parser.add_argument('--qual', default = 24.0, type = float, help = 'A float specifying the minimum average quality of reads to keep for filtering of results')
    parser.add_argument('--counts', default = 100, type = int, help = 'minimum average number of counts for a codon variant to keep post filtering')
    parser.add_argument('--pos', default = "0-0", type = str, help= 'integer range of codon positions to use in analysis. Default is to use all availble. Ex: 7-100, would include codon 7 through codon 100, inclusive on both ends based upon your input fasta file. WARNING: makes assumption that codon 1 is the first 3 bp of your input fasta reference file')
    parser.add_argument('--samplename', help='name of samples, comma-separated, in the same order as the list of input data')
    parser.add_argument('--nonSynOnly', action = 'store_true', help='If setting this flag, only nonsynonymous mutations will be considered')
    parser.add_argument('--includeStop', action = 'store_true', help='If setting this flag, include STOP codons')
    parser.add_argument('--annotate', action = 'store_true', help='DO NOT SET THIS -- experimental and needs testing')
    parser.add_argument('--outdir', default = os.getcwd(), help = "Path to output directory to write plots")
    parser.add_argument('--colors', type = str, help = 'comma-separated list of hex values in the same order as --samplename')
    args = parser.parse_args()
    
    df_hash = {} # dict to store each samples data as a panda df
    colors = {} # dict to store each sample's colors
    
    # Below commment block should be replaced with for loop at bottom
    # TO DO: TESTING REQUIRED!!
    '''
    # TO DO: read in csv lists to automate and productionize
    
    all_results = pandas.read_csv('/home/tonya/Downloads/mutDNA_S15_cleaned.codon', sep = '\t')
    all_results = pandas.read_csv('/home/tonya/Downloads/wtDNA_S16_cleaned.codon', sep = '\t')
    df = all_results
    
    # TO DO: for each df of all reasults filter first and then add to dict
    
    
    df_hash = {}
    df_hash['mutant'] = df
    df_hash['wildtype'] = df
    '''
    
    
    # The colorblind friendly palette with grey:
    #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    # The colorblind friendly palette with black:
    #cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    
    if os.path.exists(args.outdir):
        pass
    else:
        os.mkdir(args.outdir)
    

    colorKeys = args.samplename.split(',')
    colorValues = args.colors.split(',')
    colors = dict(zip(colorKeys, colorValues))
    
    
    for index,samples in enumerate(args.data):
        sample_df = pandas.read_csv(samples, sep = '\t')
        sample_df = filter(df = sample_df, minQ = args.qual, minAlt = args.counts, codonRange = args.pos)
        df_hash[args.samplename.split(',')[index].strip()] = sample_df
        print(df_hash)
        
    get_per_codon_ntNum_mutational_freq(df_hash = df_hash, outdir = args.outdir)
    get_per_codon_aaTypeChange_mutational_freq(df_hash = df_hash, colors = colors, outdir = args.outdir)
    get_per_sample_aaTypeChange_mutational_freq_stackedBarPlot(df_hash = df_hash, outdir = args.outdir)
    get_per_sample_ntNum_mutational_freq_stackedBarPlot(df_hash = df_hash, nonsynOnly = args.nonSynOnly, outdir = args.outdir)
    get_per_codon_aa_mutational_information_logoplot(df_hash = df_hash, nonsynOnly = args.nonSynOnly, outdir = args.outdir)
    get_coverage_per_codon(df_hash = df_hash, outdir = args.outdir)
    get_combined_mutational_frequencies_stacked_barplot(df_hash = df_hash, nonsynOnly = args.nonSynOnly, includeStop = args.includeStop, outdir = args.outdir)
    get_aa_diversity(df_hash = df_hash, colors = colors, outdir = args.outdir)
    get_per_codon_aa_mutational_freq_logoplot(df_hash = df_hash, nonsynOnly = args.nonSynOnly, includeStop = args.includeStop, annot = args.annotate, outdir = args.outdir)