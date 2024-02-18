import pandas
import typing
import argparse

'''
def format_data()
    input is based on what Megan has provided but it needs a minimum of 3 columns (all others will be ignored)
        - POSITION: an integer of the codon position
        - AA: amino acid to plot, generally this is a single letter amino acid symbol
        - MERGE_FRAC: float [0-1] that shows the height/frequency of the amino acid symbol to plot on the logoplot
    output is a pandas dataframe that is formatted for logomaker in input into def generate_logo_plot()
'''

def format_data(data_input:pandas.DataFrame) -> pandas.DataFrame:
     pandas_df = pandas.read_csv(data_input)
     info_df = pandas_df[['POSITION', 'AA', 'MERGE_FRAC']]
     info_pivot_matrix = info_df.pivot(index = "POSITION", columns = "AA", values = "MERGE_FRAC")
     info_pivot_matrix.fillna(0, inplace=True)

     return info_pivot_matrix


def generate_logo_plot(matrix_input:pandas.DataFrame, sample_name:str) -> None:
    import json
    import logomaker
    import matplotlib.pyplot as plt
    import math
    import matplotlib.ticker as mtick
    
    '''
    annotations:
    E3: 1-64
    N link: 65-79
    A domain: 80-198
    Arch 1: 199-236
    B domain: 237-295
    Arch 2: 296-332
    C domain: 333-405 (edited) 


    '''
    # testing annotation config file
    with open("annotations_config.json", "r") as read_json:
         annot_data = json.load(read_json)
    
    #assume sorted by start position
    annot_data = pandas.read_csv("~/amc_projects/logo_plot_custom_Megan_01042024/annotations_config.csv")
    annot_data = annot_data.to_dict("index")

    # index maps
    breakpoints = {}

    for key,value in annot_data.items():
        breakpoints[key] = range(value['start'], value['end'])

    def match_annotation(position_start:int, position_end:int) -> None:
        annot_start = -1
        annot_end = -1
        for key,value in breakpoints.items():
            if position_start in value:
                annot_start = key
            if position_end in value:
                annot_end = key    
     
        if annot_start == annot_end:
            print("same_{}_{}".format(position_start, position_end))
            ax.plot([position_start-1, position_end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_start]['color'], linewidth = 10, solid_capstyle="butt")
            ax.text(((position_start-1 +position_end+1)/2), 1.1 ,annot_data[annot_start]['region_name'],fontsize=10)         
        else:
          for subannotations in enumerate(range(annot_start, annot_end + 1)):
              if subannotations[0] == 0: # meaning it is the left most annotatation
                  begin = max(annot_data[annot_start]['start'], position_start)
                  end = annot_data[annot_start]['end']
                  print('first_{}_{}'.format(begin, end))
                  ax.plot([begin-1, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_start]['color'], linewidth = 10, solid_capstyle="butt")
                  ax.text(((begin-1 + end+1)/2), 1.1 ,annot_data[annot_start]['region_name'],fontsize=10)
              elif subannotations[0] == len(range(annot_start, annot_end)):
                  print("here is the final loop")
                  begin = annot_data[annot_end]['start']
                  end = min(annot_data[annot_end]['end'], position_end)
                  print('last_{}_{}'.format(begin, end))
                  ax.plot([begin, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_end]['color'], linewidth = 10, solid_capstyle="butt")
                  ax.text(((begin + end+1)/2), 1.1 ,annot_data[annot_end]['region_name'],fontsize=10)
              else:
                  begin = annot_data[subannotations[1]]['start']
                  end = annot_data[subannotations[1]]['end']
                  print('middle_{}_{}'.format(begin, end))
                  ax.plot([begin, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[subannotations[1]]['color'], linewidth = 10, solid_capstyle="butt")
                  ax.text(((begin + end+1)/2), 1.1 ,annot_data[subannotations[1]]['region_name'],fontsize=10)                  
               

    # Counts matrix -> Information matrix
    height_per_row = 2
    width_per_col = 7
    start=9
    end=50
    y = 1.05
    increment=end-start
    endPos=max(matrix_input.index)
    num_pos=endPos-start+1
    num_rows = math.ceil(num_pos/increment)

    fig = plt.figure(figsize=[width_per_col * 1, height_per_row * num_rows])
    #fig.suptitle('All mutations for {} library \n'.format(sample_name), fontsize=15)
    
    # loop if a graph for each row of the logo plot; each loop is one plot
    # rows are dictated by calculation above
    for i in range(0, num_rows):
          print(num_rows)
          # TODO also check when the last row is plotted since we may have to scale width differently
          tmp=matrix_input.loc[start:end]
          ax = plt.subplot2grid((num_rows, 1), (i, 0))
          logomaker.Logo(tmp, ax=ax, color_scheme="skylign_protein", show_spines=True)
     
          ax.set_ylim([0,1.5]) # make sure all axises are the same; be careful to not truncate to early though; 4 bits can represent 16 values
          ax.set_yticks([0, 1])
          ax.set_yticklabels(['0.0','1.0'])
          ax.spines['left'].set_color('white')
          ax.spines['right'].set_color('white')
          ax.spines['top'].set_color('white')
          #ax.spines['bottom'].set_color('white')
          #ax.set_ylabel('aa diversity per codon')
          #ax.set_xlabel('codon position')
          ax.xaxis.set_major_locator(mtick.MaxNLocator(integer=True)) 

          print('loopLoc_{}_start_{}_end_{}'.format(i, start, end))
          match_annotation(position_start=start, position_end=end)

          start=end+1
          if (end+increment) > max(matrix_input.index):
              end = max(matrix_input.index)
          else:
               end=end+increment

    fig.tight_layout()
    fig.savefig('{}_logoplot_of_diversity_of_amino_acids_present_per_codon_position.png'.format(sample_name), dpi=600, format = 'png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates logo plot for predefined input matrix")
    parser.add_argument('--input', type = str, required=True, help = "Path to input csv matrix containing data to plot")
    parser.add_argument('--sampleName', type = str, default = "sample_1", help = "string indicating the name to give to sample")
    args = parser.parse_args()

    formatted_logo_matrix = format_data(data_input=args.input)
    generate_logo_plot(matrix_input=formatted_logo_matrix, sample_name = args.sampleName)
