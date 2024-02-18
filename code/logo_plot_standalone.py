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


def generate_logo_plot(matrix_input:pandas.DataFrame, sample_name:str, annotations:str, increment:int, min_label_len:int) -> None:
    import json
    import logomaker
    import matplotlib.pyplot as plt
    import math
    import matplotlib.ticker as mtick
    

    #assume sorted by start position
    annot_data = pandas.read_csv(annotations)
    annot_data = annot_data.to_dict("index")

    # index maps for annotation start and end positions
    breakpoints = {}
    for key,value in annot_data.items():
        breakpoints[key] = range(value['start'], value['end'])

    '''
    function to apply annotations to top of logoplots
        func name: def match_annotation()
        input:  position_start specifies the codon start postion for the graph line being contructed;
                position_end specifies the last codon position to be plotted for the graph being constructed;
                min_label_len of consecutive amino acids under an annotation bar; anything smaller (non-inclusive) than
                this value will not have text written in the bar, to help prevent text from overflowing into margins
        output: None; adds annotation texts and bars to backend plot but does not return any objects
    '''
    def match_annotation(position_start:int, position_end:int, min_label_len:int) -> None:
        annot_start = -1
        annot_end = -1

        # find the annotation of the lowest codon and the annotation of the largest position codon
        for key,value in breakpoints.items():
            if position_start in value:
                annot_start = key
            if position_end in value:
                annot_end = key    
     
        # if the start and end annotations are the same, plotting is easy as it is all the same color, no breakpoints
        if annot_start == annot_end:
            # rectangular box display control
            ax.plot([position_start-1, position_end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_start]['color'], linewidth = 10, solid_capstyle="butt")
            if (position_end - position_start) >= min_label_len: # text control
              ax.text(((position_start-1 +position_end+1)/2), 1.1 ,annot_data[annot_start]['region_name'],fontsize=10)     
        
        # when the graph has 2 or more annotation bars present on the same line, we must identify all the breakpoints
        else:
          for subannotations in enumerate(range(annot_start, annot_end + 1)):
              if subannotations[0] == 0: # meaning it is the left most annotatation
                  begin = max(annot_data[annot_start]['start'], position_start)
                  end = annot_data[annot_start]['end']
                  ax.plot([begin-1, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_start]['color'], linewidth = 10, solid_capstyle="butt")
                  if (end - begin) >= min_label_len:
                    ax.text(((begin-1 + end+1)/2), 1.1 ,annot_data[annot_start]['region_name'],fontsize=10)
              elif subannotations[0] == len(range(annot_start, annot_end)): # meaning it is the right most or last annotation
                  begin = annot_data[annot_end]['start']
                  end = min(annot_data[annot_end]['end'], position_end)
                  ax.plot([begin, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[annot_end]['color'], linewidth = 10, solid_capstyle="butt")
                  if (end - begin) >= min_label_len:
                    ax.text(((begin + end+1)/2), 1.1 ,annot_data[annot_end]['region_name'],fontsize=10)
              else: # means there are at least 3 annotations, so any annotation in -between the start and end annotations need to be plotted in full
                  begin = annot_data[subannotations[1]]['start']
                  end = annot_data[subannotations[1]]['end']
                  ax.plot([begin, end+1], [y+0.1, y+0.1], alpha= 0.5, color = annot_data[subannotations[1]]['color'], linewidth = 10, solid_capstyle="butt")
                  if (end - begin) >= min_label_len:
                    ax.text(((begin + end+1)/2), 1.1 ,annot_data[subannotations[1]]['region_name'],fontsize=10)                  
               

    # Counts matrix -> Information matrix
    height_per_row = 2
    width_per_col = 15
    start=args.codonStartPos
    end=start+increment
    y = 1.05
    endPos=max(matrix_input.index)
    num_pos=endPos-start+1
    num_rows = math.ceil(num_pos/increment)

    # sets figure size, title, and x-/y-axis labels
    fig = plt.figure(figsize=[width_per_col * 1, height_per_row * num_rows])
    fig.suptitle('All mutations for {} library \n'.format(sample_name), fontsize=15)
    fig.supylabel('amino acid diversity per codon (absence vs presence)')
    fig.supxlabel('codon position')
    
    # loop if a graph for each row of the logo plot; each loop is one plot
    # rows are dictated by calculation above
    for i in range(0, num_rows):
          tmp=matrix_input.loc[start:end]
          ax = plt.subplot2grid((num_rows, 1), (i, 0))
          logomaker.Logo(tmp, ax=ax, color_scheme="skylign_protein", show_spines=True, stack_order = "fixed") # fixed keeps order of aa based on matrix order; in this case alphabetically
     
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

          # if annotations are available, annotate the graphs
          if annotations != None:
            match_annotation(position_start=start, position_end=end, min_label_len= min_label_len)

          # calculate next set of position ranges
          start=end+1
          if (end+increment) > max(matrix_input.index):
              end = max(matrix_input.index) # last graph  may be smaller than the increment amount so scale accordingly
          else:
               end=end+increment

    fig.tight_layout()
    fig.savefig('{}_logoplot_of_diversity_of_amino_acids_present_per_codon_position.png'.format(sample_name), dpi=600, format = 'png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generates logo plot for predefined input matrix")
    parser.add_argument('--input', type = str, required=True, help = "Path to input csv matrix containing data to plot")
    parser.add_argument('--sampleName', type = str, default = "sample_1", help = "string indicating the name to give to sample")
    parser.add_argument('--annotConfig', type = str, default = None, help = "Path to csv containing annotations.  Example file located in ref folder of github repo")
    parser.add_argument('--codonStartPos', type = int, default = None, help = "The position of which codon position you want to start at (must be present in your csv matrix provided to --input; default is to plot every position in your matrix)")
    parser.add_argument('--codonEndPos', type = int, default = None, help = "The position of which codon position you want to end at (must be present in your csv matrix provided to --input; default is to plot every position in your matrix)")
    parser.add_argument('--aaSpacing', type = int, default = 65, help = 'the number of amino acids to show per line on the logo plot')
    parser.add_argument('--minAnnotLabel', type = int, default = 7, help = 'the minimum length of consecutive amino acids under an annotation bar; anything smaller (non-inclusive) than \
                this value will not have text written in the bar, to help prevent text from overflowing into margins')
    args = parser.parse_args()

    formatted_input_matrix  = format_data(data_input=args.input)
    generate_logo_plot(matrix_input=formatted_input_matrix, sample_name = args.sampleName, 
                       annotations = args.annotConfig, increment = args.aaSpacing, min_label_len = args.minAnnotLabel)
