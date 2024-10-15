import sys
import scipy.io
import matplotlib.pyplot as plt

def main(data_file):
    # Load data from the .mat file
    mat_data = scipy.io.loadmat(data_file)
    x = mat_data['x'][0]
    y = mat_data['y'][0]
    styles = mat_data['styles'][0]
    
    # Create the plot
    plt.figure()
    for i in range(len(x)):
        plt.plot(x[i].flatten(), y[i], styles[i][0])
    
    # Check and set the title if it exists
    if 'titleStr' in mat_data and mat_data['titleStr'].size != 0:
        plt.title(mat_data['titleStr'][0])
    
    # Check and set the X label if it exists
    if 'xlabelStr' in mat_data and mat_data['xlabelStr'].size != 0:
        plt.xlabel(mat_data['xlabelStr'][0])
    
    # Check and set the Y label if it exists
    if 'ylabelStr' in mat_data and mat_data['ylabelStr'].size != 0:
        plt.ylabel(mat_data['ylabelStr'][0])
    
    # Check and set the legend if it exists
    if 'legendLabels' in mat_data and mat_data['legendLabels'].size != 0:
        legend_labels = [label[0] for label in mat_data['legendLabels'][0]]
        plt.legend(legend_labels)
    
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 plot_script.py <data_file>")
        sys.exit(1)
    
    data_file = sys.argv[1]
    main(data_file)
