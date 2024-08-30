import numpy as np
import matplotlib.pyplot as plt

def plot_luminescence_summary_separate(concentrations, summary_dict):
    """
    Plot separate luminescence data plots for each condition using summary statistics.
    
    Parameters:
    concentrations (array-like): Drug concentrations
    summary_dict (dict): Dictionary where keys are condition names and values are tuples
                         of (mean_array, std_array) for each condition
    """
    concentrations = np.array(concentrations)
    
    for condition, (means, stds) in summary_dict.items():
        plt.figure(figsize=(10, 6))
        
        plt.errorbar(np.log10(concentrations), means, yerr=stds, 
                     fmt='o-', capsize=5, color='b')
        
        plt.xlabel('Log Drug Concentration', fontsize = 20)
        plt.ylabel('Mean Luminescence', fontsize = 20)
        plt.xticks(fontsize= 20)
        plt.yticks(fontsize= 20)
        plt.title(f'Luminescence vs Drug Concentration for {condition}', fontsize = 25)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        
        plt.tight_layout()
        plt.show()


# Example usage:
if __name__ == "__main__":
    concentrations = [1000, 100, 10, 1, 0.1, 0.01, 0.001]
    
    summary_dict = {
        "GL261 1e4": (
            np.array([400.6666667, 3926.6666677, 7097, 4824.666667, 7428.666667, 18988.66667, 291223]),
            np.array([277.4334755, 714.5070562, 2255.533418, 1819.202114, 4814.441435, 5528.620925, 302198.0097]) 
        ),
        "CT2A 1e4": (np.array([4498.333333, 1200.666667, 1907.666667, 976, 945, 6083, 193051.3333]), np.array([2111.521805, 772.9038319, 88.79376855, 234.0256396, 281.2774431, 2568.988712, 54235.4181])),
        "GL261 5e3": (
            np.array([126.3333333, 2892.666667, 12866.66667, 3895.666667, 7654.333333, 280484, 122238.3333]),
            np.array([22.50185178, 404.5297682, 1336.998255, 1118.602849, 1830.166477, 31320.06419, 49842.75434])   
        ),
        "CT2A 5e3": (
            np.array([538, 647.6666667, 2068.666667, 1312.333333, 1449.666667, 4085.666667, 25130.66667]),
            np.array([529.35149, 238.0217077, 894.0080164, 905.7893427, 969.9073839, 1563.189794, 28371.36712])
        )
    }
    
    plot_luminescence_summary_separate(concentrations, summary_dict)