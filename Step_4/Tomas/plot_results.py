import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_structure_data(data_dicts, title="Structure Comparison"):
    """
    Plot each parameter as a separate subplot, showing each structure as a point.
    Parameters
    ----------
    data_dicts : dict
        Format: {
            "Structure1": {"param1": val1, "param2": val2, ...},
            "Structure2": {"param1": val1, "param2": val2, ...},
            ...
        }
    title : str
        Overall title for the figure.
    """
    # Collect all parameter names
    all_params = sorted({p for struct in data_dicts.values() for p in struct.keys()})
    n_params = len(all_params)

    # Make subplots (auto grid layout)
    ncols = 2 if n_params > 1 else 1
    nrows = int(np.ceil(n_params / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows))
    axes = np.array(axes).reshape(-1)

    structures = list(data_dicts.keys())

    for i, param in enumerate(all_params):
        ax = axes[i]
        values = [data_dicts[s].get(param, np.nan) for s in structures]
        ax.plot(structures, values, 'o', label=param)
        ax.set_title(param)
        ax.set_ylabel("Value")
        ax.set_xlabel("Structure")
        ax.grid(True, linestyle=':', alpha=0.6)

    # Hide extra axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    fig.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.show()


# === Your Data ===
E_TiFe = -448.58705834  # Ry
data = {
    "TiFe": {
        "energy (Ry)": E_TiFe,
        "volume (a.u.^3)": 174.09622,
    },
    "TiFe Fx": {
        "energy (Ry)": -449.75327905,
        "volume (a.u.^3)": 191.46797,
        "energy per H (Ry)": -449.75327905 - E_TiFe / 1,
    },
    "TiFe Fxy": {
        "energy (Ry)": -450.92277302,
        "volume (a.u.^3)": 207.87022,
        "energy per H (Ry)": (-450.92277302 - E_TiFe) / 2,
    },
    "TiFe Fxyz": {
        "energy (Ry)": -452.07629745,
        "volume (a.u.^3)": 230.52229,
        "energy per H (Ry)": (-452.07629745 - E_TiFe) / 3,
    },

    "TiFe Ex": {
        "energy (Ry)": -449.71725562,
        "volume (a.u.^3)": 189.02898,
        "energy per H (Ry)": -449.71725562 - E_TiFe,
    },

    "TiFe Exy": {
        "energy (Ry)":-450.85292737 ,
        "volume (a.u.^3)": 221.8286 ,
        "energy per H (Ry)": ( -450.85292737 - E_TiFe)/2,
    },

    "TiFe Exyz": {
        "energy (Ry)": -451.87524894,
        "volume (a.u.^3)": 260.71462,
        "energy per H (Ry)": (-451.87524894 - E_TiFe) /3,
    },
    "TiFe I": {
        "energy (Ry)": -449.67859919,
        "volume (a.u.^3)": 197.59556,
        "energy per H (Ry)": -449.67859919 - E_TiFe,
    },
}

# === Print as Table ===
df = pd.DataFrame(data).T  # Transpose so structures are rows
print("\n=== Structure Data Table ===")
print(df.to_string(float_format=lambda x: f"{x:10.8f}"))

# === Plot Data ===
plot_structure_data(data, title="TiFe Structure Energies and Volumes")
