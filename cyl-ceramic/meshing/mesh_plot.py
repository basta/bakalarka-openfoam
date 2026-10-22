# Import necessary libraries
import matplotlib.pyplot as plt
import numpy as np

# --- Configuration ---
num_comsol_nodes_x = 8  # Number of COMSOL nodes in x
num_comsol_nodes_y = 8  # Number of COMSOL nodes in y
num_openfoam_nodes = 40 # Total number of OpenFOAM nodes
noise_level = 0.05      # How irregular the OpenFOAM mesh is
target_node_index = 15  # Index of the OpenFOAM node to show mapping for
num_source_arrows = 3   # Number of arrows to draw from nearby COMSOL nodes

# --- Generate Mesh Node Coordinates ---

# 1. COMSOL Mesh Nodes (Source) - Regular Grid
x_comsol = np.linspace(0, 1, num_comsol_nodes_x)
y_comsol = np.linspace(0, 1, num_comsol_nodes_y)
X_comsol, Y_comsol = np.meshgrid(x_comsol, y_comsol)
comsol_nodes = np.vstack([X_comsol.ravel(), Y_comsol.ravel()]).T

# 2. OpenFOAM Mesh Nodes (Target) - Irregular/Random
# Ensure reproducibility
np.random.seed(42)
openfoam_nodes = np.random.rand(num_openfoam_nodes, 2)
# Add some slight clustering or structure if desired (optional)
# openfoam_nodes = openfoam_nodes * 0.8 + 0.1 # Confine to a central area

# --- Generate Sample Data on COMSOL Nodes ---
# Example: A simple function like sin(x) * cos(y)
comsol_data = np.sin(comsol_nodes[:, 0] * np.pi) * np.cos(comsol_nodes[:, 1] * np.pi) * 10

# --- Identify Nodes for Mapping Illustration ---

# Select the target OpenFOAM node
target_node = openfoam_nodes[target_node_index]

# Find the closest COMSOL nodes to the target OpenFOAM node
distances = np.sqrt(np.sum((comsol_nodes - target_node)**2, axis=1))
closest_comsol_indices = np.argsort(distances)[:num_source_arrows]
source_nodes = comsol_nodes[closest_comsol_indices]

# --- Create the Plot ---
fig, ax = plt.subplots(figsize=(8, 7))

# 1. Plot COMSOL Nodes (with data as color)
scatter_comsol = ax.scatter(
    comsol_nodes[:, 0],
    comsol_nodes[:, 1],
    c=comsol_data,
    cmap='viridis', # Colormap for the data
    s=80,           # Size of the markers
    marker='o',     # Marker shape
    edgecolor='k',  # Edge color for markers
    linewidth=0.5,
    label='COMSOL Nodes (Source Data)'
)

# 2. Plot OpenFOAM Nodes
ax.scatter(
    openfoam_nodes[:, 0],
    openfoam_nodes[:, 1],
    c='red',
    s=60,
    marker='x',
    label='OpenFOAM Nodes (Target)'
)

# 3. Highlight the specific target OpenFOAM node
ax.scatter(
    target_node[0],
    target_node[1],
    c='red',
    s=150,          # Make it larger
    marker='x',
    linewidth=2,    # Make it thicker
    label='Selected Target Node'
)

# 4. Draw Arrows for Mapping Illustration
for i, source_node in enumerate(source_nodes):
    ax.annotate(
        '', # No text label on the arrow itself
        xy=target_node,           # Point arrow tip points to
        xytext=source_node,       # Point arrow base starts from
        arrowprops=dict(
            arrowstyle="->",
            color='gray',
            linewidth=1.5,
            shrinkA=10, # Shorten arrow start to avoid marker overlap
            shrinkB=10, # Shorten arrow end to avoid marker overlap
            connectionstyle="arc3,rad=0.1" # Slight curve
        )
    )
    # Optionally highlight the source nodes used for mapping
    ax.scatter(
        source_node[0],
        source_node[1],
        s=120,
        facecolors='none', # Make them hollow circles
        edgecolors='blue',
        linewidth=1.5,
        label='Contributing Source Nodes' if i == 0 else "" # Label only once
    )


# --- Add Labels, Title, Legend, Colorbar ---
ax.set_title('Concept: Mapping Field Data from COMSOL to OpenFOAM Mesh')
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')

# Create a combined legend
handles, labels = ax.get_legend_handles_labels()
# Remove duplicate labels if any (like 'Contributing Source Nodes')
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc='upper right')


# Add a colorbar for the COMSOL data
cbar = fig.colorbar(scatter_comsol, ax=ax)
cbar.set_label('Field Value (e.g., E or B component)')

# --- Final Touches ---
ax.set_aspect('equal', adjustable='box') # Ensure aspect ratio is equal
ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, 1.1)
plt.tight_layout() # Adjust layout to prevent labels overlapping

# --- Show/Save Plot ---
# plt.savefig('mesh_mapping_concept.png', dpi=300) # Uncomment to save
plt.show()
