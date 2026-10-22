using CairoMakie
using LaTeXStrings
using Printf # For formatting numbers in strings

# Activate CairoMakie backend for consistent output
CairoMakie.activate!()

# --- Configuration for the plot ---
fontsize_theme = Theme(fontsize=18, Axis=(titlesize=20, xlabelsize=18, ylabelsize=18))
set_theme!(fontsize_theme)

# --- Generate Data: A simple sine wave ---
dt = 0.01 # Smaller dt for smoother curve
signal_period = 2 * π
t_max = signal_period * 2.5 # Show a bit more than two periods
t = 0:dt:t_max
signal = sin.(t)

# --- Define Embedding Parameters ---
embedding_dimension = 5
# Choose a single, well-chosen tau (e.g., quarter period for sine wave)
tau_embedding = round(signal_period / 4, digits=2) # π/2 ≈ 1.57

# --- Choose a Reference Time Point ---
# Let's pick a point on the rising slope for better visual separation of points
t_ref = signal_period / 6 #  π/3, so s(t_ref) is sin(π/3)
s_ref = sin(t_ref)

# --- Create the Figure ---
fig = Figure(size=(1200, 750)) # Adjusted size for clarity and legend

ax = Axis(fig[1, 1],
    title="Visualizing a $(embedding_dimension)-Dimensional Delay Embedding Vector (\$s(t), s(t+\\tau), \\dots, s(t+(m-1)\\tau)\$)",
    xlabel=L"Time ($t$)",
    ylabel=L"Signal ($s(t) = \sin(t)$)")

# Plot the main sine wave and capture the plot object
signal_line_plot = lines!(ax, t, signal, color=:dodgerblue, linewidth=2.5) # Removed label here, will add to legend_elements

# Adjust y-limits to make space for delay arrows below the x-axis
arrow_base_y = -1.2
arrow_height = 0.1 # Height for the "τ" segment markers
min_plot_y = arrow_base_y - arrow_height - 0.1 # Add a little buffer
ylims!(ax, min_plot_y, 1.3)

# --- Initialize legend elements and labels ---
# Initialize with type Any to allow different plot types
legend_elements = Any[signal_line_plot]
legend_labels = [L"s(t)"]

# --- Colors for the embedding points P₀ to P₄ ---
point_colors = [:black, :crimson, :forestgreen, :darkorange, :darkviolet]

# --- Visualize the embedding points and tau intervals ---
# Plot P₀
p0_scatter = scatter!(ax, [t_ref], [s_ref], color=point_colors[1], marker=:diamond, markersize=15)
p0_text_str = "P_0"
text!(ax, t_ref + 0.05, s_ref + (s_ref > 0.8 ? -0.1 : 0.08), text=LaTeXString(p0_text_str), color=point_colors[1], fontsize=16, align=(:left, :bottom))
vlines!(ax, [t_ref], ymin=arrow_base_y - arrow_height, ymax=s_ref, color=point_colors[1], linestyle=:dash, linewidth=1.5)

push!(legend_elements, p0_scatter)
t_ref_formatted_str = @sprintf("%.2f", t_ref)
p0_legend_str = "P_0 = s(t_{ref}) \\text{ (Ref. at } t_{ref}=$(t_ref_formatted_str)\\text{)}"
push!(legend_labels, LaTeXString(p0_legend_str))


for k in 1:(embedding_dimension-1)
    current_t = t_ref + k * tau_embedding
    current_s = sin(current_t)
    prev_t = t_ref + (k - 1) * tau_embedding

    # Mark the point Pₖ
    pk_scatter = scatter!(ax, [current_t], [current_s], color=point_colors[k+1], marker=:circle, markersize=12)

    pk_label_str = "P_{$k}" # Construct string for P_k
    text!(ax, current_t + 0.05, current_s + (current_s > 0.8 ? -0.1 : 0.08),
        text=LaTeXString(pk_label_str), color=point_colors[k+1], fontsize=16, align=(:left, :bottom))

    vlines!(ax, [current_t], ymin=arrow_base_y - arrow_height, ymax=current_s, color=point_colors[k+1], linestyle=:dash, linewidth=1.5)
    lines!(ax, [prev_t, current_t], [arrow_base_y, arrow_base_y], color=point_colors[k+1], linewidth=3)

    tau_symbol_ls = L"\tau" # Simple L"" string for single symbol
    text!(ax, prev_t + tau_embedding / 2, arrow_base_y + 0.03,
        text=tau_symbol_ls, color=point_colors[k+1], align=(:center, :bottom), fontsize=16)

    # Add to legend
    push!(legend_elements, pk_scatter)
    if k == 1
        legend_str_val = "P_{$k} = s(t_{ref} + \\tau), \\quad \\tau = $(tau_embedding)"
        push!(legend_labels, LaTeXString(legend_str_val))
    else
        legend_str_val = "P_{$k} = s(t_{ref} + $(k)\\tau)"
        push!(legend_labels, LaTeXString(legend_str_val))
    end
end

tau_label_val_str = @sprintf("%.2f", tau_embedding)
general_tau_text_str = "\\text{Time delay } \\tau = $(tau_label_val_str)"
text!(ax, t[1] + 0.5, min_plot_y + 0.05, text=LaTeXString(general_tau_text_str), fontsize=16, align=(:left, :bottom))

axislegend(ax, legend_elements, legend_labels,
    position=:cb, patchsize=(30, 30), titlesize=16, labelsize=14, # Changed :outerright to :cb
    nbanks=embedding_dimension + 1)

save_dir = "figure"
mkpath(save_dir)

save(joinpath(save_dir, "5d_embedding_visualization.pdf"), fig)
save(joinpath(save_dir, "5d_embedding_visualization.png"), fig, px_per_unit=2)

println("Makie script finished. Figure saved in '$(save_dir)' directory as '5d_embedding_visualization.{pdf,png}'.")

fig
