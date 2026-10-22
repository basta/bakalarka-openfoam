using CairoMakie   # Use CairoMakie for high-quality static output (PDF, PNG, SVG)
using LaTeXStrings # For using LaTeX syntax in plot labels and titles

# Activate CairoMakie backend for consistent output
CairoMakie.activate!()

# --- Configuration for the plot ---
fontsize_theme = Theme(fontsize=16, Axis=(titlesize=18, xlabelsize=16, ylabelsize=16))
set_theme!(fontsize_theme)

# --- Generate Data: A simple sine wave ---
dt = 0.05                           # Time step
t_max = 6 * π                       # Maximum time (ensures a few periods of sin(t))
t = 0:dt:t_max                      # Time vector
signal = sin.(t)                    # The signal s(t)

# --- Define Time Delays (τ) for Embedding ---
k_small = 5
τ_s_val = round(k_small * dt, digits=2) # τ_s value

k_good = round(Int, (π / 2) / dt)
τ_g_val = round(k_good * dt, digits=2) # τ_g value

k_large = round(Int, (2 * π) * 1.1 / dt)
τ_l_val = round(k_large * dt, digits=2) # τ_l value


# --- Create Delay Embeddings ---
function delay_embed(data::AbstractVector, k_delay::Int)
    if k_delay <= 0
        error("Delay k_delay (number of steps) must be a positive integer.")
    end
    if length(data) <= k_delay
        error("Time series is too short for the given delay k_delay. Need length > k_delay.")
    end
    x_current = data[1:end-k_delay]
    x_delayed = data[k_delay+1:end]
    return x_current, x_delayed
end

s_curr_small, s_del_small = delay_embed(signal, k_small)
s_curr_good, s_del_good = delay_embed(signal, k_good)
s_curr_large, s_del_large = delay_embed(signal, k_large)


# --- Create the Figure ---
fig = Figure(size=(950, 850)) # Increased height slightly for annotations

# Row 1: Original Time Series
ax_ts = Axis(fig[1, 1:3],
    title=L"(A) Original Time Series: $s(t) = \sin(t)$",
    xlabel=L"$t$",
    ylabel=L"$s(t)$")
lines!(ax_ts, t, signal, color=:dodgerblue, linewidth=2, label=L"s(t)")


# --- Add Delay Visualizations to Time Series Plot (Plot A) ---

# Adjust y-limits to make space for delay arrows below the x-axis
ylims!(ax_ts, -1.7, 1.3) # Original signal is -1 to 1. Made space below.

# Define a reference time point (e.g., first peak)
t_ref = π / 2
s_ref = sin(t_ref)

# Mark the reference point P₀ = s(t_ref)
scatter!(ax_ts, [t_ref], [s_ref], color=:black, markersize=10)
text!(ax_ts, t_ref + 0.15, s_ref + 0.05, text=L"P_0", color=:black, fontsize=14, align=(:left, :bottom))

# Define y-levels for arrows and text to avoid overlap
arrow_y_levels = Dict(
    τ_s_val => -1.25,
    τ_g_val => -1.4,
    τ_l_val => -1.55
)

# Common vertical line from P₀ down to below the lowest arrow's horizontal plane
vlines!(ax_ts, [t_ref], ymin=minimum(values(arrow_y_levels)) - 0.02, ymax=s_ref, color=:dimgray, linestyle=:dash, linewidth=1.5)

delays_for_visualization = [
    (τ_s_val, :crimson, "s"),
    (τ_g_val, :forestgreen, "g"),
    (τ_l_val, :darkorchid, "l")
]

for (tau_val, color_val, subscript_char) in delays_for_visualization
    t_delayed = t_ref + tau_val

    # Ensure sin can be evaluated even if t_delayed is outside 't' vector range for plotting
    # (though for this example, t_max is large enough)
    s_delayed = sin(t_delayed)
    current_arrow_y = arrow_y_levels[tau_val]

    if t_delayed <= t[end] + dt # Check if the delayed point is reasonably within or near plotted t range
        # Mark the delayed point Pᵢ = s(t_ref + τᵢ)
        scatter!(ax_ts, [t_delayed], [s_delayed], color=color_val, markersize=8)
        text!(ax_ts, t_delayed + 0.15, s_delayed + 0.05, text=L"P_{%$(subscript_char)}", color=color_val, fontsize=14, align=(:left, :bottom))

        # Vertical line from delayed point Pᵢ down to its arrow's horizontal plane
        vlines!(ax_ts, [t_delayed], ymin=current_arrow_y, ymax=s_delayed, color=color_val, linestyle=:dash, linewidth=1.5)

        # Draw horizontal arrow representing τ, starting from (t_ref, current_arrow_y)
        arrows!(ax_ts,
            [t_ref], [current_arrow_y], # Arrow start point (x,y)
            [tau_val], [0],             # Arrow vector (dx, dy) meaning length tau_val horizontally
            color=color_val,
            linewidth=2.5,
            arrowsize=12)

        # Add text label for τ value near the arrow
        text!(ax_ts,
            t_ref + tau_val / 2,          # Centered along the arrow's length
            current_arrow_y - 0.03,       # Slightly below the arrow
            text=L"\tau_{%$(subscript_char)} = %$(tau_val)",
            color=color_val,
            align=(:center, :top),      # Text centered, with its top at the y-coord
            fontsize=14)
    end
end
# Add a legend for the main sine wave if desired (currently only P0, Ps etc are implicitly identified by color/text)
# axislegend(ax_ts, [ax_ts.lines[1]], [L"s(t) = \sin(t)"], position=:rt)
# Or, could add labels to scatter calls and use a full legend, but might be too cluttered.
# The direct text labels P0, Ps, Pg, Pl are likely clearer for this specific purpose.


# Row 2: Delay Embedding Scatter Plots
ax_embed_small = Axis(fig[2, 1],
    title=L"(B) Small Delay: $\tau = %$(τ_s_val)$",
    xlabel=L"$s(t)$", ylabel=L"$s(t+\tau)$", aspect=DataAspect())
ax_embed_good = Axis(fig[2, 2],
    title=L"(C) Appropriate Delay: $\tau = %$(τ_g_val)$",
    xlabel=L"$s(t)$", ylabel=L"$s(t+\tau)$", aspect=DataAspect())
ax_embed_large = Axis(fig[2, 3],
    title=L"(D) Poor Large Delay: $\tau = %$(τ_l_val) \approx 2\pi$",
    xlabel=L"$s(t)$", ylabel=L"$s(t+\tau)$", aspect=DataAspect())

scatter!(ax_embed_small, s_curr_small, s_del_small, markersize=3, color=(:crimson, 0.5))
scatter!(ax_embed_good, s_curr_good, s_del_good, markersize=3, color=(:forestgreen, 0.5))
scatter!(ax_embed_large, s_curr_large, s_del_large, markersize=3, color=(:darkorchid, 0.5))

common_lims = (-1.1, 1.1)
for ax_embed in [ax_embed_small, ax_embed_good, ax_embed_large]
    xlims!(ax_embed, common_lims)
    ylims!(ax_embed, common_lims)
end

# --- Layout Adjustments ---
rowgap!(fig.layout, 30)
colgap!(fig.layout, 25)

# You might need to adjust row height for the first row if annotations take more space
# rowsize!(fig.layout, 1, Auto()) # or Fixed(some_height)

# --- Save the Figure ---
# save("delay_embedding_illustration_with_tau_visualization.pdf", fig)
save("figures/delay_embedding_illustration_with_tau_visualization.png", fig, px_per_unit=2)

# --- Display the Figure ---
println("Makie script finished. If uncommented, figure saved.")
fig # Display the figure if in an interactive environment
