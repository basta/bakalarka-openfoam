# Interactive Makie Plot Editor: API Scheme Reference

## 1. Introduction

This document provides a detailed specification for the Application Programming Interface (API) and communication protocol for the "Interactive Makie Plot Editor." It is intended for developers working on either the frontend or backend components of the application. The primary mode of communication is via WebSockets, with messages formatted as JSON objects.

## 2. General Message Structure

All WebSocket messages, in both directions, will adhere to the following basic structure:

```json
{
  "type": "MESSAGE_TYPE_STRING",
  "payload": {
    // Message-specific data, varies by 'type'
  },
  "messageId": "UNIQUE_MESSAGE_ID" // Optional: for request-response tracking
}
```

## 3. Frontend to Backend Messages

These messages are initiated by the frontend (user interface) and sent to the Julia/Makie backend.

### 3.1. `REQUEST_INITIAL_PLOT`

* **Type:** `REQUEST_INITIAL_PLOT`
* **Payload:**
    * `plotIdentifier` (String, Optional): An identifier for a predefined plot or dataset if the application supports loading them.
    * `initialData` (Object, Optional): Data to be plotted if not loading a predefined plot. The structure of this object will depend on the plot types supported.
* **Description:** Sent by the frontend upon loading or when a user requests a new base plot. The backend is expected to respond with an `INITIAL_PLOT_STATE` message.

### 3.2. `SET_PLOT_ATTRIBUTE`

* **Type:** `SET_PLOT_ATTRIBUTE`
* **Payload:**
    * `targetId` (String, Required): A unique identifier for the plot element (e.g., "axis_1", "line_series_A", "figure_title") to be modified. These IDs are provided by the backend in the `INITIAL_PLOT_STATE` message.
    * `attributeClass` (String, Optional but Recommended): The class of the target element (e.g., "Axis", "Lines", "Scatter", "Legend"). This helps the backend interpret the attribute correctly.
    * `attributeName` (String, Required): The name of the Makie attribute to be changed (e.g., "title", "color", "linewidth", "xticklabelrotation").
    * `newValue` (Any, Required): The new value for the attribute. The data type depends on the specific attribute (e.g., string, number, boolean, array, object).
* **Description:** Sent when the user modifies a plot property through the UI. The backend should attempt to apply this change and respond with a `PLOT_UPDATE` message or an `ERROR_MESSAGE` if the change is invalid.

#### 3.2.1. Example `SET_PLOT_ATTRIBUTE` Payloads:

* **Changing Axis Title:**
    ```json
    {
      "targetId": "axis_main",
      "attributeClass": "Axis",
      "attributeName": "title",
      "newValue": "Temperature Over Time"
    }
    ```

* **Changing Line Color of a Series:**
    ```json
    {
      "targetId": "series_temp_A",
      "attributeClass": "Lines",
      "attributeName": "color",
      "newValue": "rgba(255, 0, 0, 1.0)" // Or "red"
    }
    ```

* **Changing X-Axis Label Font Size:**
    ```json
    {
      "targetId": "axis_main",
      "attributeClass": "Axis",
      "attributeName": "xlabelsize",
      "newValue": 14
    }
    ```

* **Changing Scatter Marker Style:**
    ```json
    {
      "targetId": "scatter_data_points",
      "attributeClass": "Scatter",
      "attributeName": "marker",
      "newValue": ":diamond" // Makie marker symbol
    }
    ```

* **Toggling Legend Visibility:**
    ```json
    {
      "targetId": "legend_primary",
      "attributeClass": "Legend",
      "attributeName": "visible",
      "newValue": false
    }
    ```

#### 3.2.2. Detailed Attribute Communications by Plot Element Type

The `SET_PLOT_ATTRIBUTE` message will be used to modify a wide range of attributes. The `attributeName` and `newValue` will correspond to Makie.jl's attribute system. Below are categories of elements and common attributes. The `targetId` in the message will point to an instance of one of these elements.

**A. Figure/Scene Attributes:**
*(Corresponds to Makie `Figure` and `Scene` attributes)*
* `backgroundcolor`: Color string (e.g., `"rgba(255,255,255,1)"`, `"#FEFEFE"`) or Makie color type.
* `resolution` / `size`: Array `[width, height]`.
* `figure_padding`: Number (uniform) or Array `[left, right, bottom, top]`.
* `font`: Default font for the figure.
* `fontsize`: Default font size.
* Lights (for 3D scenes, part of `Scene`): Array of light objects (e.g., `AmbientLight`, `PointLight`, `DirectionalLight`). Each light would have its own properties like `color`, `intensity`, `position`, `direction`.
* Example for adding/modifying a light:
```json
{
  "targetId": "scene_1_light_ambient", // Or a general scene targetId
  "attributeClass": "Scene", // Or "AmbientLight" if lights are addressable
  "attributeName": "lights", // Or specific light property
  "newValue": { "type": "AmbientLight", "color": "rgba(0.1,0.1,0.1,1)" }
}
```

**B. Axis Attributes (`Axis`, `Axis3`, `PolarAxis`):**
* **Titles & Subtitles:** `title`, `subtitle`, `titlefont`, `titlesize`, `titlecolor`, `titlealign`, `titlegap`, `subtitlefont`, `subtitlesize`, `subtitlecolor`, `subtitlealign`, `subtitlegap`, `titlevisible`, `subtitlevisible`.
* **Axis Labels (X, Y, Z, R, Theta prefixes):** `xlabel`, `ylabel`, `zlabel`, `xlabelfont`, `ylabelfont`, `zlabelfont`, `xlabelsize`, `ylabelsize`, `zlabelsize`, `xlabelcolor`, `ylabelcolor`, `zlabelcolor`, `xlabelpadding`, `ylabelpadding`, `xlabeloffset` (Axis3), `xlabelrotation`, `ylabelrotation`, `zlabelrotation`, `xlabelvisible`, `ylabelvisible`, `zlabelvisible`.
* **Ticks & Tick Labels (X, Y, Z, R, Theta prefixes):** `xticks`, `yticks`, `zticks` (values, WilkinsonTicks object, specific tick generators), `xtickformat`, `ytickformat`, `ztickformat` (format string or function), `xticklabelalign`, `yticklabelalign`, `zticklabelalign`, `xticklabelcolor`, `yticklabelcolor`, `zticklabelcolor`, `xticklabelfont`, `yticklabelfont`, `zticklabelfont`, `xticklabelpad`, `yticklabelpad`, `zticklabelpad`, `xticklabelrotation`, `yticklabelrotation`, `zticklabelrotation`, `xticklabelsize`, `yticklabelsize`, `zticklabelsize`, `xticklabelsvisible`, `yticklabelsvisible`, `zticklabelsvisible`, `xtickalign`, `ytickalign`, `ztickalign`, `xtickcolor`, `ytickcolor`, `ztickcolor`, `xticksize`, `yticksize`, `zticksize`, `xticksvisible`, `yticksvisible`, `zticksvisible`, `xtickwidth`, `ytickwidth`, `ztickwidth`.
* **Minor Ticks & Grids (X, Y, Z, R, Theta prefixes):** `xminorticks`, `yminorticks`, `zminorticks` (e.g., `IntervalsBetween(5)`), and corresponding `xminorgridcolor`, `xminorgridstyle`, `xminorgridvisible`, `xminorgridwidth`, `xminortickalign`, `xminortickcolor`, `xminorticksize`, `xminorticksvisible`, `xminortickwidth`.
* **Grid Lines (X, Y, Z, R, Theta prefixes):** `xgridcolor`, `ygridcolor`, `zgridcolor`, `xgridstyle`, `ygridstyle`, `zgridstyle` (e.g., `:dash`), `xgridvisible`, `ygridvisible`, `zgridvisible`, `xgridwidth`, `ygridwidth`, `zgridwidth`.
* **Spines (`Axis`):** `bottomspinecolor`, `leftspinecolor`, `rightspinecolor`, `topspinecolor`, `bottomspinevisible`, `leftspinevisible`, `rightspinevisible`, `topspinevisible`, `spinewidth`.
* **Spines (`Axis3`):** `xspinesvisible`, `yspinesvisible`, `zspinesvisible`, `xspinewidth`, `yspinewidth`, `zspinewidth`, and individual spine colors (e.g., `xspinecolor_1`, `yspinecolor_2`).
* **Axis Panels (`Axis3`):** `xypanelcolor`, `yzpanelcolor`, `xzpanelcolor`, `xypanelvisible`, etc.
* **General Axis Properties:** `backgroundcolor`, `aspect` (Number, `DataAspect()`, `AxisAspect()`), `autolimitaspect`, `limits` (Tuple or `Makie.Automatic()`), `xreversed`, `yreversed`, `zreversed`, `xscale`, `yscale`, `zscale` (e.g., `log10`, `Makie.logit`).
* **`PolarAxis` Specific:** `rlimits`, `thetalimits`, `rtickformat`, `thetatickformat`, `radius_at_origin`, `direction`, `theta_as_x`.
* **`Axis3` Specific:** `viewmode` (`:fit`, `:stretch`, etc.), `perspectiveness`, `azimuth`, `elevation`, `show_axis`, `protrusions`.

**C. Plot Object Attributes (e.g., `Lines`, `Scatter`, `BarPlot`, `Heatmap`, `Surface`, `Contourf`, `Mesh`, `Volume`, `Arrows`):**
* **Common Visuals:** `color` (single, array, or number for colormapping), `colormap` (name, `Makie.Categorical`, or array of colors), `colorrange` (tuple `(min, max)`), `colorscale` (`identity`, `log10`), `alpha` (0-1), `visible` (boolean), `fxaa` (boolean), `depth_shift` (number), `space` (`:data`, `:pixel`), `inspectable`.
* **Line-Specific (`Lines`, `ScatterLines`, `Contour`, `ABLines`, `HLines`, `VLines`):** `linestyle` (`:solid`, `:dot`, `Makie.Linestyle` object), `linewidth`, `linecap`, `joinstyle`.
* **Marker-Specific (`Scatter`, `ScatterLines`, `MeshScatter`):** `marker` (Symbol, Char, `Makie.BezierPath`, `Makie.FastPixel()`), `markersize` (Number, `Vec`), `markercolor` (often from `color`), `strokecolor` (for marker outline), `strokewidth` (for marker outline), `marker_offset` (`Vec`), `rotation` (Number, `Vec3f`, `Quaternion`, `Billboard()`), `transform_marker`.
* **`BarPlot` Specific:** `fillto`, `width`, `gap`, `dodge`, `stack`, `direction` (`:x`, `:y`), `strokecolor`, `strokewidth`, `bar_labels`, `label_font`, `label_size`, etc.
* **`Heatmap`, `Image`, `Surface` Specific:** `interpolate`, `nan_color`, `highclip`, `lowclip`, `uv_transform` (for texture mapping).
* **`Volume` Specific:** `algorithm` (`:mip`, `:iso`), `isovalue`, `isorange`, `absorption`.
* **`Mesh` Specific:** `shading`, `diffuse`, `specular`, `shininess`, `backlight`, `matcap`.
* **`Arrows` Specific:** `arrowhead`, `arrowsize`, `arrowcolor`, `linecolor`, `lengthscale`, `normalize`, `align`.
* **`Annotations`/`text` Specific:** `text` (string or array), `font`, `fontsize`, `align`, `rotation`, `justification`, `offset`, `word_wrap_width`.

**D. Legend Attributes (`Legend`):**
* `title`, `titlefont`, `titlesize`, `titlecolor`, `titlealign`, `titlevisible`, `titlegap`.
* `labels` (Array of strings, though often derived).
* `labelfont`, `labelsize`, `labelcolor`.
* `bgcolor` (or `backgroundcolor`), `framecolor`, `framewidth`, `framevisible`.
* `orientation` (`:horizontal`, `:vertical`).
* `tellwidth`, `tellheight`.
* `halign`, `valign` (for layout).
* `padding`, `margin`, `colgap`, `rowgap`, `groupgap`.
* Entry-specific attributes (templates): `patchcolor`, `patchsize`, `linecolor`, `linestyle`, `markercolor`, `marker`.

**E. Colorbar Attributes (`Colorbar`):**
* `colormap`, `colorrange`, `colorscale`.
* `label`, `labelfont`, `labelsize`, `labelcolor`, `labelpadding`, `labelrotation`, `labelvisible`.
* Ticks & Labels: `ticks`, `tickformat`, `ticklabelalign`, `ticklabelcolor`, `ticklabelfont`, `ticklabelsize`, `ticklabelsvisible`.
* `vertical`, `flipaxis`, `highclip`, `lowclip`, `spinewidth`.
* `halign`, `valign`, `width`, `height`.

### 3.3. `APPLY_THEME`

* **Type:** `APPLY_THEME`
* **Payload:**
    * `themeName` (String): Identifier for a predefined Makie theme (e.g., "theme_dark()", "theme_ggplot2()").
    * OR
    * `themeAttributes` (Object): A JSON object representing a custom theme, mirroring Makie's `Theme` structure.
* **Description:** Instructs the backend to apply a new theme to the plot. The backend should update the plot and respond with `PLOT_UPDATE` and potentially `INITIAL_PLOT_STATE` if inspectable attributes change significantly.

### 3.4. `EXPORT_PLOT`

* **Type:** `EXPORT_PLOT`
* **Payload:**
    * `format` (String, Required): The desired export format ("PNG", "SVG", "PDF", "JULIA_CODE").
    * `options` (Object, Optional): Format-specific options.
        * For image formats ("PNG"): `resolution` (e.g., `[1920, 1080]`), `px_per_unit`.
        * For vector formats ("SVG", "PDF"): `pt_per_unit` (for CairoMakie).
* **Description:** Requests the backend to export the current state of the plot. The backend should respond with `EXPORT_RESULT` or `ERROR_MESSAGE`.

### 3.5. Future Messages (as per Design Document)

* **`UNDO_ACTION`**:
    * **Type:** `UNDO_ACTION`
    * **Payload:** (Potentially empty, or could include a transaction ID if finer control is needed)
    * **Description:** Requests the backend to revert the last user-initiated change. Requires backend to maintain a history of states or reversible commands.
* **`REDO_ACTION`**:
    * **Type:** `REDO_ACTION`
    * **Payload:** (Similar to `UNDO_ACTION`)
    * **Description:** Requests the backend to reapply the last undone action.
* **`LOAD_DATA`**:
    * **Type:** `LOAD_DATA`
    * **Payload:**
        * `dataSourceType`: (String, e.g., "CSV_URL", "JSON_INLINE", "JULIA_VARIABLE_NAME")
        * `sourceDetails`: (Object, content depends on `dataSourceType`)
            * If "CSV_URL": `{ "url": "http://example.com/data.csv" }`
            * If "JSON_INLINE": `{ "data": [{ "x": 1, "y": 2 }, ...] }`
        * `plotMapping`: (Object, Optional) How to map loaded data to plot arguments (e.g., `{ "x_column": "column_A", "y_column": "column_B" }`).
    * **Description:** Requests the backend to load new data and update the plot. Backend responds with `PLOT_UPDATE` and `INITIAL_PLOT_STATE`.

## 4. Backend to Frontend Messages

These messages are initiated by the Julia/Makie backend and sent to the frontend.

### 4.1. `INITIAL_PLOT_STATE`

* **Type:** `INITIAL_PLOT_STATE`
* **Payload:**
    * `plotData` (String or Object, Required): The visual representation of the plot. This could be an SVG string, a base64 encoded PNG, or a more complex scene description if using WGLMakie/JSServe.
    * `elements` (Array, Required): An array of objects, each describing an interactable plot element. This is vital for populating the frontend's inspector panel.
        ```json
        [
          {
            "id": "figure_main",         // Unique ID for this element
            "type": "Figure",            // Makie type/class (e.g., "Figure", "Axis", "Lines")
            "name": "Main Figure",       // User-friendly display name
            "attributes": [              // Array of modifiable attributes for this element
              {
                "name": "backgroundcolor", // Technical attribute name
                "value": "rgba(255,255,255,1)", // Current value
                "type": "color",         // UI hint: "string", "number", "boolean", "color", "enum", "font"
                "label": "Background Color", // User-friendly label for UI
                "description": "Sets the background color of the entire figure." // Optional tooltip
              },
              {
                "name": "resolution",
                "value": [800, 600],
                "type": "dimensions_wh", // e.g. [width, height]
                "label": "Resolution (WxH)"
              }
              // ... more attributes
            ],
            "children": ["axis_1", "global_legend"] // IDs of child elements for tree view construction
          },
          {
            "id": "axis_1",
            "type": "Axis",
            "name": "Primary Data Axis",
            "attributes": [
              { "name": "title", "value": "Sales Data", "type": "string", "label": "Axis Title" },
              { "name": "xlabel", "value": "Month", "type": "string", "label": "X-Axis Label" },
              {
                "name": "xlabelsize", "value": 12, "type": "number", "label": "X Label Size",
                "min": 6, "max": 30, "step": 1 // UI hints for a slider/number input
              },
              {
                "name": "xgridvisible", "value": true, "type": "boolean", "label": "Show X Grid"
              },
              {
                "name": "xtickformat", "value": "{:.0f}", "type": "string", "label": "X Tick Format"
              },
              {
                "name": "xscale", "value": "identity", "type": "enum",
                "options": ["identity", "log10", "sqrt"], "label": "X-Axis Scale"
              }
              // ... many more attributes for Makie.Axis
            ],
            "children": ["line_plot_sales", "scatter_plot_outliers"]
          },
          {
            "id": "line_plot_sales",
            "type": "Lines", // Could be more specific like "Makie.Lines"
            "name": "Monthly Sales Trend",
            "attributes": [
              { "name": "color", "value": "blue", "type": "color", "label": "Line Color" },
              {
                "name": "linewidth", "value": 2.0, "type": "number", "label": "Line Width",
                "min": 0.5, "max": 10, "step": 0.1
              },
              {
                "name": "linestyle", "value": "solid", "type": "enum",
                "options": ["solid", "dot", "dash", "dashdot"], "label": "Line Style"
              }
              // ... other attributes from Makie.Lines
            ]
            // No children for a primitive plot object
          }
          // ... other elements like Legend, Colorbar, individual Lights etc.
        ]
        ```
    * `availableThemes` (Array of Strings, Optional): A list of theme names the backend supports (e.g., `["theme_light", "theme_dark"]`).
    * `availableColormaps` (Array of Strings or Objects, Optional): A list of colormap names or representations.
* **Description:** Sent by the backend in response to `REQUEST_INITIAL_PLOT` or when a significant structural change occurs (e.g., after `APPLY_THEME` if it redefines inspectable elements). The `attributes` array within each element should ideally include metadata for the UI (like `type`, `label`, `min`, `max`, `step`, `options` for enums) to help the frontend render appropriate controls.

### 4.2. `PLOT_UPDATE`

* **Type:** `PLOT_UPDATE`
* **Payload:**
    * `plotData` (String or Object, Required): The updated visual representation of the plot (e.g., SVG string, base64 PNG, or WGLMakie/JSServe scene modifications).
    * `updatedAttributes` (Object, Optional): A more granular update for the inspector panel. This object maps `targetId`s to an array of their attributes that have changed. This is more efficient than resending the entire `INITIAL_PLOT_STATE` for minor changes.
        ```json
        {
          "axis_1": [ // targetId
            { "name": "title", "value": "New Plot Title Updated" },
            { "name": "xlabelsize", "value": 16 }
          ],
          "line_plot_sales": [ // targetId
            { "name": "color", "value": "green" }
          ]
        }
        ```
* **Description:** Sent by the backend whenever the plot's visual appearance or an inspectable attribute value changes as a result of a `SET_PLOT_ATTRIBUTE`, `APPLY_THEME`, or other backend-initiated operations.

### 4.3. `EXPORT_RESULT`

* **Type:** `EXPORT_RESULT`
* **Payload:**
    * `format` (String, Required): The format of the exported data (e.g., "PNG", "SVG", "PDF", "JULIA_CODE").
    * `data` (String, Required): The exported plot data.
        * For "PNG": Base64 encoded image string.
        * For "SVG", "JULIA_CODE": Plain text string.
        * For "PDF": Base64 encoded string or a URL for direct download if too large.
    * `filename` (String, Optional): A suggested filename for the download (e.g., "my_plot.png").
* **Description:** Sent by the backend after successfully processing an `EXPORT_PLOT` request.

### 4.4. `STATUS_MESSAGE`

* **Type:** `STATUS_MESSAGE`
* **Payload:**
    * `message` (String, Required): A human-readable status update (e.g., "Plotting complete.", "Applying theme...").
    * `level` (String, Optional): Severity level ("INFO", "WARNING", "SUCCESS"). Default is "INFO".
    * `duration` (Number, Optional): Time in milliseconds the operation took, if relevant.
* **Description:** Used for general feedback from the backend to the frontend, not necessarily tied to an error.

### 4.5. `ERROR_MESSAGE`

* **Type:** `ERROR_MESSAGE`
* **Payload:**
    * `message` (String, Required): A human-readable description of the error.
    * `details` (String, Optional): More technical details about the error, or a stack trace (potentially summarized).
    * `targetId` (String, Optional): If the error is related to a specific plot element from a `SET_PLOT_ATTRIBUTE` message.
    * `attributeName` (String, Optional): If the error is related to a specific attribute from a `SET_PLOT_ATTRIBUTE` message.
* **Description:** Sent by the backend when an error occurs (e.g., invalid attribute value, failure in plot generation, export error).

## 5. Element Identification (`targetId`)

* **Generation:** The backend is responsible for generating unique and stable IDs for all significant, modifiable plot elements. This includes the Figure, Scenes, Axes, individual plot objects (lines, scatters), Legends, Colorbars, and potentially even sub-components like Titles or specific Labels if they are to be independently targetable.
* **Provision:** These IDs are provided to the frontend primarily through the `elements` array in the `INITIAL_PLOT_STATE` message.
* **Usage:** The frontend uses these `targetId`s in messages like `SET_PLOT_ATTRIBUTE` to specify which element's property needs to be changed.
* **Hierarchy:** While IDs can be flat strings (e.g., UUIDs), a hierarchical naming convention (e.g., "figure1.axis1.series_A") could be considered for easier debugging. The `children` property in the `INITIAL_PLOT_STATE.elements` array allows the frontend to reconstruct the hierarchical structure for display (e.g., in a tree view inspector) regardless of the ID format.

## 6. Role of JSServe.jl / WGLMakie.jl (Considerations)

As outlined in the project's design document:

* **Potential for Reduced Messaging:** If JSServe.jl and WGLMakie.jl are utilized, JSServe's built-in synchronization of Julia `Observable`s with frontend JavaScript elements can handle many direct attribute modifications. This could reduce the volume of explicit `SET_PLOT_ATTRIBUTE` / `PLOT_UPDATE` messages for simple, directly mapped attributes.
* **Continued Importance of Explicit Protocol:** The explicit WebSocket protocol defined herein remains crucial for:
    * **Initial State Synchronization:** Loading the initial plot and comprehensively populating the frontend inspector panel (`REQUEST_INITIAL_PLOT`, `INITIAL_PLOT_STATE`).
    * **Complex Operations:** Actions that are not simple attribute assignments, such as applying entire themes (`APPLY_THEME`) or handling plot exports (`EXPORT_PLOT`).
    * **Inspector Panel Management:** Ensuring the inspector panel accurately reflects the current state of all attributes, especially those not directly tied to a simple `Observable` or those requiring backend validation.
    * **Non-Standard Interactions:** Custom interactions or modifications that go beyond what JSServe's default synchronization handles.
    * **Error Handling and Status Reporting:** Providing clear feedback to the user (`ERROR_MESSAGE`, `STATUS_MESSAGE`).
    * **Code Generation:** The future goal of generating Julia Makie code.
* **WGLMakie Efficiency:** WGLMakie.jl can make plot updates more efficient by sending scene modifications (diffs) rather than full plot re-renders (images/SVG) for many interactions, which would be reflected in the `plotData` field of `PLOT_UPDATE` messages.

This API scheme aims to be comprehensive enough to support the described features while allowing for optimization if JSServe.jl/WGLMakie.jl handles some aspects of real-time synchronization.
