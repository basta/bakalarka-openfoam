// frontend/js/ui.js

// Define icons as SVG strings
const ICONS = {
    chevronRight: '<svg class="icon-chevron" viewBox="0 0 24 24"><path d="M9 18l6-6-6-6"/></svg>',
    figure: '<svg class="icon-item" viewBox="0 0 24 24"><path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm-1-13h2v6h-2zm0 8h2v2h-2z"/></svg>',
    axes: '<svg class="icon-item" viewBox="0 0 24 24"><path d="M3 3v18h18V3H3zm16 2v1.5L14.5 11l-3-3-6 6V5h14zm-14 12.5V13l6-6 3 3 4.5-4.5V17H5z"/></svg>',
    elements: '<svg class="icon-item" viewBox="0 0 24 24"><path d="M12 2L2 7l10 5 10-5-10-5zM2 17l10 5 10-5M2 12l10 5 10-5"></path></svg>',
    legend: '<svg class="icon-item" viewBox="0 0 24 24"><path d="M3 3h18v18H3V3zm2 2v14h14V5H5zm2 2h10v2H7V7zm0 4h10v2H7v-2zm0 4h10v2H7v-2z"/></svg>',
    title: '<svg class="icon-item" viewBox="0 0 24 24"><path d="M5 4v3h5.5v12h3V7H19V4H5z"/></svg>'
};

// Sample tree data structure. This would be dynamically generated or fetched in a real app.
const TREE_DATA = [
    {
        id: 'figureNode',
        label: 'Figure Options',
        icon: ICONS.figure,
        expanded: true,
        contentHTML: `
            <div class="form-group">
                <label for="fig-resolution-w" class="form-label">Resolution (W x H)</label>
                <div class="input-group-sm">
                    <input type="number" id="fig-resolution-w" placeholder="W" class="form-input">
                    <input type="number" id="fig-resolution-h" placeholder="H" class="form-input">
                </div>
            </div>
            <div class="form-group">
                <label for="fig-bg-color" class="form-label">Background Color</label>
                <input type="color" id="fig-bg-color" value="#333333" class="form-color-picker">
            </div>
            <div class="form-group">
                <label for="fig-font" class="form-label">Global Font</label>
                <select id="fig-font" class="form-select">
                    <option>Inter</option><option>Arial</option><option>System Default</option>
                </select>
            </div>
        `,
        children: []
    },
    {
        id: 'axesNode',
        label: 'Axes Configuration',
        icon: ICONS.axes,
        expanded: false,
        contentHTML: `
            <h4 class="font-medium text-inspector-text-secondary border-b border-inspector-border pb-1 mb-2 text-xs uppercase">X-Axis</h4>
            <div class="form-group">
                <label for="xaxis-label" class="form-label">Label</label>
                <input type="text" id="xaxis-label" placeholder="X-Axis Label" class="form-input">
            </div>
            <div class="form-group">
                <input type="checkbox" id="xaxis-grid-visible" class="form-checkbox">
                <label for="xaxis-grid-visible" class="checkbox-label">Grid Lines Visible</label>
            </div>
            <h4 class="font-medium text-inspector-text-secondary border-b border-inspector-border pb-1 my-2 text-xs uppercase">Y-Axis</h4>
            `,
        children: [
            {
                id: 'xAxisNode',
                label: 'X-Axis Details',
                icon: ICONS.axes,
                expanded: false,
                contentHTML: `
                    <div class="form-group">
                        <label for="xaxis-min" class="form-label">Min Limit</label>
                        <input type="number" id="xaxis-min" class="form-input">
                    </div>
                    <div class="form-group">
                        <label for="xaxis-max" class="form-label">Max Limit</label>
                        <input type="number" id="xaxis-max" class="form-input">
                    </div>
                `,
                children: []
            },
            {
                id: 'yAxisNode',
                label: 'Y-Axis Details',
                icon: ICONS.axes,
                expanded: false,
                contentHTML: `<div>Y-Axis specific deep nested controls...</div>`,
                children: []
            }
        ]
    },
    {
        id: 'elementsNode',
        label: 'Plot Elements',
        icon: ICONS.elements,
        expanded: false,
        contentHTML: `
            <h4 class="font-medium text-inspector-text-secondary border-b border-inspector-border pb-1 mb-2 text-xs uppercase">Line 1 (Series A)</h4>
            <div class="form-group">
                <label for="line1-color" class="form-label">Color</label>
                <input type="color" id="line1-color" value="#007ACC" class="form-color-picker">
            </div>
            <div class="form-group">
                <label for="line1-width" class="form-label">Linewidth</label>
                <input type="number" id="line1-width" value="1.5" step="0.1" class="form-input">
            </div>
        `,
        children: []
    },
    {
        id: 'legendNode',
        label: 'Legend',
        icon: ICONS.legend,
        expanded: false,
        contentHTML: `
            <div class="form-group">
                <input type="checkbox" id="legend-visible" class="form-checkbox" checked>
                <label for="legend-visible" class="checkbox-label">Visible</label>
            </div>
            <div class="form-group">
                <label for="legend-position" class="form-label">Position</label>
                <select id="legend-position" class="form-select">
                    <option value="northwest">North West</option>
                    <option value="center">Center</option>
                </select>
            </div>
        `,
        children: []
    },
    {
        id: 'titleNode',
        label: 'Plot Title',
        icon: ICONS.title,
        expanded: false,
        contentHTML: `
            <div class="form-group">
                <label for="plot-title-text" class="form-label">Text</label>
                <input type="text" id="plot-title-text" class="form-input" placeholder="Enter plot title">
            </div>
             <div class="form-group">
                <label for="plot-title-size" class="form-label">Font Size</label>
                <input type="number" id="plot-title-size" class="form-input" value="16">
            </div>
        `,
        children: []
    }
];

function createTreeNodeHTML(nodeData) {
    const hasChildren = nodeData.children && nodeData.children.length > 0;
    const li = document.createElement('li');
    li.classList.add('tree-node');
    if (nodeData.expanded) {
        li.classList.add('expanded');
    }
    li.dataset.nodeId = nodeData.id;

    const labelDiv = document.createElement('div');
    labelDiv.classList.add('tree-node-label');
    labelDiv.innerHTML = `
        ${ICONS.chevronRight}
        ${nodeData.icon || ''}
        <span>${nodeData.label}</span>
    `;
    labelDiv.addEventListener('click', () => {
        li.classList.toggle('expanded');
        // Potentially: highlight selected node, update a "current selection" state
        document.querySelectorAll('.tree-node.selected').forEach(n => n.classList.remove('selected'));
        li.classList.add('selected');
        console.log(`Selected node: ${nodeData.id} - ${nodeData.label}`);
    });

    li.appendChild(labelDiv);

    if (nodeData.contentHTML || hasChildren) {
        const childrenDiv = document.createElement('div');
        childrenDiv.classList.add('tree-node-children');

        if (nodeData.contentHTML) {
            const contentSection = document.createElement('div');
            contentSection.classList.add('control-section-content');
            contentSection.innerHTML = nodeData.contentHTML;
            childrenDiv.appendChild(contentSection);
        }

        if (hasChildren) {
            const nestedUl = document.createElement('ul');
            nodeData.children.forEach(childNode => {
                nestedUl.appendChild(createTreeNodeHTML(childNode));
            });
            childrenDiv.appendChild(nestedUl);
        }
        li.appendChild(childrenDiv);
    }
    return li;
}

function buildTree(data, parentElementId) {
    const parentElement = document.getElementById(parentElementId);
    if (!parentElement) {
        console.error(`Parent element with ID '${parentElementId}' not found for tree.`);
        return;
    }
    parentElement.innerHTML = ''; // Clear existing tree
    data.forEach(nodeData => {
        parentElement.appendChild(createTreeNodeHTML(nodeData));
    });
}

function setupSidebarToggle() {
    const controlPanel = document.getElementById('controlPanel');
    const toggleSidebarButton = document.getElementById('toggleSidebarButton');
    const toggleSidebarLink = document.getElementById('toggleSidebarLink');
    const sidebarToggleIcon = document.getElementById('sidebarToggleIcon');
    const panelRightCloseIconSVG = '<rect width="18" height="18" x="3" y="3" rx="2"></rect><path d="M15 3v18"></path><path d="m8 9 3 3-3 3"></path>';
    const panelRightOpenIconSVG = '<rect width="18" height="18" x="3" y="3" rx="2"></rect><path d="M15 3v18"></path><path d="m10 9-3 3 3 3"></path>';

    function toggleInspectorPanel() {
        if (!controlPanel) return;
        controlPanel.classList.toggle('collapsed');
        if (sidebarToggleIcon) {
            if (controlPanel.classList.contains('collapsed')) {
                sidebarToggleIcon.innerHTML = panelRightOpenIconSVG;
                if(toggleSidebarButton) toggleSidebarButton.title = "Open Inspector Panel";
            } else {
                sidebarToggleIcon.innerHTML = panelRightCloseIconSVG;
                if(toggleSidebarButton) toggleSidebarButton.title = "Close Inspector Panel";
            }
        }
    }

    if (toggleSidebarButton) {
        toggleSidebarButton.addEventListener('click', toggleInspectorPanel);
    }
    if (toggleSidebarLink) {
        toggleSidebarLink.addEventListener('click', (e) => {
            e.preventDefault();
            toggleInspectorPanel();
            const parentMenu = toggleSidebarLink.closest('.menu-bar-item');
            if (parentMenu) parentMenu.blur(); // Close dropdown
        });
    }
}

function setupMenuDropdowns() {
    document.addEventListener('click', function(event) {
        const openDropdowns = document.querySelectorAll('.menu-bar-item .dropdown-menu[style*="display: block"]');
        openDropdowns.forEach(dropdown => {
            if (!dropdown.parentElement.contains(event.target)) {
                dropdown.style.display = 'none';
            }
        });
    });

    document.querySelectorAll('.menu-bar-item').forEach(item => {
        item.addEventListener('click', function(event) {
            if (!event.target.closest('.dropdown-menu')) {
                const dropdown = this.querySelector('.dropdown-menu');
                if (dropdown) {
                    const currentDisplay = window.getComputedStyle(dropdown).display;
                    document.querySelectorAll('.menu-bar-item .dropdown-menu').forEach(d => {
                        if (d !== dropdown) d.style.display = 'none';
                    });
                    dropdown.style.display = currentDisplay === 'none' ? 'block' : 'none';
                }
            }
            event.stopPropagation();
        });
    });
}

// Expose functions to be called by app.js
window.UiModule = {
    buildTree,
    TREE_DATA,
    setupSidebarToggle,
    setupMenuDropdowns
};
