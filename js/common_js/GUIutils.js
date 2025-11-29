"use strict";



export const GUIutils = {

    /**
     * Create a container div with flex layout
     */
    box: function(parent, vertical = true) {
        const div = document.createElement('div');
        div.style.display = 'flex';
        div.style.flexDirection = vertical ? 'column' : 'row';
        div.style.gap = '4px';
        if (parent) parent.appendChild(div);
        return div;
    },

    /**
     * Create a label
     */
    label: function(parent, text) {
        const el = document.createElement('span');
        el.textContent = text;
        el.className = 'compact-label'; // Use existing class if available, or style inline
        if (parent) parent.appendChild(el);
        return el;
    },

    /**
     * Create a number input (spin box)
     */
    spinBox: function(parent, value, min, max, step, decimals = 2, callback = null) {
        const input = document.createElement('input');
        input.type = 'number';
        input.value = value;
        if (min !== undefined) input.min = min;
        if (max !== undefined) input.max = max;
        if (step !== undefined) input.step = step;
        input.className = 'compact-input';
        
        // Style for consistency
        input.style.width = '60px';

        if (callback) {
            input.addEventListener('input', () => callback(parseFloat(input.value)));
            input.addEventListener('change', () => callback(parseFloat(input.value)));
        }

        if (parent) parent.appendChild(input);
        return input;
    },

    /**
     * Create a text input
     */
    text: function(parent, value, callback = null) {
        const input = document.createElement('input');
        input.type = 'text';
        input.value = value;
        input.className = 'compact-input';
        input.style.width = '80px';
        
        if (callback) {
            input.addEventListener('change', () => callback(input.value));
        }

        if (parent) parent.appendChild(input);
        return input;
    },

    /**
     * Create a checkbox
     */
    checkBox: function(parent, text, checked, callback = null) {
        const label = document.createElement('label');
        label.style.display = 'flex';
        label.style.alignItems = 'center';
        label.style.gap = '4px';
        
        const input = document.createElement('input');
        input.type = 'checkbox';
        input.checked = checked;
        
        if (callback) {
            input.addEventListener('change', () => callback(input.checked));
        }

        label.appendChild(input);
        label.appendChild(document.createTextNode(text));

        if (parent) parent.appendChild(label);
        return { label, input };
    },

    /**
     * Create a button
     */
    button: function(parent, text, callback = null) {
        const btn = document.createElement('button');
        btn.textContent = text;
        btn.className = 'small-btn'; // Use existing class
        
        if (callback) {
            btn.addEventListener('click', callback);
        }

        if (parent) parent.appendChild(btn);
        return btn;
    },

    /**
     * Create a select dropdown
     */
    select: function(parent, options, selected, callback = null) {
        const sel = document.createElement('select');
        sel.className = 'compact-select';

        for (const key in options) {
            const opt = document.createElement('option');
            opt.value = options[key]; // Value
            opt.textContent = key;    // Display text
            if (options[key] === selected) opt.selected = true;
            sel.appendChild(opt);
        }

        if (callback) {
            sel.addEventListener('change', () => callback(sel.value));
        }

        if (parent) parent.appendChild(sel);
        return sel;
    },

    /**
     * Create a group box (fieldset-like)
     */
    group: function(parent, title) {
        const div = document.createElement('div');
        div.className = 'panel-section'; // Use existing class for consistent styling
        
        if (title) {
            const h3 = document.createElement('h3');
            h3.textContent = title;
            h3.style.margin = '2px 0 5px 0';
            h3.style.fontSize = '1em';
            div.appendChild(h3);
        }
        
        if (parent) parent.appendChild(div);
        return div;
    },

    /**
     * Data-driven GUI creation
     * @param {HTMLElement} parent - Container element
     * @param {Object} param_specs - Dictionary of parameter specifications
     * @param {Function} on_change - Global callback when any parameter changes (optional)
     * @returns {Object} - Object containing references to widgets and a getValues() function
     */
    create_gui: function(parent, param_specs, on_change = null) {
        const widgets = {};
        const groups = {};
        
        // Helper to get or create group
        const getGroup = (groupName) => {
            if (!groupName) return parent;
            if (!groups[groupName]) {
                groups[groupName] = this.group(parent, groupName);
            }
            return groups[groupName];
        };

        for (const name in param_specs) {
            const spec = param_specs[name];
            const container = getGroup(spec.group);
            
            // Row for label + widget
            const row = this.box(container, false);
            row.style.alignItems = 'center';
            row.style.justifyContent = 'space-between';
            row.style.marginBottom = '2px';

            this.label(row, name + ":");

            let widget;
            const handleChange = (val) => {
                if (on_change) on_change(name, val);
            };

            if (spec.widget === 'double' || spec.widget === 'int') {
                const decimals = spec.widget === 'int' ? 0 : (spec.decimals || 2);
                widget = this.spinBox(row, spec.value, spec.range[0], spec.range[1], spec.step, decimals, handleChange);
            } else if (spec.widget === 'bool') {
                const res = this.checkBox(row, '', spec.value, handleChange);
                widget = res.input;
            } else if (spec.widget === 'text') {
                widget = this.text(row, spec.value, handleChange);
            }

            widgets[name] = widget;
        }

        return {
            widgets: widgets,
            getValues: () => {
                const vals = {};
                for (const name in widgets) {
                    const w = widgets[name];
                    if (w.type === 'checkbox') vals[name] = w.checked;
                    else if (w.type === 'number') vals[name] = parseFloat(w.value);
                    else vals[name] = w.value;
                }
                return vals;
            },
            setValues: (vals) => {
                for (const name in vals) {
                    if (widgets[name]) {
                        const w = widgets[name];
                        if (w.type === 'checkbox') w.checked = vals[name];
                        else w.value = vals[name];
                    }
                }
            }
        };
    }
};

// Optional browser global for legacy code / debugging
if (typeof window !== 'undefined') {
    window.GUIutils = GUIutils;
}
