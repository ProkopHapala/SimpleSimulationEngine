// gui-builder.js - Dynamic GUI builder for generator parameters
export function buildParamControls(container, params, onChange, values = {}) {
    params.forEach(p => {
        const row = document.createElement('div');
        row.style.marginBottom = '8px';
        const label = document.createElement('label');
        label.textContent = p.label + ':';
        label.style.display = 'block';
        label.style.marginBottom = '2px';
        label.style.fontSize = '12px';
        row.appendChild(label);
        if (p.type === 'int') {
            const input = document.createElement('input');
            input.type = 'range';
            input.min = p.min; input.max = p.max; input.step = 1;
            input.value = values[p.name] ?? p.default;
            input.style.width = '100%';
            const valSpan = document.createElement('span');
            valSpan.textContent = input.value;
            valSpan.style.fontSize = '11px';
            valSpan.style.float = 'right';
            label.appendChild(valSpan);
            input.addEventListener('input', () => { valSpan.textContent = input.value; onChange(p.name, parseInt(input.value)); });
            row.appendChild(input);
        } else if (p.type === 'float') {
            const input = document.createElement('input');
            input.type = 'range';
            input.min = p.min; input.max = p.max; input.step = 0.1;
            input.value = values[p.name] ?? p.default;
            input.style.width = '100%';
            const valSpan = document.createElement('span');
            valSpan.textContent = input.value;
            valSpan.style.fontSize = '11px';
            valSpan.style.float = 'right';
            label.appendChild(valSpan);
            input.addEventListener('input', () => { valSpan.textContent = input.value; onChange(p.name, parseFloat(input.value)); });
            row.appendChild(input);
        }
        container.appendChild(row);
    });
}
