import { computeCoefficients } from './AeroSurface.js';
export class UI {
    constructor(callbacks) {
        this.callbacks = callbacks;
        this.params = {
            nParticles: 16,
            dt: 0.01,
            gravity: -9.81,
            rotScatterDeg: 20,
            speedBase: 1.0,
            speedSpread: 1.0,
            comOffset: 0.5,
            showGPU: true,
            showCPU: true,
            showCPUMesh: true,
            showCPULines: true,
            headSize: 0.02,
            verbosity: 0,
            paused: false,
            running: false,
            windTunnelMode: false,
            substeps: 1,
            plotType: 'aoa',
            plotWing: 0
        };

        this.elements = {
            nParticles: document.getElementById('n-particles'),
            dt: document.getElementById('dt'),
            gravity: document.getElementById('gravity'),
            speedBase: document.getElementById('speed-base'),
            speedSpread: document.getElementById('speed-spread'),
            comOffset: document.getElementById('com-offset'),
            showGPU: document.getElementById('show-gpu'),
            showCPU: document.getElementById('show-cpu'),
            showCPUMesh: document.getElementById('show-cpu-mesh'),
            showCPULines: document.getElementById('show-cpu-lines'),
            headSize: document.getElementById('head-size'),
            renderMode: document.getElementById('render-mode'),
            rotScatter: document.getElementById('rot-scatter'),
            verbosity: document.getElementById('verbosity'),
            btnRun: document.getElementById('btn-run'),
            btnReset: document.getElementById('btn-reset'),
            btnPause: document.getElementById('btn-pause'),
            windTunnelMode: document.getElementById('wind-tunnel-mode'),
            wingGeometry: document.getElementById('wing-geometry'),
            aeroParams: document.getElementById('aero-params'),
            aeroPlot: document.getElementById('aero-plot'),
            substeps: document.getElementById('substeps'),
            plotType: document.getElementById('plot-type'),
            plotWing: document.getElementById('plot-wing'),
            fps: document.getElementById('fps-val')
        };

        this.init();
    }

    setWingOptions(count) {
        const select = this.elements.plotWing;
        if (!select) return;
        select.innerHTML = '';
        for (let i = 0; i < count; i++) {
            const opt = document.createElement('option');
            opt.value = i;
            opt.textContent = `Wing ${i + 1}`;
            select.appendChild(opt);
        }
        if (count > 0) {
            this.params.plotWing = Math.min(this.params.plotWing, count - 1);
            select.value = this.params.plotWing;
        }
    }

    init() {
        this.elements.nParticles.addEventListener('change', (e) => {
            this.params.nParticles = parseInt(e.target.value);
            if (this.callbacks.onReset) this.callbacks.onReset();
        });
        this.elements.dt.addEventListener('input', (e) => {
            this.params.dt = parseFloat(e.target.value);
        });
        this.elements.gravity.addEventListener('input', (e) => {
            this.params.gravity = parseFloat(e.target.value);
        });
        this.elements.speedBase.addEventListener('input', (e) => {
            this.params.speedBase = parseFloat(e.target.value);
        });
        this.elements.speedSpread.addEventListener('input', (e) => {
            this.params.speedSpread = parseFloat(e.target.value);
        });
        this.elements.comOffset.addEventListener('input', (e) => {
            this.params.comOffset = parseFloat(e.target.value);
            if (this.callbacks.onReset) this.callbacks.onReset();
        });
        this.elements.showGPU.addEventListener('change', (e) => {
            this.params.showGPU = e.target.checked;
            if (this.callbacks.onToggleShowGPU) this.callbacks.onToggleShowGPU(this.params.showGPU);
        });
        this.elements.showCPU.addEventListener('change', (e) => {
            this.params.showCPU = e.target.checked;
            if (this.callbacks.onToggleShowCPU) this.callbacks.onToggleShowCPU(this.params.showCPU);
        });
        this.elements.showCPUMesh.addEventListener('change', (e) => {
            this.params.showCPUMesh = e.target.checked;
            if (this.callbacks.onToggleCPUDebug) this.callbacks.onToggleCPUDebug();
        });
        this.elements.showCPULines.addEventListener('change', (e) => {
            this.params.showCPULines = e.target.checked;
            if (this.callbacks.onToggleCPUDebug) this.callbacks.onToggleCPUDebug();
        });
        this.elements.headSize.addEventListener('input', (e) => {
            this.params.headSize = parseFloat(e.target.value);
            if (this.callbacks.onHeadSizeChange) this.callbacks.onHeadSizeChange(this.params.headSize);
        });
        this.elements.verbosity.addEventListener('input', (e) => {
            this.params.verbosity = parseInt(e.target.value);
        });
        this.elements.rotScatter.addEventListener('input', (e) => {
            this.params.rotScatterDeg = parseFloat(e.target.value);
        });
        this.elements.substeps.addEventListener('input', (e) => {
            const v = parseInt(e.target.value);
            this.params.substeps = Math.max(1, isNaN(v) ? 1 : v);
        });
        this.elements.renderMode.addEventListener('change', (e) => {
            if (this.callbacks.onRenderModeChange) this.callbacks.onRenderModeChange(e.target.value);
        });
        this.elements.btnRun.addEventListener('click', () => {
            this.params.running = !this.params.running;
            this.elements.btnRun.textContent = this.params.running ? "Run (On)" : "Run (Off)";
            if (this.callbacks.onRunToggle) this.callbacks.onRunToggle(this.params.running);
        });
        this.elements.btnReset.addEventListener('click', () => {
            if (this.callbacks.onReset) this.callbacks.onReset();
        });
        this.elements.btnPause.addEventListener('click', () => {
            this.params.paused = !this.params.paused;
            this.elements.btnPause.textContent = this.params.paused ? "Resume" : "Pause";
        });

        this.elements.windTunnelMode.addEventListener('change', (e) => {
            this.params.windTunnelMode = e.target.checked;
        });

        this.elements.plotType.addEventListener('change', (e) => {
            this.params.plotType = e.target.value;
        });
        this.elements.plotWing.addEventListener('change', (e) => {
            const idx = parseInt(e.target.value);
            this.params.plotWing = isNaN(idx) ? 0 : idx;
        });

        // Initialize foldable panels
        document.querySelectorAll('.fold-toggle').forEach(button => {
            button.addEventListener('click', () => {
                const content = button.nextElementSibling;
                content.style.display = content.style.display === 'block' ? 'none' : 'block';
            });
        });

        this.initTextAreas();
    }

    initTextAreas() {
        const defaultGeom = `# px py pz  tx ty tz  nx ny nz
# Main Wings (slightly behind CoM)
-1.5 0.0 -0.2   0.0 0.0 1.0   0.0 1.0 0.0
 1.5 0.0 -0.2   0.0 0.0 1.0   0.0 1.0 0.0
# Elevator (back)
 0.0 0.0 -2.0   0.0 0.0 1.0   0.0 1.0 0.0
# Rudder (back, perpendicular)
 0.0 0.0 -2.0   0.0 0.0 1.0   1.0 0.0 0.0`;
        const defaultAero = `# area CD0 dCD dCDS dCL dCLS sStall wStall
1.0 0.02 0.9 0.9 6.28 2.82 0.16 0.08
1.0 0.02 0.9 0.9 6.28 2.82 0.16 0.08
0.6 0.02 0.9 0.9 6.28 2.82 0.16 0.08
0.4 0.02 0.9 0.9 6.28 2.82 0.16 0.08`;
        
        this.elements.wingGeometry.value = defaultGeom;
        this.elements.aeroParams.value = defaultAero;

        this.elements.wingGeometry.addEventListener('change', () => {
            if (this.callbacks.onReset) this.callbacks.onReset();
        });
        this.elements.aeroParams.addEventListener('change', () => {
            if (this.callbacks.onReset) this.callbacks.onReset();
        });
    }

    updateFPS(fps) {
        this.elements.fps.textContent = Math.round(fps);
    }

    drawPlots(bodyDef) {
        const canvas = this.elements.aeroPlot;
        const ctx = canvas.getContext('2d');
        const w = canvas.width;
        const h = canvas.height;

        ctx.clearRect(0, 0, w, h);
        
        if (!bodyDef || !bodyDef.points || bodyDef.points.length === 0) return;

        const wingIdx = Math.min(this.params.plotWing || 0, bodyDef.points.length - 1);
        const p = bodyDef.points[wingIdx];
        const params = p.params;

        const pad = 8;
        ctx.strokeStyle = '#444';
        ctx.beginPath();
        ctx.moveTo(pad, h - pad); ctx.lineTo(w - pad, h - pad); // X axis (positive right)
        ctx.moveTo(pad, h - pad); ctx.lineTo(pad, pad);         // Y axis (positive up)
        ctx.stroke();

        ctx.lineWidth = 2;
        const steps = w;

        const computeLD = (AoA) => computeCoefficients(params, AoA);

        if (this.params.plotType === 'polar') {
            // Polar: x = CD*10, y = CL (positive quadrant)
            const samples = [];
            let maxX = 0, maxY = 0;
            for (let i = 0; i < steps; i++) {
                const AoA = ((i / (steps - 1)) - 0.5) * Math.PI;
                const { CL, CD } = computeLD(AoA);
                const x = Math.max(0, CD * 10);
                const y = Math.max(0, CL);
                samples.push({ x, y });
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
            }
            maxX = maxX || 1;
            maxY = maxY || 1;

            ctx.strokeStyle = '#0ff';
            ctx.beginPath();
            samples.forEach((s, i) => {
                const px = pad + (s.x / maxX) * (w - 2 * pad);
                const py = (h - pad) - (s.y / maxY) * (h - 2 * pad);
                if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
            });
            ctx.stroke();

            ctx.fillStyle = '#fff';
            ctx.font = '10px monospace';
            ctx.fillText('Polar: x=CD*10, y=CL', pad + 2, pad + 10);
        } else {
            // AoA plots: x = |AoA| mapped to [0..pi/2], y = CL/CD (clamped to +)
            const samplesCL = [];
            const samplesCD = [];
            let maxY = 0;
            for (let i = 0; i < steps; i++) {
                const AoA = ((i / (steps - 1)) - 0.5) * Math.PI; // -90..90
                const { CL, CD } = computeLD(AoA);
                const x = Math.abs(AoA) / (0.5 * Math.PI); // 0..1 for 0..90deg
                const yCL = Math.max(0, CL);
                const yCD = Math.max(0, CD);
                samplesCL.push({ x, y: yCL });
                samplesCD.push({ x, y: yCD });
                maxY = Math.max(maxY, yCL, yCD);
            }
            maxY = maxY || 1;

            const drawSeries = (data, color) => {
                ctx.strokeStyle = color;
                ctx.beginPath();
                data.forEach((s, i) => {
                    const px = pad + s.x * (w - 2 * pad);
                    const py = (h - pad) - (s.y / maxY) * (h - 2 * pad);
                    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
                });
                ctx.stroke();
            };

            drawSeries(samplesCL, '#44f'); // CL
            drawSeries(samplesCD, '#f44'); // CD

            ctx.fillStyle = '#fff';
            ctx.font = '10px monospace';
            ctx.fillText('AoA vs CL/CD (positive quadrant)', pad + 2, pad + 10);
        }
    }
}
