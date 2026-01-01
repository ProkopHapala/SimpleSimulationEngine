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
            verbosity: 0,
            paused: false,
            running: false
        };

        this.elements = {
            nParticles: document.getElementById('n-particles'),
            dt: document.getElementById('dt'),
            gravity: document.getElementById('gravity'),
            speedBase: document.getElementById('speed-base'),
            speedSpread: document.getElementById('speed-spread'),
            renderMode: document.getElementById('render-mode'),
            rotScatter: document.getElementById('rot-scatter'),
            verbosity: document.getElementById('verbosity'),
            btnRun: document.getElementById('btn-run'),
            btnReset: document.getElementById('btn-reset'),
            btnPause: document.getElementById('btn-pause'),
            fps: document.getElementById('fps-val')
        };

        this.init();
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
        this.elements.verbosity.addEventListener('input', (e) => {
            this.params.verbosity = parseInt(e.target.value);
        });
        this.elements.rotScatter.addEventListener('input', (e) => {
            this.params.rotScatterDeg = parseFloat(e.target.value);
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
    }

    updateFPS(fps) {
        this.elements.fps.textContent = Math.round(fps);
    }
}
