export class UI {
    constructor(callbacks) {
        this.callbacks = callbacks;
        this.params = {
            nParticles: 1,
            dt: 0.01,
            gravity: 0.0,
            paused: false
        };

        this.elements = {
            nParticles: document.getElementById('n-particles'),
            dt: document.getElementById('dt'),
            gravity: document.getElementById('gravity'),
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
