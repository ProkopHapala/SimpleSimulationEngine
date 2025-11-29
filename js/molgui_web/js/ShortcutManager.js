export class ShortcutManager {
    constructor(editor) {
        console.log("[ShortcutManager] Constructor called");
        this.editor = editor;
        this.shortcuts = {};
        this.init();
    }

    init() {
        console.log("[ShortcutManager] Init called - attaching listener");
        window.addEventListener('keydown', this.onKeyDown.bind(this));

        // Define default shortcuts
        this.register('g', () => this.editor.toggleGizmo(), 'Toggle Gizmo');
        this.register('t', () => this.editor.setGizmoMode('translate'), 'Translate Mode');
        this.register('r', () => this.editor.setGizmoMode('rotate'), 'Rotate Mode');
        this.register('s', () => this.editor.setGizmoMode('scale'), 'Scale Mode');
        this.register('Escape', () => this.editor.clearSelection(), 'Clear Selection');
        this.register('Delete', () => this.editor.deleteSelection(), 'Delete Selection');
        this.register('Backspace', () => this.editor.deleteSelection(), 'Delete Selection');
        this.register('b', () => this.editor.recalculateBonds(), 'Recalculate Bonds');
        this.register('a', () => this.editor.addAtom(), 'Add Atom');
    }

    register(key, action, description) {
        this.shortcuts[key.toLowerCase()] = { action, description };
    }

    onKeyDown(e) {
        // Ignore if typing in an input field
        if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;

        const key = e.key.toLowerCase();

        // Emergency Debug
        console.log(`[ShortcutManager] RAW KEYDOWN: ${e.key}`);

        window.logger.debug(`[ShortcutManager] Key pressed: '${e.key}' (mapped to '${key}')`);

        if (this.shortcuts[key]) {
            window.logger.info(`[ShortcutManager] Executing action: ${this.shortcuts[key].description}`);
            this.shortcuts[key].action();
        } else {
            window.logger.debug(`[ShortcutManager] No shortcut found for '${key}'`);
        }
    }
}
