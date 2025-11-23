// Global Verbosity Level
// Global Verbosity Level
// window.VERBOSITY_LEVEL removed in favor of Logger instance properties

class Logger {
    static NONE = 0;
    static ERROR = 1;
    static WARN = 2;
    static INFO = 3;
    static DEBUG = 4;

    constructor() {
        this.domElement = null; // Will be set by GUI
        this.consoleVerbosity = 4; // Default: DEBUG
        this.uiVerbosity = 4;      // Default: DEBUG
    }

    setContainer(element) {
        this.domElement = element;
    }

    setConsoleVerbosity(level) {
        this.consoleVerbosity = level;
    }

    setUIVerbosity(level) {
        this.uiVerbosity = level;
    }

    shouldLog(level) {
        return level <= Math.max(this.consoleVerbosity, this.uiVerbosity);
    }

    clear() {
        if (this.domElement) {
            this.domElement.textContent = '';
        }
    }

    log(message, level = Logger.INFO) {
        // Caller info (simple stack trace parsing)
        let caller = "";
        try {
            throw new Error();
        } catch (e) {
            // Stack: Error \n at Logger.log \n at Logger.info \n at Caller
            const lines = e.stack.split('\n');
            if (lines.length >= 4) {
                const match = lines[3].match(/at\s+([^\s]+)/);
                if (match) caller = match[1];
            }
        }

        // Map level to string for display
        let levelStr = 'INFO';
        let color = '#ffffff';

        if (level === Logger.ERROR) { levelStr = 'ERROR'; color = '#ff5555'; }
        else if (level === Logger.WARN) { levelStr = 'WARN'; color = '#ffaa00'; }
        else if (level === Logger.DEBUG) { levelStr = 'DEBUG'; color = '#55ff55'; }

        const formattedMsg = caller ? `[${levelStr}] [${caller}] ${message}` : `[${levelStr}] ${message}`;

        // Console
        if (level <= this.consoleVerbosity) {
            if (level === Logger.ERROR) console.error(formattedMsg);
            else if (level === Logger.WARN) console.warn(formattedMsg);
            else console.log(formattedMsg);
        }

        // DOM
        if (this.domElement && level <= this.uiVerbosity) {
            const line = document.createElement('div');
            line.textContent = formattedMsg;
            line.style.color = color;
            line.style.fontFamily = 'monospace';
            line.style.whiteSpace = 'pre-wrap';
            line.style.borderBottom = '1px solid #333';
            line.style.padding = '2px 0';

            this.domElement.appendChild(line);
            this.domElement.scrollTop = this.domElement.scrollHeight;
        }
    }

    debug(msg) { this.log(msg, Logger.DEBUG); }
    info(msg) { this.log(msg, Logger.INFO); }
    warn(msg) { this.log(msg, Logger.WARN); }
    error(msg) { this.log(msg, Logger.ERROR); }
}

// Global instance
window.logger = new Logger();
