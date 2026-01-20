#!/usr/bin/env bash

# Usage: ./run_web_server.sh [PORT] [DOCROOT]
# Defaults: PORT=8000, DOCROOT=.
# - Kills any existing process listening on PORT
# - Then starts `python3 -m http.server PORT` in DOCROOT




set -euo pipefail

PORT="${1:-8000}"
DOCROOT="${2:-.}"
echo "running serever on http://localhost:${PORT}"
echo "serving from ${DOCROOT}"

# Kill anything already on PORT (prefer lsof; fall back to fuser)
kill_listen(){
  local pids
  pids="$(lsof -tiTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null || true)"
  if [[ -n "${pids}" ]]; then
    echo "Killing existing server(s) on port ${PORT}: ${pids}" >&2
    kill ${pids} 2>/dev/null || true
    return 0
  fi
  if command -v fuser >/dev/null 2>&1; then
    # fuser prints PIDs to stderr
    if fuser -k "${PORT}/tcp" >/dev/null 2>&1; then
      echo "Killed via fuser on port ${PORT}" >&2
      return 0
    fi
  fi
  return 1
}

for _ in {1..5}; do
  kill_listen || true
  sleep 0.2
  if ! lsof -tiTCP:"${PORT}" -sTCP:LISTEN >/dev/null 2>&1; then
    break
  fi
done

if lsof -tiTCP:"${PORT}" -sTCP:LISTEN >/dev/null 2>&1; then
  echo "Port ${PORT} still busy after kill attempts. Aborting." >&2
  exit 1
fi

# Get the directory of this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Create a python script to run the server with no-cache headers
cat <<EOF > "$SCRIPT_DIR/server.py"
import http.server
import socketserver
import sys

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 8000

class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

with socketserver.TCPServer(("", PORT), NoCacheHTTPRequestHandler) as httpd:
    print(f"Serving at port {PORT} with no-cache headers")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
EOF

echo "Starting custom no-cache server on port ${PORT} in ${DOCROOT}" >&2
cd "${DOCROOT}"
exec python3 "$SCRIPT_DIR/server.py" "${PORT}"
