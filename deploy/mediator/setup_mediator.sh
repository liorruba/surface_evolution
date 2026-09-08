#!/usr/bin/env bash
# Set up "mediator mode" on the compute machine (the arrangement used for moon.liorruba.com):
# the web UI runs here from this working copy and a reverse SSH tunnel publishes it on the droplet,
# whose nginx forwards the public host name to the tunnel. Run as your normal user, from anywhere:
#
#   WEB_USER=regolith WEB_PASSWORD='moon' PORT=8030 DROPLET=root@192.241.128.158 bash deploy/mediator/setup_mediator.sh
#
# RUNS_DIR chooses where the runs are stored (default runs/web in the repository; on Haworth
# /Data/liorr/regolit/runs on the 20 TB data disk) and MAX_DISK_GB the disk budget for them.
# Re-running updates the units and the environment (credentials are kept if not given).
# Requirements: python3-venv, make, g++, an SSH key on this machine authorized on the droplet.
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PORT="${PORT:-8030}"
DROPLET="${DROPLET:-root@192.241.128.158}"
WEB_USER="${WEB_USER:-}"
WEB_PASSWORD="${WEB_PASSWORD:-}"
UNIT_DIR="$HOME/.config/systemd/user"
ENV_FILE="$HOME/.config/regolit/web.env"
RUNS_DIR="${RUNS_DIR:-$APP_DIR/runs/web}"
MAX_DISK_GB="${MAX_DISK_GB:-200}"

echo "==> model binary and Python environment ($APP_DIR)"
make -s -C "$APP_DIR"
[[ -d "$APP_DIR/venv" ]] || python3 -m venv "$APP_DIR/venv"
"$APP_DIR/venv/bin/pip" install --quiet --upgrade pip
"$APP_DIR/venv/bin/pip" install --quiet -r "$APP_DIR/web/requirements.txt"
mkdir -p "$RUNS_DIR"

echo "==> environment ($ENV_FILE)"
mkdir -p "$(dirname "$ENV_FILE")"
if [[ -n "$WEB_USER" || ! -f "$ENV_FILE" ]]; then
  {
    echo "REGOLIT_BINARY=$APP_DIR/build/apps/regolit_main.run"
    echo "REGOLIT_WEB_RUNS=$RUNS_DIR"
    echo "REGOLIT_WEB_CONCURRENCY=4"
    echo "REGOLIT_WEB_TIMEOUT=43200"
    echo "REGOLIT_WEB_MAX_RUNS=500"
    echo "REGOLIT_WEB_MAX_DISK_GB=$MAX_DISK_GB"
    if [[ -n "$WEB_USER" ]]; then
      echo "REGOLIT_WEB_USER=$WEB_USER"
      echo "REGOLIT_WEB_PASSWORD=$WEB_PASSWORD"
    fi
  } > "$ENV_FILE"
  chmod 600 "$ENV_FILE"
fi

echo "==> user services"
mkdir -p "$UNIT_DIR"
sed -e "s|__APP_DIR__|$APP_DIR|g" -e "s|__PORT__|$PORT|g" "$APP_DIR/deploy/mediator/regolit-web.service" > "$UNIT_DIR/regolit-web.service"
sed -e "s|__PORT__|$PORT|g" -e "s|__DROPLET__|$DROPLET|g" "$APP_DIR/deploy/mediator/regolit-tunnel.service" > "$UNIT_DIR/regolit-tunnel.service"
systemctl --user daemon-reload
systemctl --user enable regolit-web.service regolit-tunnel.service
systemctl --user restart regolit-web.service regolit-tunnel.service

# Keep the user services running when no session is open:
if [[ "$(loginctl show-user "$USER" -p Linger --value 2>/dev/null)" != "yes" ]]; then
  loginctl enable-linger "$USER" 2>/dev/null || echo "NOTE: run 'sudo loginctl enable-linger $USER' once so the services survive logouts and reboots."
fi

sleep 3
echo "==> status"
systemctl --user --no-pager --lines=0 status regolit-web.service regolit-tunnel.service | grep -E 'regolit|Active' || true
curl -s -o /dev/null -w "local health: HTTP %{http_code}\n" "http://127.0.0.1:$PORT/health" || true
echo
echo "On the droplet, nginx must proxy the public host name to 127.0.0.1:$PORT (setup_droplet.sh MODE=proxy),"
echo "and this machine's key ($HOME/.ssh/id_ed25519.pub) must be in the droplet's authorized_keys for $DROPLET."
echo "Tunnel log: journalctl --user -u regolit-tunnel -f"
