#!/usr/bin/env bash
# Install or update the REGOLIT web UI on an Ubuntu droplet (tested layout: Ubuntu 22.04/24.04 with
# nginx, the same arrangement as moon.liorruba.com). Safe to re-run: it pulls the latest code,
# rebuilds the model, reinstalls the Python environment and restarts the service.
#
# Usage (as root or with sudo):
#   DOMAIN=regolit.liorruba.com PORT=8010 WEB_USER=lior WEB_PASSWORD='secret' bash deploy/setup_droplet.sh
#
# Variables (all optional):
#   DOMAIN        host name served by nginx            (default regolit.liorruba.com)
#   PORT          local port of the uvicorn service     (default 8010)
#   APP_DIR       installation directory                (default /opt/regolit)
#   REPO_URL      git repository to clone               (default https://github.com/liorruba/surface_evolution.git)
#   BRANCH        branch to deploy                      (default master)
#   SERVICE_USER  system user running the service       (default regolit, created if missing)
#   WEB_USER / WEB_PASSWORD  HTTP basic-auth credentials; leave empty for an open site
#   CERTBOT_EMAIL if set, obtains a Let's Encrypt certificate non-interactively
#   SKIP_GIT=1    do not clone or pull: deploy the code already present in APP_DIR (e.g. after rsync)
#   AUTO_UPDATE   1 (default) installs a systemd timer that redeploys when the GitHub branch moves; 0 disables it
set -euo pipefail

DOMAIN="${DOMAIN:-regolit.liorruba.com}"
PORT="${PORT:-8010}"
APP_DIR="${APP_DIR:-/opt/regolit}"
REPO_URL="${REPO_URL:-https://github.com/liorruba/surface_evolution.git}"
BRANCH="${BRANCH:-master}"
SERVICE_USER="${SERVICE_USER:-regolit}"
WEB_USER="${WEB_USER:-}"
WEB_PASSWORD="${WEB_PASSWORD:-}"
CERTBOT_EMAIL="${CERTBOT_EMAIL:-}"
SKIP_GIT="${SKIP_GIT:-}"
AUTO_UPDATE="${AUTO_UPDATE:-1}"

if [[ $EUID -ne 0 ]]; then echo "run as root (sudo)"; exit 1; fi

echo "==> packages"
apt-get update -qq
DEBIAN_FRONTEND=noninteractive apt-get install -y -qq g++ make git python3 python3-venv nginx >/dev/null
if [[ -n "$CERTBOT_EMAIL" ]]; then DEBIAN_FRONTEND=noninteractive apt-get install -y -qq certbot python3-certbot-nginx >/dev/null; fi

echo "==> service user and code"
id -u "$SERVICE_USER" >/dev/null 2>&1 || useradd --system --create-home --home-dir "/home/$SERVICE_USER" --shell /usr/sbin/nologin "$SERVICE_USER"
if [[ -n "$SKIP_GIT" ]]; then
  [[ -f "$APP_DIR/makefile" ]] || { echo "SKIP_GIT is set but $APP_DIR does not contain the code"; exit 1; }
elif [[ -d "$APP_DIR/.git" ]]; then
  # The tree may be owned by another user from an earlier install; tell git it is trusted.
  git -c safe.directory="$APP_DIR" -C "$APP_DIR" fetch --quiet origin
  git -c safe.directory="$APP_DIR" -C "$APP_DIR" checkout --quiet "$BRANCH"
  git -c safe.directory="$APP_DIR" -C "$APP_DIR" reset --quiet --hard "origin/$BRANCH"
else
  git clone --quiet --branch "$BRANCH" "$REPO_URL" "$APP_DIR"
fi

echo "==> model binary"
make -C "$APP_DIR" -s
"$APP_DIR/build/apps/regolit_main.run" --help >/dev/null 2>&1 || true

echo "==> python environment"
[[ -d "$APP_DIR/venv" ]] || python3 -m venv "$APP_DIR/venv"
"$APP_DIR/venv/bin/pip" install --quiet --upgrade pip
"$APP_DIR/venv/bin/pip" install --quiet -r "$APP_DIR/web/requirements.txt"
# The code and the virtual environment stay owned by root; the service only needs to write runs/.
chown -R root:root "$APP_DIR"
mkdir -p "$APP_DIR/runs/web"
chown -R "$SERVICE_USER:$SERVICE_USER" "$APP_DIR/runs"

echo "==> configuration (/etc/regolit-web.env)"
{
  echo "REGOLIT_BINARY=$APP_DIR/build/apps/regolit_main.run"
  echo "REGOLIT_WEB_RUNS=$APP_DIR/runs/web"
  echo "REGOLIT_WEB_CONCURRENCY=2"
  echo "REGOLIT_WEB_TIMEOUT=120"
  echo "REGOLIT_WEB_MAX_RUNS=200"
  if [[ -n "$WEB_USER" ]]; then
    echo "REGOLIT_WEB_USER=$WEB_USER"
    echo "REGOLIT_WEB_PASSWORD=$WEB_PASSWORD"
  fi
} > /etc/regolit-web.env
chmod 600 /etc/regolit-web.env

echo "==> systemd service"
sed -e "s|__USER__|$SERVICE_USER|g" -e "s|__APP_DIR__|$APP_DIR|g" -e "s|__PORT__|$PORT|g" \
  "$APP_DIR/deploy/regolit-web.service" > /etc/systemd/system/regolit-web.service
systemctl daemon-reload
systemctl enable --quiet regolit-web
systemctl restart regolit-web
sleep 2
systemctl --no-pager --lines=5 status regolit-web || true

echo "==> auto-update timer"
if [[ "$AUTO_UPDATE" == "1" && -z "$SKIP_GIT" ]]; then
  sed -e "s|__USER__|$SERVICE_USER|g" -e "s|__APP_DIR__|$APP_DIR|g" -e "s|__BRANCH__|$BRANCH|g" \
    "$APP_DIR/deploy/regolit-update.service" > /etc/systemd/system/regolit-update.service
  cp "$APP_DIR/deploy/regolit-update.timer" /etc/systemd/system/regolit-update.timer
  systemctl daemon-reload
  systemctl enable --quiet --now regolit-update.timer
  echo "pushes to $BRANCH are deployed within about two minutes (journalctl -u regolit-update)"
else
  systemctl disable --quiet --now regolit-update.timer 2>/dev/null || true
fi

echo "==> nginx site"
sed -e "s|__DOMAIN__|$DOMAIN|g" -e "s|__PORT__|$PORT|g" "$APP_DIR/deploy/nginx-regolit.conf" > /etc/nginx/sites-available/regolit
ln -sf /etc/nginx/sites-available/regolit /etc/nginx/sites-enabled/regolit
nginx -t
systemctl reload nginx

if [[ -n "$CERTBOT_EMAIL" ]]; then
  echo "==> TLS certificate"
  certbot --nginx --non-interactive --agree-tos --redirect -m "$CERTBOT_EMAIL" -d "$DOMAIN" || echo "certbot failed; run it manually once DNS for $DOMAIN points here"
fi

echo
echo "Done. Check: curl -s http://127.0.0.1:$PORT/health ; then open https://$DOMAIN"
echo "Logs: journalctl -u regolit-web -f"
