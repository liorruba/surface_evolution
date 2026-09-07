#!/usr/bin/env bash
# Redeploy REGOLIT when the GitHub branch has moved. Run periodically by regolit-update.timer
# (installed by setup_droplet.sh), or by hand: bash /opt/regolit/deploy/update_droplet.sh
# Exits quietly when nothing changed. Local edits in APP_DIR are discarded on update.
set -euo pipefail

APP_DIR="${APP_DIR:-/opt/regolit}"
BRANCH="${BRANCH:-master}"
SERVICE_USER="${SERVICE_USER:-regolit}"

git -c safe.directory="$APP_DIR" -C "$APP_DIR" fetch --quiet origin
current=$(git -c safe.directory="$APP_DIR" -C "$APP_DIR" rev-parse HEAD)
latest=$(git -c safe.directory="$APP_DIR" -C "$APP_DIR" rev-parse "origin/$BRANCH")
if [[ "$current" == "$latest" ]]; then
  exit 0
fi

echo "updating ${current:0:7} -> ${latest:0:7}"
git -c safe.directory="$APP_DIR" -C "$APP_DIR" reset --quiet --hard "origin/$BRANCH"

if ! make -C "$APP_DIR" -s; then
  echo "build failed; the running service was left untouched (old binary and old Python code may now be mixed)"
  exit 1
fi
"$APP_DIR/venv/bin/pip" install --quiet -r "$APP_DIR/web/requirements.txt"

chown -R root:root "$APP_DIR"
mkdir -p "$APP_DIR/runs/web"
chown -R "$SERVICE_USER:$SERVICE_USER" "$APP_DIR/runs"

systemctl restart regolit-web
echo "deployed ${latest:0:7}: $(git -c safe.directory="$APP_DIR" -C "$APP_DIR" log -1 --format=%s)"
