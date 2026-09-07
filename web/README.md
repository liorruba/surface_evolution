# REGOLIT web UI

A FastAPI service (`web/server.py`) with a single static page (`web/static/index.html`), laid out
like the lunar thermal model dashboard: a rail with the stored runs and a "New run" composer, and a
stage with the selected run. The page posts a parameter set, the server runs the model in
`runs/web/<id>/` and renders on demand:

- maps of every saved step: shaded relief (hillshade of the elevation, sun direction selectable),
  elevation, and surface or depth-integrated composition, with a step slider and animation;
- a layered subsurface cross-section along a west-east or south-north line that the user places by
  clicking the map or with a slider, drawn from the full-resolution layer stacks of the final state
  (`isPrintSubsurface 2`, stored compactly as `subsurface_XX.npz`), with the surface of the selected
  step overlaid;
- crater and impactor histograms, the parameter set, and downloads (zip of the raw output, crater
  list, animation, log).

Results can be shared with the page link (`?run=<id>`).

## Run locally

```bash
make                                   # the model binary
python3 -m venv venv && venv/bin/pip install -r web/requirements.txt
venv/bin/uvicorn web.server:app --reload --port 8010
# open http://127.0.0.1:8010
```

## Mediator mode (recommended, as for moon.liorruba.com)

The model and the web server run on the compute machine straight from the working copy, and the
droplet only forwards: a reverse SSH tunnel from the compute machine publishes the server on the
droplet's `127.0.0.1:8030`, where nginx proxies `regolit.liorruba.com` to it. Python changes are
live immediately (uvicorn `--reload`); after a C++ change run `systemctl --user restart regolit-web`,
which rebuilds the binary. Nothing is pulled from GitHub.

On the compute machine, as your user (the machine's SSH key must be authorized on the droplet):

```bash
WEB_USER=regolith WEB_PASSWORD='moon' PORT=8030 DROPLET=root@192.241.128.158 bash deploy/mediator/setup_mediator.sh
```

This installs two systemd user units, `regolit-web` and `regolit-tunnel`, enables them at boot
(`loginctl enable-linger` keeps them running without a login session) and stores the credentials in
`~/.config/regolit/web.env`. On the droplet, as root, install only the nginx site:

```bash
MODE=proxy DOMAIN=regolit.liorruba.com PORT=8030 CERTBOT_EMAIL=you@example.com bash /opt/regolit/deploy/setup_droplet.sh
```

Logs: `journalctl --user -u regolit-web -f` and `journalctl --user -u regolit-tunnel -f` on the
compute machine.

## Deploy on the droplet itself (alternative)

The layout is the one used for moon.liorruba.com: uvicorn bound to localhost as a systemd service,
nginx in front, TLS from Let's Encrypt, HTTP basic authentication inside the app.

1. DNS: add an `A` record `regolit.liorruba.com` pointing at the droplet's IP (the moon
   application's droplet is 192.241.128.158).
2. On the droplet, as root:

   ```bash
   git clone https://github.com/liorruba/surface_evolution.git /opt/regolit
   DOMAIN=regolit.liorruba.com PORT=8030 WEB_USER=regolith WEB_PASSWORD='moon' \
     CERTBOT_EMAIL=you@example.com bash /opt/regolit/deploy/setup_droplet.sh
   ```

   The script installs `g++`, `make`, `python3-venv` and `nginx`, builds the model, creates a
   virtual environment, writes `/etc/regolit-web.env`, installs and starts the `regolit-web`
   systemd service, enables the nginx site and, when `CERTBOT_EMAIL` is given, obtains the
   certificate. Ports 8010, 8011, 8020 and 8021 on the droplet belong to the moon app's tunnels, so use 8030 or another free port (`ss -ltnp`).
   The firewall must allow 80 and 443 (`ufw allow 'Nginx Full'`), which is already the case if
   nginx serves the moon site.
3. Updating: the installer enables a systemd timer (`regolit-update.timer`) that checks GitHub
   every two minutes and, when `master` has new commits, pulls, rebuilds, reinstalls the Python
   requirements and restarts the service. A push therefore goes live within a few minutes; watch it
   with `journalctl -u regolit-update -f`. To update immediately run
   `bash /opt/regolit/deploy/update_droplet.sh`; to disable automatic updates run
   `systemctl disable --now regolit-update.timer` (or install with `AUTO_UPDATE=0`). Because the
   updater resets the checkout to `origin/master`, do not edit files in `/opt/regolit` on the droplet.

Logs: `journalctl -u regolit-web -f`. Health check: `curl http://127.0.0.1:8010/health`.

## Limits

Requests are validated against the ranges in `PARAMETERS` in `server.py`, and a run is refused if it
would exceed 500 x 500 cells, 500 output steps, 800 MB of output maps or about 200,000 impacts. At most two runs execute at
a time (others wait, up to eight), each with a 120 s timeout. The full layer stacks are not written
for web runs. Old runs are deleted once more than 200 are stored. All of these are environment
variables or constants at the top of `server.py`.
