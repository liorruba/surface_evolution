# REGOLIT web UI

A FastAPI service (`web/server.py`) with a single static page (`web/static/index.html`). The page
posts a parameter set, the server runs the model in `runs/web/<id>/` and renders maps,
cross-sections, histograms and an animation on demand. Results can be shared with the page link
(`?run=<id>`) and downloaded as a zip of the raw output.

## Run locally

```bash
make                                   # the model binary
python3 -m venv venv && venv/bin/pip install -r web/requirements.txt
venv/bin/uvicorn web.server:app --reload --port 8010
# open http://127.0.0.1:8010
```

## Deploy on the droplet

The layout is the one used for moon.liorruba.com: uvicorn bound to localhost as a systemd service,
nginx in front, TLS from Let's Encrypt, HTTP basic authentication inside the app.

1. DNS: add an `A` record `regolit.liorruba.com` pointing at the droplet's IP (the moon
   application's droplet is 192.241.128.158).
2. On the droplet, as root:

   ```bash
   git clone https://github.com/liorruba/surface_evolution.git /opt/regolit
   DOMAIN=regolit.liorruba.com PORT=8010 WEB_USER=lior WEB_PASSWORD='choose-a-password' \
     CERTBOT_EMAIL=you@example.com bash /opt/regolit/deploy/setup_droplet.sh
   ```

   The script installs `g++`, `make`, `python3-venv` and `nginx`, builds the model, creates a
   virtual environment, writes `/etc/regolit-web.env`, installs and starts the `regolit-web`
   systemd service, enables the nginx site and, when `CERTBOT_EMAIL` is given, obtains the
   certificate. Pick another `PORT` if 8010 is taken by the moon service (`ss -ltnp | grep 8010`).
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
would exceed 500 x 500 cells, 50 output steps or about 200,000 impacts. At most two runs execute at
a time (others wait, up to eight), each with a 120 s timeout. The full layer stacks are not written
for web runs. Old runs are deleted once more than 200 are stored. All of these are environment
variables or constants at the top of `server.py`.
