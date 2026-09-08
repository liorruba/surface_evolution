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

The run setup offers, besides the raw model parameters:

- **Target body** presets (Moon, Mercury, Mars, Ceres, Vesta, custom) that set gravity, impact
  velocity and the target properties (density, effective strength `Ybar`, strength constant `K2`)
  after Holsapple (1993) as tabulated by Williams et al. (2014): lunar regolith for the airless
  bodies, dry soil for Mars. Editing one of those fields switches the preset to custom. The angle of
  repose is never changed by a preset (default 35 deg).
- **Production functions** for the impactor flux. The model samples a single power law
  `N(>d) = c d^-b` of impactor diameters:
  - *Williams et al. (2014)* (default): the annual flux of terrestrial fireballs (Brown et al. 2002,
    `log10 N(>E) = 0.5677 - 0.90 log10 E`) converted to impactor diameters with the impactor density
    and the body's mean impact velocity and scaled by the body/Earth flux ratio (Moon 0.725,
    Mars 1.885). It is an impactor distribution already, with cumulative slope 2.7.
  - *Neukum et al. (2001)*, lunar crater production function, and *Daubar et al. (2013)*, the
    present-day martian rate (1.65e-6 craters km^-2 yr^-1 with D >= 3.9 m, differential slope
    -2.45, i.e. cumulative slope 1.45). Crater functions are converted to impactor sizes through the
    model's own crater scaling and fitted by the power law over the run's impactor size range; the
    fitted `c`, `b` are shown.
  - *Power law (manual)*: type `c` and `b` yourself.
  All constants live in `python/regolit/scaling.py`.
- An **automatic basement thickness**: the initial layers plus three times the depth of the
  largest crater expected in the run, with a checkbox to override it.
- A **secondary craters** panel that unfolds when secondaries are switched on: minimum landing
  speed, largest secondary / primary diameter, fragment size–speed exponent, fragment size-distribution
  slope, depth/diameter, the per-primary budget, and the distant-primaries switch with its range
  (secondaries are ejecta fragments; see the main README).
- A live **estimate** line: expected number of impacts, largest impactor and crater, smallest crater,
  and the basement suggestion, recomputed as you edit (`POST /api/estimate`).
- A **Tests** tab on the setup page with predefined scenarios (`TESTS` in `server.py`); the first
  forms a single 1 km crater with its secondaries on a 10 km domain. A test runs like any other run
  (its card and title say so), and "Edit & run again" keeps the test crater (a notice says so in the
  setup form; Reset defaults returns to a normal run).
- **Saved settings**: the bar at the bottom of the setup page stores the parameters, layers and
  presets on the server under a name (`/api/settings`, one JSON file per name in
  `REGOLIT_WEB_SETTINGS`, by default the `settings` directory next to the runs directory). Pick a
  name from the combo to load it; saving under an existing name replaces it. A run's Setup tab also
  records its presets, and "Edit & run again" restores them.

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

Add `RUNS_DIR=/Data/liorr/regolit/runs MAX_DISK_GB=2000` to keep the runs on the data disk (this is
the live configuration on Haworth). This installs two systemd user units, `regolit-web` and `regolit-tunnel`, enables them at boot
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

## How a run executes

`POST /api/runs` validates the request, creates the run directory with a provisional
`summary.json` (status `queued`) and returns at once; the page switches to the run and polls
`GET /api/runs/{id}`, which adds the model's progress (parsed from its log) while the status is
`queued` or `running`. Up to `REGOLIT_WEB_CONCURRENCY` runs execute at a time; each one is a
detached worker process (`python -m web.worker <id>`, its own session) that runs the model, compacts
the layer stacks and writes the final summary itself. The worker therefore survives a reload or
restart of the server, a closed browser tab and any proxy timeout. On startup the server re-queues
runs that were still `queued`, leaves `running` workers alone, and marks runs whose worker died
as `failed`. Maps and downloads answer 409 until the run is `done`; a running run cannot be
deleted.

The worker is started with the standard library's `subprocess` rather than the event loop's
subprocess support: under uvloop a child inherits every inheritable descriptor, including uvicorn's
listening socket, and a long-running worker would then keep the port busy across server restarts.
Should that ever happen with an old worker, `LOCAL_PORT=8031 bash deploy/mediator/setup_mediator.sh`
moves uvicorn to another local port while the tunnel keeps publishing the same droplet port.

## Limits

Requests are validated against the ranges in `PARAMETERS` in `server.py`. A run is refused if it
would exceed 8000 x 8000 cells, 500 output steps, 500 GB of output maps or about twenty million
impacts. The region width itself is bounded only through the cell count (up to 1000 km at
1000 m/pixel). These are sized for the compute machine (16 cores, 125 GB): measured with 8 threads
over 100 Ma, 2000 x 2000 cells run in 18 s and need 0.8 GB, 4000 x 4000 cells run in 3 min and
need 3.3 GB, so 8000 x 8000 cells need about 15 GB and four of them fit at a time. Four runs execute
at a time (others wait, up to twelve), each with a 12 h timeout, and the slope relaxation of each run uses
`OMP_NUM_THREADS` threads (default: physical cores divided by the concurrency). The full layer
stacks are stored compactly for the cross-sections. Old runs are deleted once more than 500 are
stored or the runs directory exceeds 200 GB. All of these are environment variables
(`/etc/regolit-web.env` on the droplet, `~/.config/regolit/web.env` in mediator mode) or constants
at the top of `server.py`.
