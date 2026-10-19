# @sideratlas · NO.002 · THE PULSAR

Calculated **0.1–100 GeV gamma-ray photons** sweep through space and periodically reach a fixed observer. The tiny central point marks the neutron star. The blue radiation pattern and white direction markers represent the same calculated photon population

- [32-second Reel](https://github.com/Loganpis/sideratlas-pulsar/releases/download/photon-v1.1.0/sideratlas-002-pulsar-reel.mp4) · 1080×1920 · 30 fps
- [Instagram caption](instagram-description.text)
- [Model, equations, assumptions and study credits](docs/MODEL-NOTES.md)
- [Complete downloadable bundle](https://github.com/Loganpis/sideratlas-pulsar/releases/tag/photon-v1.1.0)

![Calculated gamma-ray radiation around the pulsar](media/preview-00.png)

This is a **Vela-informed visualization of an approximate calculation**. Its pulse and spectrum differ from Vela's measurements. The displayed rays follow calculated emission positions and laboratory directions; the observer pulse includes the original photon travel-time delays. Gamma-ray colors are explicitly false colors. No artificial cones or luminous magnetic-field guides are added

## Render the saved calculation

Python 3.12 and FFmpeg with libx264/zscale are required

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python render_photon_reel.py
.venv/bin/python verify_reel.py media
```

The supplied `scene/photons.npz` contains the calculated source positions, photon directions, spectral weights and original observer profile. Rendering requires **no field solve or particle simulation**. One physical 89.371 ms rotation is displayed over eight seconds and repeated four times. The camera stays external and fixed; the title is horizontally centered above the pulsar

## Find the relevant files

| Location | Purpose |
|---|---|
| `render_photon_reel.py`, `reel_export.py` | Photon propagation display, labels and video export |
| `scene/` | Ready photon scene and original calculation checks |
| `src/pulsar/` | Photon transport, curvature radiation and supporting numerical routines |
| `media/` | Poster, timing diagram and measured export checks; rendered videos are release downloads |
| `instagram-description.text` | Caption ready to copy |
| `docs/MODEL-NOTES.md` | Current scientific and visualization assumptions, equations and citations |
| [`physics/`](physics/README.md) | The force-free field, particle calculation, selected snapshot, source patches and provenance used for this animation |
| `tests/`, `physics/tests/` | Relevant software and numerical checks |

```sh
.venv/bin/python -m pytest -q tests physics/tests
```

Directions, light-speed propagation, original arrival delays and loop periodicity are checked before export. Separate verification inspects the encoded video and repeated cycles. These checks establish implementation properties; they do not establish observational agreement with Vela

The code is GPL-2.0. Font, upstream software and media terms are identified in [NOTICE](NOTICE). The full study list and source roles appear in the [model notes](docs/MODEL-NOTES.md)
