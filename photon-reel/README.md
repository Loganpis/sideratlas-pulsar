# SiderAtlas NO.002 — outgoing gamma-ray photon Reel

The original 0.1–100 GeV physics is preserved. This renderer replaces v1's
emission-location volume with photons propagating from actual source cells
in their calculated laboratory directions. A fixed observer indicator and
the full original delayed pulse provide the lighthouse connection.

This is a Vela-informed visualization of an approximate calculation, not a
validated Vela reconstruction or a literal gamma-ray camera image.

## Delivered files

- [32-second Reel](https://github.com/Loganpis/sideratlas-pulsar/releases/download/photon-v1.0.0/sideratlas-002-pulsar-reel.mp4): 1080×1920, 30 fps
- `release/photon-v2/instagram-description.text`: upload caption, credits and false-color key
- `release/photon-v2/MODEL-NOTES.md`: geometry, timing, physical and visualization assumptions
- `outputs/sideratlas-pulsar-photons/scene/photons.npz`: ready-to-render calculated photon scene
- `outputs/sideratlas-pulsar-photons/media/observer-frame-flux.csv`: exact frame-averaged observer flux
- Calculation, direction/timing and encoded-video verification JSON files
- `pulse-timing-check.png` and `contact-sheet.jpg`: timing comparison and eight loop phases

## Rerender from saved results

Python 3.12+, NumPy, SciPy, Pillow and FFmpeg with libx264/zscale are required.
The focused tests additionally use pytest. From this `photon-reel` directory (or the downloaded bundle root):

```sh
python -m pip install numpy scipy Pillow pytest
PYTHONPATH=src python -m pytest -q tests/test_photon_scene.py tests/test_emission_quadrature.py tests/test_outer_transport.py tests/test_radiation.py
python render_photon_reel.py --output outputs/rerender
python verify_reel.py outputs/rerender
```

No field solve or particle integration is needed. The renderer verifies
directions, propagation and timing before exporting. It encodes one 8-second
cycle and repeats the encoded cycle four times without another lossy encode.

`prepare_photon_scene.py` is supplied for audit/repreparation in the original
project, where the six original trajectory-cache chunks are available. It
reuses those trajectories and the original spectral integration; the scene
in this package removes any need to retrieve or rerun that large cache.
The original protocol/profile/result receipts are included. Preparation
refuses to overwrite an existing scene; rerendering may overwrite its chosen
output folder.

The original project and full field-solver provenance are at
https://github.com/Loganpis/sideratlas-pulsar

Code is GPL-2.0; see LICENSE, NOTICE and bundled font licenses. The original
media and caption retain the original project's media copyright notice.

The exact standalone bundle, video and checksums are available in the [photon release](https://github.com/Loganpis/sideratlas-pulsar/releases/tag/photon-v1.0.0)
