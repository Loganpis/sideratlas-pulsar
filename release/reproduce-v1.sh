#!/bin/sh
# Run from the repository root, in the Python environment described by README.
# Choose fresh output directories: completed scientific runs are not overwritten.
set -eu
python run_emission_quadrature.py recorded-inputs/snapshot-001020.npz \
  --output outputs/reproduced-v1-calculation \
  --azimuths 128 --rings 6 --dt .002 --widths 2 --nphase 1000 \
  --cap-boundary outgoing_sphere
python prepare_reel_scene.py outputs/reproduced-v1-calculation \
  --output outputs/reproduced-v1-scene --grid 176 --angular-display-samples 9
python render_reel.py outputs/reproduced-v1-scene/scene.npz \
  --output outputs/reproduced-v1-media --width 1080 --fps 30 --repeats 4
python verify_reel.py outputs/reproduced-v1-media
