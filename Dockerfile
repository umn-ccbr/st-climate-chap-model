# Dockerfile — st-climate-chap-model
#
# Builds on the chapkit-r-inla base image (R 4.5 + INLA + spatial stack),
# installs CARBayesST and yaml, then installs the Python chapkit service layer.
#
# Build:
#   docker build -t malaria-spatiotemporal .
#
# Run directly:
#   docker run --rm -p 9090:8000 malaria-spatiotemporal

FROM ghcr.io/umn-ccbr/st-climate-chap-model-base:latest

WORKDIR /work

COPY pyproject.toml               ./pyproject.toml
COPY uv.lock                      ./uv.lock

# ── Install Python deps (chapkit + spatial utils) ───────────────────────────
RUN uv pip install chapkit geopandas

# ── Copy model source files ────────────────────────────────────────────────
COPY utils.py                     ./utils.py
COPY train.r                      ./train.r
COPY predict.r                    ./predict.r
COPY model_helpers.R              ./model_helpers.R
COPY main.py                      ./main.py
COPY scripts/                     ./scripts/

LABEL org.opencontainers.image.source=https://github.com/umn-ccbr/st-climate-chap-model

CMD ["fastapi", "run", "main.py", "--host", "0.0.0.0", "--port", "8000"]
