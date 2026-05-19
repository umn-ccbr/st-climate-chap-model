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

FROM ghcr.io/dhis2-chap/chapkit-r-inla:latest

WORKDIR /work

# ── Install R packages not already in the base image ──────────────────────
RUN R -e "install.packages(c('CARBayesST', 'yaml'), repos = 'https://cloud.r-project.org')"

# ── Install Python deps (chapkit + spatial utils) ───────────────────────────
RUN uv pip install chapkit geopandas

# ── Copy model source files ────────────────────────────────────────────────
COPY utils.py                     ./utils.py
COPY train.r                      ./train.r
COPY predict.r                    ./predict.r
COPY model_helpers.R              ./model_helpers.R
COPY main.py                      ./main.py
COPY scripts/                     ./scripts/

CMD ["fastapi", "run", "main.py", "--host", "0.0.0.0", "--port", "8000"]
