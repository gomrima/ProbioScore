# ProbioScore I45 portable container.
#
# Build:
#   docker build -t probioscore-i45:1.0 .
#
# Run with a V4 outputs directory and an output directory mounted:
#   docker run --rm \
#     -v /path/to/v4_outputs:/data/v4_outputs:ro \
#     -v /path/to/out:/data/out \
#     probioscore-i45:1.0 \
#     probioscore --v4-tsv-dir /data/v4_outputs --out-dir /data/out \
#                     --mode prospective_frozen --skip-visuals

FROM python:3.12-slim

LABEL org.opencontainers.image.title="ProbioScore I45"
LABEL org.opencontainers.image.description="Portable AHP plus FCE probiogenomic triage calculator (I31 to I45 consolidated)"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.version="1.0"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /opt/probioscore

# Copy the calculator source tree. The .dockerignore controls exclusions.
COPY . /opt/probioscore

# Install the package and its Python dependencies.
RUN python -m pip install --upgrade pip \
 && python -m pip install .

# Create a non-root user for container execution.
RUN useradd --create-home --shell /bin/bash probioscore \
 && mkdir -p /data/v4_outputs /data/out \
 && chown -R probioscore:probioscore /data

USER probioscore
WORKDIR /home/probioscore

# Build-time smoke check on the console entry point.
RUN probioscore --help >/dev/null

ENTRYPOINT []
CMD ["probioscore", "--help"]
