# Multi-stage Dockerfile for WisecondorX using Pixi
# Build stage
FROM ghcr.io/prefix-dev/pixi:latest AS build

# Set working directory
WORKDIR /app

# Copy project files
COPY . .

# Install runtime dependencies into the default environment using the lock file
RUN pixi install --locked -e default

# Runtime stage: Use a clean Ubuntu base for minimal image size
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Set up environment path
ENV PATH="/app/.pixi/envs/default/bin:$PATH"
WORKDIR /app

# Copy the environment from the build stage
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
# Copy the source code as it is installed in editable mode in the environment
COPY --from=build /app/src /app/src
COPY --from=build /app/pyproject.toml /app/pyproject.toml
