FROM ghcr.io/prefix-dev/pixi:latest

WORKDIR /app

# Copy the environment lock / config
COPY pixi.lock pyproject.toml ./

# Copy the rest of the project required by pypi builds
COPY src/ ./src/
COPY README.md LICENSE.md ./

# Install dependencies and the project
RUN pixi install --locked

# Set the entrypoint to run the wisecondorx CLI
ENTRYPOINT ["pixi", "run", "wisecondorx"]