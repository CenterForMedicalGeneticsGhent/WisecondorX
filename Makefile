# ====================================================================================
# Variables
# ====================================================================================

# The name of your application's binary
BINARY_NAME=wisecondorx

# The main package of your application
MAIN_PACKAGE=.

# The directory where build artifacts will be stored
BUILD_DIR=./bin

# Go commands
GO=go
GOBUILD=$(GO) build
GOCLEAN=$(GO) clean
GOTOOLS=$(GO) tool
GOFMT=$(GO) fmt ./...
GOLINT=golangci-lint run
GOTIDY=$(GO) mod tidy
GOVET=$(GO) vet ./...

# List of platforms and architectures to build for
# Format: GOOS/GOARCH
PLATFORMS ?= linux/amd64 linux/arm64 darwin/amd64 darwin/arm64 windows/amd64 windows/arm64

# ====================================================================================
# Default Target
# ====================================================================================

# The default target executed when you run `make`
all: lint build

# ====================================================================================
# Build Targets
# ====================================================================================

# Build the application for the host system's OS and architecture
build-local:
	@echo "Building for local development..."
	@$(GOBUILD) -v -o $(BUILD_DIR)/$(BINARY_NAME) $(MAIN_PACKAGE)

# Build the application for all specified platforms
build:
	@echo "Building for all target platforms..."
	@$(foreach platform, $(PLATFORMS), $(call build_platform,$(platform)))

# Helper function to build for a specific platform
define build_platform
	$(eval GOOS_GOARCH = $(subst /, ,$(1)))
	$(eval GOOS = $(word 1,$(GOOS_GOARCH)))
	$(eval GOARCH = $(word 2,$(GOOS_GOARCH)))
	@echo "--> Building for $(GOOS)/$(GOARCH)"
	@GOOS=$(GOOS) GOARCH=$(GOARCH) $(GOBUILD) -v -o $(BUILD_DIR)/$(BINARY_NAME)-$(GOOS)-$(GOARCH) $(MAIN_PACKAGE)
endef


# ====================================================================================
# Linting and Formatting
# ====================================================================================

# Run the linter
lint:
	@echo "Running linter..."
	@$(GOLINT)

# Run go vet
vet:
	@echo "Running go vet..."
	@$(GOVET)

# Format the source code
fmt:
	@echo "Formatting source code..."
	@$(GOFMT)

# Tidy the go.mod file
tidy:
	@echo "Tidying go.mod file..."
	@$(GOTIDY)

# ====================================================================================
# Housekeeping Targets
# ====================================================================================

# Clean up build artifacts
clean:
	@echo "Cleaning up build artifacts..."
	@$(GOCLEAN)
	@rm -rf $(BUILD_DIR)

# Install golangci-lint if it's not already installed
install-linter:
	@echo "Checking for golangci-lint..."
	@which golangci-lint || (echo "Installing golangci-lint..."; $(GO) install github.com/golangci/golangci-lint/cmd/golangci-lint@latest)

# ====================================================================================
# Help Target
# ====================================================================================

# Display help information about the available targets
help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all              - Run lint and build for all platforms (default)."
	@echo "  build            - Build the application for all specified platforms."
	@echo "  build-local      - Build the application for the local OS/architecture."
	@echo "  lint             - Run the golangci-lint linter."
	@echo "  vet              - Run go vet."
	@echo "  fmt              - Format the Go source code."
	@echo "  tidy             - Tidy the go.mod and go.sum files."
	@echo "  clean            - Remove all build artifacts."
	@echo "  install-linter   - Install the golangci-lint tool."
	@echo "  help             - Show this help message."
	@echo ""
	@echo "Variables:"
	@echo "  PLATFORMS        - Override the list of platforms to build for."
	@echo "                   Example: make build PLATFORMS=\"linux/amd64 darwin/amd64\""

# Phony targets prevent conflicts with file names
.PHONY: all build build-local lint vet fmt tidy clean install-linter help

