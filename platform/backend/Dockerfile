# Backend Dockerfile for CriOS Dr. Crowe Coder
# NOTE: This Dockerfile must be built from the repository root
# Build context: repository root
# Docker build command: docker build -f platform/backend/Dockerfile -t crios-backend .

FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies for RDKit and compilation
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    librdkit-dev \
    libpq-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy root requirements first (for caching)
COPY requirements.txt /app/root-requirements.txt
COPY platform/backend/requirements.txt /app/requirements.txt

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /app/root-requirements.txt && \
    pip install --no-cache-dir -r /app/requirements.txt

# Copy source code from repository root
COPY src/ /app/src/
COPY crios/ /app/crios/

# Copy backend application code
COPY platform/backend/ /app/backend/

# Set up Python path
ENV PYTHONPATH=/app

# Create data directory
RUN mkdir -p /app/data

# Change to backend directory for execution
WORKDIR /app/backend

# Expose port (Railway sets PORT env var dynamically)
EXPOSE ${PORT:-8000}

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
  CMD curl -f http://localhost:${PORT:-8000}/health || exit 1

# Run the application
# Use PORT from environment if set (Railway requirement)
CMD uvicorn app:app --host 0.0.0.0 --port ${PORT:-8000}