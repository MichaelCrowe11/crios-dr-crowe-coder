# Backend Dockerfile for CriOS Dr. Crowe Coder - Railway Production
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies for RDKit and scientific computing
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libboost-all-dev \
    librdkit-dev \
    libpq-dev \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# Copy parent directory files needed for imports
COPY ../../src /app/src
COPY ../../crios /app/crios
COPY ../../pyproject.toml /app/pyproject.toml
COPY ../../requirements.txt /app/requirements-root.txt

# Install root dependencies
RUN pip install --no-cache-dir -r /app/requirements-root.txt

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /app/data /app/logs /app/uploads

# Set environment variables
ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1
ENV PORT=8000

# Expose port (Railway will set PORT env var)
EXPOSE ${PORT}

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:${PORT}/health || exit 1

# Run the application with production settings
CMD uvicorn app:app --host 0.0.0.0 --port ${PORT} --workers 4 --log-level info