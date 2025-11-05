# CriOS Deployment Scripts

This directory contains deployment scripts and configurations for various platforms.

## Available Deployment Scripts

### Railway Deployment

Deploy to Railway.app with automated setup:

```bash
./deploy/railway-deploy.sh
```

**Prerequisites:**
- Railway CLI installed (`npm i -g @railway/cli`)
- Railway account
- Anthropic API key

**What it does:**
1. Checks Railway CLI installation
2. Authenticates with Railway
3. Creates/links Railway project
4. Adds PostgreSQL and Redis services
5. Sets up environment variables
6. Deploys all services

### Manual Deployment

See detailed deployment guides:
- [Railway Deployment Guide](../RAILWAY_DEPLOYMENT.md)
- [General Deployment Guide](../DEPLOYMENT.md)

## Files

- `railway-deploy.sh` - Automated Railway deployment script
- `production.yml` - Production Docker Compose configuration
- `deploy.sh` - General deployment helper script

## Quick Start

### Railway (Recommended)

```bash
# Install Railway CLI
npm i -g @railway/cli

# Run deployment script
./deploy/railway-deploy.sh

# Or deploy manually
railway login
railway init
railway add --database postgresql
railway add --database redis
railway up
```

### Docker Compose

```bash
# Local development
docker-compose up --build

# Production
docker-compose -f docker-compose.yml -f deploy/production.yml up -d
```

## Environment Variables

Copy and configure environment variables:

```bash
# For Railway
cp .env.railway.example .env.railway
# Edit .env.railway with your values

# For Docker
cp .env.production .env
# Edit .env with your values
```

## Support

See main documentation:
- [Railway Deployment](../RAILWAY_DEPLOYMENT.md) - Complete Railway guide
- [Deployment Guide](../DEPLOYMENT.md) - General deployment
- [README](../README.md) - Project overview

---

© 2025 Crowe BioSystems
