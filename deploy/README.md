# CriOS Deployment Scripts

This directory contains deployment scripts and configurations for various platforms.

## Available Deployment Scripts

### Vercel Deployment (Frontend)

Deploy frontend to Vercel for global CDN and fast performance:

```bash
./deploy/vercel-deploy.sh
```

**Prerequisites:**
- Vercel CLI installed (`npm i -g vercel`)
- Vercel account
- For hybrid deployment: Backend deployed to Railway

**What it does:**
1. Checks Vercel CLI installation
2. Authenticates with Vercel
3. Guides you through deployment options:
   - **Hybrid**: Frontend on Vercel + Backend on Railway (recommended)
   - **Pure Vercel**: Everything on Vercel (limited RDKit support)
4. Configures environment variables
5. Deploys to production

**Deployment Options:**

**Option 1: Hybrid (Recommended)**
```bash
# 1. Deploy backend to Railway first
./deploy/railway-deploy.sh

# 2. Deploy frontend to Vercel
./deploy/vercel-deploy.sh
# Choose option 1 (Hybrid)
```

**Option 2: Pure Vercel**
```bash
./deploy/vercel-deploy.sh
# Choose option 2 (Pure Vercel)
# Note: Limited chemistry features, no RDKit
```

### Railway Deployment (Full Stack)

Deploy complete stack to Railway with RDKit support:

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
6. Deploys all services (backend, frontend, workers)

### Manual Deployment

See detailed deployment guides:
- [Vercel Deployment Guide](../VERCEL_DEPLOYMENT.md) - Frontend on Vercel
- [Railway Deployment Guide](../RAILWAY_DEPLOYMENT.md) - Full stack on Railway
- [General Deployment Guide](../DEPLOYMENT.md) - Docker and other options

## Files

- `vercel-deploy.sh` - Automated Vercel deployment script
- `railway-deploy.sh` - Automated Railway deployment script
- `production.yml` - Production Docker Compose configuration
- `deploy.sh` - General deployment helper script
- `README.md` - This file

## Quick Start

### Recommended: Hybrid Vercel + Railway

Best performance and full functionality:

```bash
# 1. Deploy backend to Railway
npm i -g @railway/cli
./deploy/railway-deploy.sh

# 2. Deploy frontend to Vercel
npm i -g vercel
./deploy/vercel-deploy.sh
```

**Result:**
- Frontend: https://your-project.vercel.app (global CDN)
- Backend API: https://your-backend.railway.app (full RDKit)
- **Cost:** ~$30-55/month

### Alternative: Railway Only

Full stack on one platform:

```bash
npm i -g @railway/cli
./deploy/railway-deploy.sh
```

**Result:**
- Complete stack on Railway
- **Cost:** ~$25-45/month

### Development: Docker Compose

Local development with all services:

```bash
# Start all services
docker-compose up --build

# Production mode
docker-compose -f docker-compose.yml -f deploy/production.yml up -d
```

## Environment Variables

Copy and configure environment variables:

```bash
# For Vercel
cp .env.vercel.example .env.vercel
# Edit .env.vercel with your values
# Or set in Vercel Dashboard: https://vercel.com/dashboard

# For Railway
cp .env.railway.example .env.railway
# Edit .env.railway with your values
# Or set in Railway Dashboard: https://railway.app/dashboard

# For Docker
cp .env.production .env
# Edit .env with your values
```

## Deployment Comparison

| Feature | Vercel (Hybrid) | Railway (Full) | Docker Compose |
|---------|----------------|----------------|----------------|
| **Setup Time** | 10 minutes | 10 minutes | 5 minutes |
| **RDKit Support** | ✅ (via Railway backend) | ✅ Full | ✅ Full |
| **Global CDN** | ✅ Vercel | ❌ | ❌ |
| **Auto SSL** | ✅ Both platforms | ✅ | Manual |
| **Scaling** | ✅ Automatic | ✅ Easy | Manual |
| **Cost** | $30-55/month | $25-45/month | Infrastructure cost |
| **Best For** | Production | Full stack simplicity | Development |

## Support

See detailed documentation:

### Deployment Guides
- [Vercel Deployment](../VERCEL_DEPLOYMENT.md) - Complete Vercel guide
- [Railway Deployment](../RAILWAY_DEPLOYMENT.md) - Complete Railway guide
- [Deployment Guide](../DEPLOYMENT.md) - General deployment

### Project Documentation
- [README](../README.md) - Project overview
- [Architecture](../docs/crowe-discovery-framework.md) - System architecture

### Platform Documentation
- [Vercel Docs](https://vercel.com/docs)
- [Railway Docs](https://docs.railway.app)

---

**CriOS Discovery Engine** • Deploy anywhere • Science before status. Discovery before profit.

© 2025 Crowe BioSystems
